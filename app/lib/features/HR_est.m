function [ BP , HR ] = HR_est( PPG , sf , Acc , HR_old , HR_old2 )
% Function HR_est uses PPG signals to estimate HR
%Needs to make a version in case there is not Acc data
% Date=180724;Hour=154702;Video1=3.01;Video2=3.1;Video3=3.05;Video4=3.07;

% Initial Parameters
e=15;           % Error allowed to not discard SSA results
nFFT=4096;      % Number of FFT points used to compute periodograms
delta_s=10;     % Maximum variation between time windows 
delta_sup=6;    % Upper rectifier amplitude
delta_inf=3;    % Lower rectifier amplitude (negative)
% The new estimated HR is a weighted average of the current estimate and
% the last 2 HR
alpha=0.9;      % Current estimate weight
beta=0.05;      % HR(k-1) weight
gamma=0.05;     % HR(k-2) weight


% RLS filtering is only performed if there is acceleration data, otherwise
% this step is skipped
% SSA is performed regardless the presence of acceleration data. If there
% is not acceleration data, it removes the non principal components of the
% series, or the ones that do not present harmonics. If acceleration data
% is given, the algorithm removes the components with similar structure to
% the acceleration data.
if isempty(Acc)
    BP=filtSSA(PPG,[],sf,HR_old,e,nFFT);
else
    BP=filtRLS(PPG,Acc)+filtSSA(PPG,Acc,sf,HR_old,e,nFFT);
end

% Finally, the position of the peaks in the periodogram is evaluated to
% provide the estimation of HR for an specific time window
[S_BP,F_BP]=periodogram(BP,[],nFFT,sf);
% If previous HR data is not provided, the algorithm assumes that the first
% two windows of the series are being analyzed
if isempty(HR_old)
    % The algorithm assumes that the series starts in a resting state,
    % therefore, initial spectral peaks far from 60 BPM are penalized
    % before looking for the maximum peak in the spectrum. 
    HRmin=round(40*nFFT/(60*sf))+1;
    HRmid=round(60*nFFT/(60*sf))+1;
    HRmax=round(180*nFFT/(60*sf))+1;
    G=[linspace(0.85,1,HRmid+1-HRmin) linspace(1,0.1,HRmax-HRmid)]';
    S_BP(HRmin:HRmax)=S_BP(HRmin:HRmax).*G;
    % Looks for the higher peak in the 40-180 BPM spectrum. If the highest
    % peak is far from 60 BPM even after penalization, that peak is taken
    % as the real HR for that window
    HR=F_BP(S_BP==max(findpeaks(S_BP(HRmin:HRmax))))*60;
else
    % When the preivous HRs are provided, the algorithm limits the search
    % to the vicinity of the previous HR, this way it avoids noise in other
    % frequencies that could bias the final estimation.
    HRmin=round((HR_old-delta_s)*nFFT/(60*sf))+1;
    HRmax=round((HR_old+delta_s)*nFFT/(60*sf))+1;
    MP=max(findpeaks(S_BP(HRmin:HRmax)));
    % If no peaks are found in the analyzed range, the algorithm copies the
    % previous HR
    if isempty(MP)
        HR=HR_old;
    else
        HR=F_BP(S_BP==MP)*60;
    end
    % Then, the weighted average is performed for the new estimated HR
    HR=alpha*HR+beta*HR_old+gamma*HR_old2;
end
% Finally, if the difference in the HR is higher or lower than the superior
% and inferior saturation points, then the new HR is set to match that
% saturation point
if (HR-HR_old)>delta_sup
    HR=HR_old+delta_sup;
elseif (HR_old-HR)>delta_inf
    HR=HR_old-delta_inf;
end
end

function ppg_new = filtRLS(PPG,Acc)
if isempty(Acc)
    disp('Error!: Acceleration data and PPG data should be the same length')
else
    % Filter Parameters:
    L=55;               % Length of the filter
    p0 = 100 * eye(L);  % Initial inverse covariance matrix
    lambda = 0.999;     % Forgetting factor
    % Now it checks the amount of channels provided by the acceleration
    % matrix. It is recommended to pass all three channels separated and
    % not a linear combination of them
    n_channels=size(Acc,2);
    % The matrix that will contain the filtered signals are created
    % The RLS filter will try to model the system that transforms the
    % acceleration into motion artifacts. Therefore, we do not want the
    % predicted output of the systems but the error, which corresponds to
    % the real PPG data.
    Filt_PPG=zeros(size(Acc,1),size(Acc,2)+1);
    Filt_Acc=zeros(size(Acc));
    Filt_PPG(:,1)=PPG;
    for i=1:n_channels
        % We need a different system for each channel, for this reason,
        % after the signal is filtered, we delete the filter a create a new
        % one for the next channel
        rls = dsp.RLSFilter(L,'ForgettingFactor',lambda,...
            'InitialInverseCovariance',p0);
        % [y,e]=rls(x,d); y is filter output, e is error, x is input signal and d
        % is the desired signal or contaminated signal
        [Filt_Acc(:,i),Filt_PPG(:,i+1)]=rls(Acc(:,i),Filt_PPG(:,i));
        clear rls
    end
    ppg_new=Filt_PPG(:,end); % The last column of the filtered signal is the function output
end
end
function ppg_new = filtSSA(PPG,Acc,sf,HR_old,e,nFFT)
if isempty(Acc)
    disp('No acceleration data, algorithm will not account for motion artifacts')
    % If there is not acceleration data, the algorithm will only remove non
    % principal components and keep the ones that present harmonics
    L=round(length(PPG)*0.4);   % L is close to half the length of the signal
    K=length(PPG)-L+1;          % K is calculated
    ppg_1=PPG(1:L);             % signal is 
    ppg_2=PPG(L:end);           % decomposed to
    PPG_H=hankel(ppg_1,ppg_2);  % form Hankel Matrix
    
    d=min(L,K);                 % Number of rank-one matrices
    PPG_set=zeros([size(PPG_H) d]);   
    [U,S,V]=svd(PPG_H);         % Performs SVD on PPG
    s=diag(S);
    for i=1:d
        PPG_set(:,:,i)=s(i)*U(:,i)*(V(:,i)');
    end
    % It looks for components with almost the same singular value, which
    % correpsond to a group of harmonics. 
    dds=diff(s);                            % Diff of singular values
    dds(2:end)=dds(2:end)/sum(s(2:end));    % normalize singular values but the first one
    [pks,locs]=findpeaks(dds,'MinPeakProminence',0.0005);
    g=length(pks); % the number of final groups is equal to the number of clusters
    % Now we reconstruct the signal with the provided components
    PPG_I_set=zeros(size(PPG_set,1),size(PPG_set,2));
    for j=1:g
        PPG_I_set=sum(PPG_set(:,:,locs(j):locs(j)+1),3)+PPG_I_set;
    end

    % 4. Reconstruction of time series
    ppg_recon=zeros(length(PPG),1);
    PPG_I_set=flipud(PPG_I_set);
    for j=1:L
        ppg_recon(j)=mean(diag(PPG_I_set,j-L));
    end
    for j=1:K-1
        ppg_recon(j+L)=mean(diag(PPG_I_set,j));
    end
    ppg_new=ppg_recon;    
    
else
    
    L=round(length(PPG)*0.4);   % L is close to half the length of the signal
    K=length(PPG)-L+1;        % K is calculated
    ppg_1=PPG(1:L);             % signal is 
    ppg_2=PPG(L:end);           % decomposed to
    PPG_H=hankel(ppg_1,ppg_2);      % form Hankel Matrix

    d=min(L,K);                 % Number of rank-one matrices
    PPG_set=zeros([size(PPG_H) d]);   
    [U,S,V]=svd(PPG_H);             % SVD on X

    s=diag(S);
    for i=1:d
        PPG_set(:,:,i)=s(i)*U(:,i)*(V(:,i)');
    end

    % 3. Grouping
    % Analysis of singular values: Usually signals can be clustered if they
    % have close singular values. Therefore, the difference of adjacent singular values
    % (which is always negative because the singular values are monotonically
    % decreasing) will have peaks for close adjacent values. The first singular
    % value corresponds to the trend of the signal and might be very large with
    % respect to the other values, for this reason, we do not take it into
    % account when calculating the %VAF of each singular value to the complete
    % signal. Also, to ensure that there is a clear peak (similar values), we
    % restrict the prominence of the peaks to be at least .05%, otherwise the
    % peaks are not considered
    dds=diff(s);                            % Diff of singular values
    dds(2:end)=dds(2:end)/sum(s(2:end));    % normalize singular values but the first one
    [pks,locs]=findpeaks(dds,'MinPeakProminence',0.0005);
    g=length(pks); % the number of final groups is equal to the number of clusters
    % Now we separate the groups:
    PPG_I_set=zeros(size(PPG_set,1),size(PPG_set,2),g);
    for j=1:g
        PPG_I_set(:,:,j)=sum(PPG_set(:,:,locs(j):locs(j)+1),3);
    end

    % 4. Reconstruction of time series
    ppg_recon=zeros([length(PPG) g]);
    PPG_I_set=flipud(PPG_I_set);
    % We need to find a way to improve the reconstruction in terms of
    % computating time
    % This reonstructs the time series for each identified group 
    for i=1:g
        for j=1:L
            ppg_recon(j,i)=mean(diag(PPG_I_set(:,:,i),j-L));
        end
        for j=1:K-1
            ppg_recon(j+L,i)=mean(diag(PPG_I_set(:,:,i),j));
        end
    end
    % //Improve this loop
    % If there is previous HR data, 
    if ~isempty(HR_old)
        % The algorithm looks for the dominant peaks in the acceleration
        % spectrum (at least 50% of the maximum peak)
        Spect_ppg=periodogram(ppg_recon,[],nFFT,sf);
        S_Acc=periodogram(Acc,[],nFFT,sf);
        [~,locs_x]=findpeaks(S_Acc(:,1),'MinPeakHeight',0.5*max(S_Acc(:,1)));
        [~,locs_y]=findpeaks(S_Acc(:,2),'MinPeakHeight',0.5*max(S_Acc(:,2)));
        [~,locs_z]=findpeaks(S_Acc(:,3),'MinPeakHeight',0.5*max(S_Acc(:,3)));
        HR_old_pos=round(HR_old*nFFT/(60*sf));
        locs_acc=sort([locs_x;locs_y;locs_z]);
        ppg_new=zeros(size(ppg_recon,1),1);
        % If the location of the dominant peaks in the acceleration
        % spectrum matches peaks in one of the principal components of the
        % PPG signal, that component is removed. However, if that component
        % also has a peak around HR_old, then that component is kept.
        for i=1:g
            [~,locs]=findpeaks(Spect_ppg(:,i),'MinPeakHeight',0.5*max(Spect_ppg(:,i)));
            
            if(sum(abs(locs-HR_old_pos)<6))
                ppg_new=ppg_new+ppg_recon(:,i);
            elseif ~sum(sum(abs(locs_acc-locs')<4))
                ppg_new=ppg_new+ppg_recon(:,i);
            end
        end
        % Finally, if the new PPG signal does not have peaks around HR_old
        % (allowing an error band) the signal is discarded
        [S2,F]=periodogram(ppg_new,[],nFFT,sf);
        if abs(F(S2==max(S2))*60-HR_old)>e
            ppg_new=0;
        end
    else
        % If there is not HR_old data, then the signal is the sum of all
        % the groups
        ppg_new=sum(ppg_recon,2);
    end
end
end