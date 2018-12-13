function [ u_k,du_k,AR_Model,PRV_features ] = PRV_est( PPG, sf, varargin )
% This function receives the PPG sequence as an input and returns the
% modelled PP sequence and PRV features
% Inputs:
%           - PPG: The PPG signal as a column vector
%           - sf: The sampling frequency
% Outputs:
%           - PP sequence
% ----------------------------------------------------------------------------------
options={{'decimation_ratio' 10 'ratio to decimate data bt'} ...
         {'maxhr' 180 'Maximum HR accepted'} ...
         {'shift' 1/3 'shift of the ddt2 peak locations to search for PPG peaks'} ...
         {'bw_ppg' [0.5 5] 'Bandwidth of the PPG signal for filtering'} ...
         {'pkfinder_mode' 1 'Mode 1: derivative, Method 2: Hilbert transform'} ...
         {'ensemble_filter' 'yes' '{yes, no} using ensemble filter'} ...
         {'wl_ens' 0.25 'size of the window for ensemble filter'}...
         {'th_corr' 0.45 'Correlation threshold for ensemble filter'}...
         {'median_filter' 'yes' '{yes, no} using median filter'} ...
         {'th_med' 5 'Threshold for median filtering'}...
         {'missing_beat_correction' 'yes' '{yes, no} correcting missing beats'} ...
         {'ar_p' 5 'Linear AR order'}...
         {'ar_w' 90 'Window size for AR modelling'}...
         {'ar_weight' 0.98 'weight of past values for AR model'}...
         
     };
 if arg_parse(options,varargin)
     return
 end
% -----------------------------------------------------------------------------------------------------
% Initialization
thr_rrlength=ceil(sf*60/maxhr);
datalength=length(PPG)/sf;
PPG=nldat(PPG,'domainIncr',1/sf);

% Initial peak search
disp('Algorithm is searching for the peaks in the signal')
tic
switch pkfinder_mode
    case 1
        u_k=pkfinder1(PPG,thr_rrlength,decimation_ratio,shift);
    case 2
        u_k=pkfinder2(PPG,sf,thr_rrlength,bw_ppg,shift);
end
t1=toc;
disp(['Elapsed time is ' num2str(t1) ' seconds.'])

% Filtering the PP interval series
disp('Algorithm is correcting the PP interval time series')
tic
switch ensemble_filter
    case 'yes'
        errors_corr=ens_filt(u_k,PPG,wl_ens,th_corr);
    case 'no'
        errors_corr=zeros(size(u_k));
end
switch median_filter
    case 'yes'
        errors_med=med_filt(u_k,th_med);
    case 'no'
        errors_med=zeros(size(u_k));
end
u_k(errors_med|errors_corr)=[];
switch missing_beat_correction
    case 'yes'
        [u_k,du_k]=missbeat_correction(u_k,thr_rrlength);
    case 'no'
        du_k=diff(u_k);
end
t2=toc;
disp(['Elapsed time is ' num2str(t2) ' seconds.'])

% Fitting Linear Autoregressive model of heart rate dynamics
% receives u_k vector in seconds
disp('Algorithm is fitting AR model to the sequence of heartbeat intervals')
tic
[Thetap,Mu,Kappa,L,opt] = pplikel(u_k*(1/sf),'delta',(1/sf),'P',ar_p,'W',ar_w,'weight',ar_weight);
Var = opt.meanRR.^3 ./ Kappa; % variance of an inverse Gaussian
Var = 1e6 * Var; % from [s^2] to [ms^2]
[powLF, powHF, bal,Crat] = hrv_indices(Thetap, Var, 1./opt.meanRR);
AR_Model.Theta=Thetap;
AR_Model.Theta0=Thetap;
AR_Model.Kappa=Kappa;
AR_Model.LogLikel=opt.LogLikel;
AR_Model.L=L;

PRV_features.PP=Mu;
PRV_features.Var=Var;
PRV_features.Mean_PP=opt.meanRR;
PRV_features.LF=powLF;
PRV_features.HF=powHF;
PRV_features.Coherence=Crat;
PRV_features.LFHF_ratio=bal;
t3=toc;
disp(['Elapsed time is ' num2str(t3) ' seconds.'])

disp(['The algorithm took ' num2str(t1+t2+t3) ' seconds to process ' num2str(datalength) ' seconds of data'])
end

function u_k=pkfinder1(PPG,thr_rrlength,decimation_ratio,shift)
% This pkfinder uses the derivatives of the PPG signal to find the
% beginning of each heartbeat cycle

% -------------------------------------------------------------------------------------------
% WARNING: When using this method, it is important to set the decimation
% ratio accordingly with the sampling frequency. Given the nature of the
% PPG signal, it is preferred that the decimated sampling frequency is at
% least 30 Hz. Therefore, when using the empatica (sf=64) the decimation
% ratio cannot be more than 2.
% -------------------------------------------------------------------------------------------

% Decimates the derivative of the PPG signal to reduce noise
DDT=decimate(ddt(PPG),decimation_ratio);
% Now it takes the derivative of the half-wave rectified velocity
DDT=double(ddt(DDT+abs(DDT)));
% Look for the peaks of the new signal, such that they are large enough to
% correspond to the beginning of a cycle.
[~,locs]=findpeaks(DDT,'MinPeakDistance',thr_rrlength/decimation_ratio,'MinPeakProminence',std(DDT));
% Multiply the locations by the decimation ratio to obtain location in the
% original signal
locs=locs*decimation_ratio;
% Now we want to find the real negative peaks in the PPG signal, using the
% derivative peaks as a reference of the length of each cycle. The second 
% derivative peaks are located close to the points we are searching for.
% Thus, we need to shift the obtained locations to make them fall in the
% middle of a cycle before searching for the real peaks.
L_shift=round((shift)*min(diff(locs)));
locs=locs+L_shift;
% Now the signal is filtered with a short average filter to avoid tagging
% peaks resulting from noise
PPG=smo(PPG,3);

% Use the found local minima of ddt2 to find the local minima in the original
% signal and prevent time shift from filter
for i=0:length(locs)-1
    % For the first location, look for the peak that is located from the
    % first value of the vector until there.
    if i==0
        tempPPG=double(PPG(1:locs(1)));
        [~,loctemp,~,prom]=findpeaks(-tempPPG,'MinPeakDistance',min([thr_rrlength;length(tempPPG)-2]));
        locs2=loctemp(prom==max(prom));
        % For the other locations, look for the peaks in between each
        % location
    else
        tempPPG=double(PPG(locs(i):min([locs(i+1);length(PPG)])));
        [~,loctemp,~,prom]=findpeaks(-tempPPG,'MinPeakDistance',thr_rrlength);
        Ltemp=loctemp(prom==max(prom));
        % Even though it shouldn't be the case, if we do not find a peak in
        % the given interval, it might mean that the one of the extremes of
        % the search interval has fallen exactly in the location of a peak.
        % Therefore, we expand our search interval in one point on each
        % extreme
        if isempty(prom)
            tempPPG=double(PPG(locs(i)-1:locs(i+1)+1));
            [~,loctemp,~,prom]=findpeaks(-tempPPG,'MinPeakDistance',thr_rrlength);
            Ltemp=loctemp(prom==max(prom));
            % If we find a peak, we copy its location in the locs vector
            if ~isempty(prom)
                locs(i)=Ltemp+locs(i)-2;
            % Otherwise, we extrapolate the location using the previous two
            % locations.
            else
                locs(i)=locs(i-1)+(locs(i-1)-locs(i-2));
            end
        % If we found a peak initially, then the location is saved
        else
            locs(i)=Ltemp+locs(i)-1;
        end
    end
end
% The vector of peak locations is constructed 
u_k=[locs2;locs(1:end-1)];

end
function u_k=pkfinder2(PPG,sf,thr_rrlength,BW_PPG,shift)
% Some signals (like the ones from the empatica band) will have a very
% clear local minima and little to not baseline drifting. In these cases,
% we can use the hilbert transform to find the position of the local minima

% Filter the signal to remove noise
n_filt=2*sf;                       % Filter order (transient will be 1 s)
b=fir1(n_filt,2/sf*BW_PPG);     % Generates fir BP filter
a=1;

% Filter the PPG signal
PPG2=filter(PPG,b,a);
PPG2=double(PPG2);
% Remove the transient
PPG2(1:2*sf)=[];

% Use the Hilbert transform to find local minima in the filtered signal
% The Hilbert transform estimates the intantaneous phase of a signal. The
% phase slips when there is a local minima, these slips can be used to find
% the location of the minimum points in the PPG signal.
Hbvp = hilbert(PPG2);
pha=angle(Hbvp);
% Since the peak slip is from pi to -pi, then, we count the peaks that are
% at least larger than 5 rads.
[~,locs]=findpeaks(diff(-pha), 'MinPeakDistance',thr_rrlength,'Threshold',5);
% The location of the peaks found here correspond to the filtered signal,
% which has a delay due to filtering. Here we account for that delay:
locs=locs+sf;
L_shift=round((shift)*min(diff(locs)));
locs=locs+L_shift;
% Use the found local minima of Hilbert phase to find the local minima in the original
% signal and avoid time shift from filtering
for i=0:length(locs)-1
    % For the first location, look for the peak that is located from the
    % first value of the vector until there.
    if i==0
        tempPPG=double(PPG(sf+1:locs(1)));
        [~,loctemp,~,prom]=findpeaks(-tempPPG,'MinPeakDistance',min([thr_rrlength;length(tempPPG)-2]));
        locs2=loctemp(prom==max(prom))+sf;
        % For the other locations, look for the peaks in between each
        % location
    else
        tempPPG=double(PPG(locs(i):min([locs(i+1);length(PPG)])));
        [~,loctemp,~,prom]=findpeaks(-tempPPG,'MinPeakDistance',thr_rrlength);
        Ltemp=loctemp(prom==max(prom));
        % Even though it shouldn't be the case, if we do not find a peak in
        % the given interval, it might mean that the one of the extremes of
        % the search interval has fallen exactly in the location of a peak.
        % Therefore, we expand our search interval in one point on each
        % extreme
        if isempty(prom)
            tempPPG=double(PPG(locs(i)-1:locs(i+1)+1));
            [~,loctemp,~,prom]=findpeaks(-tempPPG,'MinPeakDistance',thr_rrlength);
            Ltemp=loctemp(prom==max(prom));
            % If we find a peak, we copy its location in the locs vector
            if ~isempty(prom)
                locs(i)=Ltemp+locs(i)-2;
            % Otherwise, we extrapolate the location using the previous two
            % locations.
            else
                locs(i)=locs(i-1)+(locs(i-1)-locs(i-2));
            end
        % If we found a peak initially, then the location is saved
        else
            locs(i)=Ltemp+locs(i)-1;
        end
    end
end
% The vector of peak locations is constructed 
u_k=[locs2;locs(1:end-1)];
end
function errors_corr=ens_filt(u_k,PPG,wl_ens,th_corr)
% defines the size of the window length for the ensemble filter
wl=round(wl_ens*median(diff(u_k))); % wave length to compare = 2*wl+1
ppg=double(PPG);
% separates the waves around each peak
init_idx=u_k(1)-wl;
end_idx=u_k(end)+wl;
if (init_idx<1)||(end_idx>length(ppg))
    wave_set=ppg((-wl:wl)'+u_k(2:end-1)');
    if (init_idx<1)
        wave_set=[[zeros(1-init_idx,1);ppg((1:wl+u_k(1)))] wave_set];
    else
        wave_set=[ppg((-wl:wl)'+u_k(1)) wave_set];
    end
    if (end_idx>length(ppg))
        wave_set=[wave_set [ppg(u_k(end)-wl:end);zeros(end_idx-length(ppg),1)]];
    else
        wave_set=[wave_set ppg((-wl:wl)'+u_k(end))];
    end
else
    wave_set=ppg((-wl:wl)'+u_k');
end
% Finds the average wave
Av_wave=mean(wave_set,2);
% Normalizes the waves and the average wave for correlation
Av_wave=(Av_wave-mean(Av_wave))/std(Av_wave);
wave_set=(wave_set-mean(wave_set,1))./std(wave_set,0,1);
% Obtains the static correlation coefficient
CorrCoef=(1/(2*wl))*(wave_set'*Av_wave);
% Generates vector with erroneous peaks
errors_corr=CorrCoef<th_corr;
end
function errors_med=med_filt(u_k,th_med)
% This filter searches interval lengths that are too far away from the
% median length and marks them to be later removed
d=diff(u_k);
sM=median(d);
med=median(abs(d-sM));
D=abs(d-sM)./(1.483*med);
errors_med=[0;D>th_med];
end
function [u_k,du_k]=missbeat_correction(u_k,thr_rrlength)
% This function identifies missing beats based on their neighbours and
% fills the blank space with beats of random size depending on the
% information of the previous beats
du_k=diff(u_k);
i=1;
% I use while instead of for, because the size of the vector u_k might
% change if we need to add beats.
while i<length(du_k)
    % If the current interval is too large compared to the previous
    % intervals, we might be dealing with missing beats
    if ((1.5*du_k(i))>(2*(mean(du_k(max([i-4;1]):i))))) 
        % We estimate the number of beats that should be inside that
        % interval by dividing the size of the current interval by the size
        % of the mean interval of the last 10 beats
        Nbeats=round(du_k(i)/mean(du_k(max([i-10;1]):i-1))-0.1);
        % We use the mean and std of the last 10 beats to make our new
        % interval sequence
        M=mean(du_k(max([i-10;1]):i-1));
        S=std(du_k(max([i-10;1]):i-1));
        % The new interval lengths are drawn from a uniform distribution in
        % the interval given by the mean interval length +/- the std.
        du_k_new=randi([round(M-S); round(M+S)],Nbeats,1);
        % We normalize the obtained intervals such that the sum of them is
        % equal to the long initial interval. This way, the position of the
        % peaks in the whole sequence will not be affected.
        du_k_new=round(du_k_new*(du_k(i))/sum(du_k_new));
        % We now need to insert the new sequence of peak locations inside
        % the initial vector of locations
        Nu_k=cumsum([u_k(i);du_k_new(1:end-1)]);
        u_k=[u_k(1:i-1);Nu_k;u_k(i+1:end)];
        % And the new sequence of interval lengths is defined by:
        du_k=diff(u_k);
    % If there was no missing beats, we check if the beat we are dealing
    % with is too short compared to the maximum HR possible.
    elseif (du_k(i))<thr_rrlength
        u_k(i+1)=[];
        du_k=diff(u_k);
        i=i-1;
    end
    i=i+1;
end
end