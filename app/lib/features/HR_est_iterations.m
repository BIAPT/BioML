function [ pos,HRest ] = HR_est_iterations( PPG,Acc )
% This function iterates the HR_est function when the analysis is done
% offline

% Initial parameters

t_window=8;         % Size of time window (in seconds)
t_step=1;       % Size of time step (in seconds)
sf=75;              % Sampling Frequency

% BandPass filter the whole signal
BW=[0.5 3.5];               % Hz - Bandwidth
n_filt=74;                  % Filter order (transient will be 1 s)
b=fir1(n_filt,2/sf*BW);     % Generates fir BP filter
a=1;
PPG=filter(b,a,PPG);        % Filters the signal
PPG=PPG(n_filt+1:end);      % Removes the transient

% It is important to notice that the filter is two-sided, but the signal
% was filtered as if it was one-sided. This means that the output is
% delayed by half the size of the filter.

% This calculates the number of windows that can be generated in the whole
% signal
N_windows=round((length(PPG)-t_window*sf)/(t_step*sf)+1);
% This initializes the vectors that will contain the estimated HR and their
% position in time

HRest=zeros(N_windows,1);
pos=zeros(N_windows,1);

if isempty(Acc)
    for i=1:N_windows
        % Now we estimate the position of every time window
        lower_index=1+(i-1)*t_step*sf;
        higher_index=min([((i-1)*t_step+t_window)*sf;length(PPG)]);
        ppg=PPG(lower_index:higher_index);
        acc=[];
        pos(i)=lower_index+t_window*sf/2+n_filt/2;
        if i<3
            [ ~,HRest(i) ] = HR_est( ppg, sf , acc, [],[]);
        else
            [ ~,HRest(i) ] = HR_est( ppg, sf , acc, HRest(i-1),HRest(i-2));
        end
    end
else
    for i=1:N_windows
        lower_index=1+(i-1)*t_step*sf;
        higher_index=min([((i-1)*t_step+t_window)*sf;length(PPG)]);
        ppg=PPG(lower_index:higher_index);
        acc=Acc(lower_index:higher_index,:);
        pos(i)=lower_index+t_window*sf/2+n_filt/2;
        if i<3
            [ ~,HRest(i) ] = HR_est( ppg, sf , acc ,[],[]);
        else
            [ ~,HRest(i) ] = HR_est( ppg, sf , acc , HRest(i-1),HRest(i-2));
        end
    end
end
end

