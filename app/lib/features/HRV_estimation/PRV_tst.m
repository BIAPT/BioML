%% Load input
% load('PPG_ank_coh_str_coh')
% load('PPG_emp_coh_str_coh') 

load('PPG_ank_coh_rel_coh')
load('PPG_emp_coh_rel_coh')
%% Preprocess ankkoro data
% {
% Filtering
n_filt=2*sf_ank;                  % Filter order (transient will be 1 s)
b=fir1(n_filt,2/sf_ank*[0.6;6]);     % Generates fir BP filter
PPG_ank=filter(b,1,PPG_ank);
PPG_ank(isnan(PPG_ank))=[];
PPG_ank=PPG_ank(sf_ank+1:end);
b=fir1(2*sf_emp,2/sf_emp*[0.9;6]);
PPG_emp=filter(b,1,PPG_emp);
PPG_ank=PPG_ank(sf_emp+1:end);
% Resampling
t=(0:length(PPG_ank)-1)'*(1/sf_ank);
t_new=(0:length(PPG_ank)*(sf_emp/sf_ank)-1)'*(1/sf_emp);
PPG_ank_new=spline(t,PPG_ank,t_new);
PPG_ank=PPG_ank_new;
sf_ank=sf_emp;
%}

%% Use the algorithm
[ u_k_ank,du_k_ank,AR_Model_ank,PRV_features_ank ] = PRV_est( PPG_ank, sf_ank, 'pkfinder_mode',1, 'decimation_ratio',2,'ar_w',60,'ar_p',5,'ar_weight',0.9903);
[ u_k_emp,du_k_emp,AR_Model_emp,PRV_features_emp ] = PRV_est( PPG_emp, sf_emp, 'pkfinder_mode',1, 'decimation_ratio',2,'ar_w',60,'ar_p',5,'ar_weight',0.9903); 


figure; plot((0:length(PPG_emp)-1)*(1/sf_emp)*(1/60),PPG_emp); hold on; plot(u_k_emp*(1/sf_emp)*(1/60), PPG_emp(u_k_emp),'r*')
figure; plot((0:length(PPG_ank)-1)*(1/sf_ank)*(1/60),PPG_ank); hold on; plot(u_k_ank*(1/sf_ank)*(1/60), PPG_ank(u_k_ank),'r*')

t_emp = AR_Model_emp.t0 + (0:length(PRV_features_emp.PP)-1) * (1/sf_emp);
figure; plot(u_k_emp(2:end)*(1/sf_emp)*1/60, 1000*diff(u_k_emp)*(1/sf_emp), 'r*'); hold on;
plot(t_emp*(1/60), 1000*PRV_features_emp.PP)
legend('RR', 'First moment of IG distribution')
xlabel('time [s]')
ylabel('[ms]')

t_ank = AR_Model_ank.t0 + (0:length(PRV_features_ank.PP)-1) * (1/sf_ank);
figure; plot(u_k_ank(2:end)*(1/sf_ank)*1/60, 1000*diff(u_k_ank)*(1/sf_ank), 'r*'); hold on;
plot(t_ank*(1/60), 1000*PRV_features_ank.PP)
legend('RR', 'First moment of IG distribution')
xlabel('time [s]')
ylabel('[ms]')

figure; plot(t_emp*(1/60),PRV_features_emp.Mean_PP);hold on; plot(t_ank*(1/60),PRV_features_ank.Mean_PP);hold off;




%% Obtain the features
LF_emp=nldat(PRV_features_emp.LF','domainIncr',1/sf_emp);
HF_emp=nldat(PRV_features_emp.HF','domainIncr',1/sf_emp);
LFHF_emp=nldat(PRV_features_emp.LFHF_ratio','domainIncr',1/sf_emp);
Coherence_emp=nldat(PRV_features_emp.Coherence','domainIncr',1/sf_emp);

LF_ank=nldat(PRV_features_ank.LF','domainIncr',1/sf_ank);
HF_ank=nldat(PRV_features_ank.HF','domainIncr',1/sf_ank);
LFHF_ank=nldat(PRV_features_ank.LFHF_ratio','domainIncr',1/sf_ank);
Coherence_ank=nldat(PRV_features_ank.Coherence','domainIncr',1/sf_ank);


%% Decimate the features (st = 1s)


LF_emp=set(decimate(LF_emp,sf_emp*5),'domainStart',0);
HF_emp=set(decimate(HF_emp,sf_emp*5),'domainStart',0);
LFHF_emp=set(decimate(LFHF_emp,sf_emp*5),'domainStart',0);
Coherence_emp=set(decimate(Coherence_emp,sf_emp*5),'domainStart',0);

LF_ank=set(decimate(LF_ank,sf_ank*5),'domainStart',0);
HF_ank=set(decimate(HF_ank,sf_ank*5),'domainStart',0);
LFHF_ank=set(decimate(LFHF_ank,sf_ank*5),'domainStart',0);
Coherence_ank=set(decimate(Coherence_ank,sf_ank*5),'domainStart',0);


%% Plot Features


figure; plot(LFHF_emp); hold on; plot(LFHF_ank); hold off
figure; plot(Coherence_emp); hold on; plot(Coherence_ank); hold off
figure; plot(LF_emp); hold on; plot(LF_ank); hold off
figure; plot(HF_emp); hold on; plot(HF_ank); hold off