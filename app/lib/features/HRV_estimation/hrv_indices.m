function [powLF, powHF, bal, powCohe, warn, powVLF, powTot] = hrv_indices(Thetap, Var, fsamp, dwsample)
% function [powLF, powHF, bal, warn] = hrv_indices(Thetap, Var, fsamp, dwsample)
%
% Evaluate Heart Rate Variability indices from the outputs of
% the point-process likelihood modeling (pplikel.m)
%
%
% Copyright (C) Luca Citi and Riccardo Barbieri, 2010-2011.
% All Rights Reserved. See LICENSE.TXT for license details.
% {lciti,barbieri}@neurostat.mit.edu
% http://users.neurostat.mit.edu/barbieri/pphrv

[P,J] = size(Thetap);

% initialize output
% powLF and powHF are NaN where they cannot be evaluated (e.g. at the beginning)
% powLF = NaN(1,J);
% powHF = NaN(1,J);
% powVLF = NaN(1,J);
% powTot = NaN(1,J);
% powCohe = NaN(1,J);
powLF = zeros(1,J);
powHF = ones(1,J);
powVLF = zeros(1,J);
powTot = zeros(1,J);
powCohe = zeros(1,J);


% warn==1 is OK
% warn~=1 warning:
%   mod(warn, 2) is the scale factor used to shrink the poles in case of instability (slightly less than 1 is fine)
%   bitand(floor(warn), 2) if powLF was negative (increase AR order?)
%   bitand(floor(warn), 4) if powHF was negative (increase AR order?)
warn = NaN(1,J);

for i = 1:J
    if isnan(Thetap(1,i))||(Var(i)<((1/9)*median(Var)))
        continue;
    end
    [tot,comps,compsp,f,pole_freq,pole_pow,pole_res,poles,mod_scale] = spectral(Thetap(:,i), Var(i), fsamp(i), 256, 0);
    warn(i) = mod_scale;
    pf = abs(pole_freq);
%     powLF(i) = sum(pole_pow((pf>0.04) & (pf<0.15)));
    powLF(i) = trapz(f((f>0.04) & (f<0.15)),tot((f>0.04) & (f<0.15)));
    if powLF(i) <= 0
        powLF(i) = powLF(i-1);
        warn(i) = warn(i) + 2;
    end
%     powHF(i) = sum(pole_pow((pf>0.15) & (pf<0.45)));
    powHF(i) = trapz(f((f>0.15) & (f<0.45)),tot((f>0.15) & (f<0.45)));
    if powHF(i) <= 0
        powHF(i) = powHF(i-1);
        warn(i) = warn(i) + 4;
    end
%     powVLF(i) = sum(pole_pow(pf<0.04));
%     powTot(i) = sum(pole_pow);
    powVLF(i) = trapz(f(f<0.04),tot(f<0.04));
    powTot(i) = trapz(f,tot);
    if ~isnan(Thetap(:,i))
        mmmm=0;
    end
    uf=f((f>0.04) & (f<0.26));
    [peaks,~]=findpeaks(tot((f>0.04) & (f<0.26)));
    if isempty(peaks)
        uf=0.1;
    else
        uf=f(tot==max(peaks));
    end
%     uf=uf(tot((f>0.04) & (f<0.26))==max(tot((f>0.04) & (f<0.26))));
    uf1=uf-0.015;
    uf2=uf+0.015;
    t1=sum(tot(f<uf1));
    t2=sum(tot((f>uf1) & (f<uf2)));
    t3=sum(tot(f>uf2));
    powCohe(i)=t2/(t1+t3);
end

b = ceil(J/10);
if b > 1
    b = hamming(min(21, b));
    b = b/sum(b);
    powLF = fliplr(filter(b,1,fliplr(filter(b,1,powLF))));
    powHF = fliplr(filter(b,1,fliplr(filter(b,1,powHF))));
    powVLF = fliplr(filter(b,1,fliplr(filter(b,1,powVLF))));
    powTot = fliplr(filter(b,1,fliplr(filter(b,1,powTot))));
end

if nargin > 3 && dwsample ~= 1
    powLF = decimate(powLF, dwsample, 'fir');
    powHF = decimate(powHF, dwsample, 'fir');
    powVLF = decimate(powVLF, dwsample, 'fir');
    powTot = decimate(powTot, dwsample, 'fir');
end

bal = powLF ./ powHF;
