function [Thetap,Mu,Kappa,L,opt] = pplikel(EKGR, varargin)
% function [Thetap,Mu,Kappa,L,opt] = pplikel(EKGR, varargin)
%
%
% Copyright (C) Luca Citi and Riccardo Barbieri, 2010-2011.
% All Rights Reserved. See LICENSE.TXT for license details.
% {lciti,barbieri}@neurostat.mit.edu
% http://users.neurostat.mit.edu/barbieri/pphrv

% Default options
opt.delta = 0.003125;%+0.015; % time increment in updating parameters (in seconds)
opt.P = 11; % RR order
opt.hasTheta0 = 1; % wether or not the AR model has a theta0 constant to account for the average mu
opt.weight = 0.9903;%+3.02; % weighting factor
opt.W = 75; % window length for local likelihood estimate (in seconds)
opt.maximize_loglikel = @likel_invnorm_mex; % use loglikelihood of inverse gaussian

% PROCESS OPTIONS
% Checks initial options (varargin)
opt_names = fieldnames(opt);
J = 1;
while J <= length(varargin)
    key = varargin{J};
    keyN = find(strcmpi(key, opt_names));
    if isempty(keyN)
        warning('PPLIKEL:UnknownOption', 'Do not know option ''%s''', key);
    else
        opt.(opt_names{keyN}) = varargin{J+1};
    end
    J = J + 2;
end

% proxies for options
delta = opt.delta;          % step size for window
W = opt.W;                  % Window size
P = opt.P;                     % Linear order
% I ADDED THIS
M = 1;                          % NL quadratic order... Hopefully, this i=order is smaller than P
ordVolt=1;
% END
maximize_loglikel = opt.maximize_loglikel;

% Time vector: the time domain starts with the first beat
opt.t0 = EKGR(1);
EKGR = EKGR(:) - opt.t0; % times are relative to EKGR(1) which is set to 0
J = floor(EKGR(end) / delta) + 1; % Approximate size of the time vector

% Now find the last beat before the end of the time window, for that checks
% the first 5W beats (accounting for a maximum of 300bpm)
lastRi = find(EKGR(1:min(end,floor(5*W))) > W, 1) - 1; % index of last R event within W (check only first 5*W beats)
% Now, separate the sequence of observations until the end of the time
% window
observ_ev = EKGR(1:lastRi);

% Initialize parameters:
% The dynamic parameters Theta contain the AR parameters and the offset
% parameter. Each parameter varies with time for the total length J.
% I MODIFIED THIS
if ordVolt==1
    Thetap = NaN(P + opt.hasTheta0 , J);
elseif ordVolt==2
    Thetap = NaN(P + opt.hasTheta0 + (sum(1:M)), J);
elseif ordVolt==3
    Thetap = NaN(P + opt.hasTheta0 + (sum(1:M)+M+sum(1:M-1)*2+size(combnk(1:M,3),1)), J);
end
% END
% The mean parameter, which will depend on Theta
Mu = NaN(1, J);
% The kappa parameter or lambda in the IG representation (shape param)
Kappa = NaN(1, J);
% The number of steps at the end of the algorithm
steps = zeros(1, J);
% Local likelihood
L = NaN(1, J);
% STILL NOT SURE
meanRR = NaN(1, J);
opt.LogLikel = NaN(1, J);
% Instantaneous estimate of parameters
thetap = [];

% Now repeat iterations for all the time points after the first
% initialization window
for j = ceil(W / delta):J
    % variable 'time' is the time at the end of the window
    time = (j-1) * delta;    
    
    % If the first observed event is out of the analyzed window, it is
    % removed
    if ~isempty(observ_ev) && (observ_ev(1) < time - W)
        observ_ev(1) = []; % remove older event (there could be only one because we assume that in any delta interval there is at most one event)
        thetap = []; % force re-evaluation of starting point for thetap
    end
    % check if there is an event that will happen inside the new time step
    event = EKGR(lastRi + 1) <= time; % whether an event happened in ((j-1)*delta,j*delta]
    if event
        % If there is an event in that step, it will be added to the
        % sequence of observations
        lastRi = lastRi + 1;
        R = EKGR(lastRi);
        observ_ev(end+1,1) = R; % append current event
        thetap = []; % force re-evaluation of starting point for thetap
    end
    
    % Until here:
    % What the loop is doing is to first see if the first event in the
    % observed_events vector is still inside the analyzed time window. If
    % it is not, then that event is removed and the parameters vector theta
    % p is deleted to be re-estimated. 
    % Then, the loop checks if there's an event happening in the new step.
    % If there is, the event is added to the vector of observed events and
    % the parameter vector is emptied for re-estimation. 
    
    
    % if thetap is empty (i.e., observ_ev has changed) re-evaluate the variables that depend on observ_ev
    if isempty(thetap)
        % The position of the events will be taken from the vector of
        % observed events. We eliminate the first P+1 events, because we
        % need to remove the transient for the model. We can only estimate
        % the new events location after the transient is accounted for.
        uk = observ_ev(P+2:end);
        % the vector of rr intervals is computed for all the observations,
        % but it is only defined starting the second event (we need to have
        % intervals for this to be defined).
        rr = diff(observ_ev);
        % the vector wn is the 'output' of the AR model. It corresponds to
        % the next rr interval length, whihc is predicted by the AR model. 
        wn = rr(P+1:end);
        % matrix of observations for regression, and matrix with the last
        % observation for prediction
        xn = []; xt = [];
        % If our model has a zero order value, we add a column of ones. 
        if opt.hasTheta0
            xn = ones(length(wn),1);
            xt = 1;
        end
        % Now the regressors are constructed
        % xn is the linear regressor that contains the column of ones and
        % the lagged matrix for the AR model. The last event is not
        % considered because it belongs in the output matrix (future value
        % to be predicted)
        xn = [xn, toeplitz(rr(P:end-1), rr(P:-1:1))];
        
        % It seems to be the input vector for the last event. This means
        % that, this will be used, along with the estimated parameters, to
        % predict the rr interval for this instant in time. 
        xt = [xt; rr(end:-1:end-P+1)];
        
        % THIS IS WHAT I AM ADDING TO MAKE IT 2ND ORDER NONLINEAR
        if ordVolt==2
            xn = [xn, secord_reg(rr(P-M+1:end-1),M )];
            xt = [xt; (secord_reg( rr(end-M+1:end),M ))'];
        elseif ordVolt==3
            xn = [xn, secord_reg(rr(P-M+1:end-1),M ),third_reg( rr(P-M+1:end-1),M )];
            xt = [xt; (secord_reg( rr(end-M+1:end),M ))';(third_reg( rr(end-M+1:end),M ))'];
        end
        % END OF MY MODIFICATION
        
        % Computes the vector of weights of past events for the
        % loglikelihood function
        eta = weights(time, uk, opt.weight);
        % Now it performs uncensored loglikelihood maximization. 
        [thetap, k, steps(j)] = maximize_loglikel(xn, wn, eta); % the uncensored loglikelihood is a good starting point
        
        
        % Variables like wn, xn, xt are computed only inside this
        % conditional since they have information of observed events. This
        % does not change until the new window catches a new event or loses
        % an old event. Then, the parameter vector will be erased again and
        % the matrices and vectors would be updated. There is not need to
        % update these vectors until the vector of observed events changes.
        % Only the vector 'eta' changes, because this weighted vector
        % depends directly on the continues time value.
    else
        % If the thetap vector is not empty, it means that there exists an
        % estimate for the current series of observed events. Then we only
        % need to compute the new series of weights. 
        eta = weights(time, uk, opt.weight);
    end
    % Now the waiting time since the last event for the current window is
    % equal to the time at the end of the window minus the time of the last
    % observed event. 
    wt = time - observ_ev(end);
    % With the waiting time information, the algorithm computes the system
    % parameters again, having wt as the waiting time without event until
    % the end of the window (if wt=0 it means that the window ends exactly
    % at a event.
    % To maximize the waiting time likelihood, we need the previously
    % computed parameters like thetap and k
    [thetap, k, stepsj, L(j), loglikel] = maximize_loglikel(xn, wn, eta, thetap, k, xt, wt);
    steps(j) = steps(j) + stepsj;
    
    % The new rr(t) interval is predicted and added to the time series
    mu = thetap' * xt;
    Mu(1,j) = mu;
    % The new parameters vector is added to the TV matrix
    Thetap(:,j,1) = thetap;
    Kappa(j) = k;
    % meanRR is the weighted mean of the RR intervals inside the vector of observed events, given the weights
    % defined by 'eta'
    meanRR(j) = eta' * wn / sum(eta);
    opt.LogLikel(:,j) = sum(loglikel);

end

% Sepparates the zero order parameter from the total Thetap vector
if opt.hasTheta0
    opt.Theta0 = Thetap(1,:,1);
    Thetap(1,:,:) = [];
end

opt.steps = steps;
opt.meanRR = meanRR;

