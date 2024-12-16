function [output] = tmcmc_par(varargin)
%% Transitional Markov Chain Monte Carlo sampler
% (Parallelized version)

% parse the information in the name/value pairs: 
pnames = {'nsamples','loglikelihood','priorpdf','priorrnd','burnin','lastburnin','beta'};
dflts =  {[],[],[],[],0,0,0.2}; % define default values
      
[nsamples,loglikelihood,priorpdf,prior_rnd,burnin,lastBurnin,beta] = ...
       internal.stats.parseArgs(pnames, dflts, varargin{:});

%% Number of cores
if ~isempty(gcp('nocreate'))
    pool = gcp;
    Ncores = pool.NumWorkers;
    fprintf('TMCMC is running on %d cores.\n', Ncores);
end

%% Obtain N samples from the prior pdf f(T)
j      = 0;                   % Initialise loop for the transitional likelihood
thetaj = prior_rnd(nsamples); % theta0 = N x D
pj     = 0;                   % p0 = 0 (initial tempering parameter)
Dimensions = size(thetaj, 2); % size of the vector theta

count = 1; % Counter
samps(:,:,count) = thetaj;
beta_j(count) = pj;
beta = 2.4./sqrt(Dimensions);
scale(count) = beta;

%% Initialization of matrices and vectors
thetaj1   = zeros(nsamples, Dimensions);

%% Main loop
while pj < 1    
    j = j+1;
    
    %% Calculate the tempering parameter p(j+1):
    log_fD_T_thetaj = zeros(nsamples, 1);
    parfor l = 1:nsamples
        log_fD_T_thetaj(l) = loglikelihood(thetaj(l,:));
    end
    if any(isinf(log_fD_T_thetaj))
        error('The prior distribution is too far from the true region');
    end
    pj1 = calculate_pj1(log_fD_T_thetaj, pj);
    fprintf('TMCMC: Iteration j = %2d, pj1 = %f\n', j, pj1);
    
    %% Compute the plausibility weight for each sample wrt f_{j+1}
    fprintf('Computing the weights ...\n');
    a       = (pj1-pj)*log_fD_T_thetaj;
    wj      = exp(a);
    wj_norm = wj./sum(wj);           % normalization of the weights
    
    %% Compute S(j) = E[w{j}] (eq 15)
    S(j) = mean(wj);
    
    %% Resampling and Metropolis-Hastings
    fprintf('Markov chains ...\n\n');
    idx = randsample(nsamples, nsamples, true, wj_norm);
    log_posterior = @(t) log(priorpdf(t)) + pj1*loglikelihood(t);

    % Weighted mean and covariance
    mu = zeros(1, Dimensions);
    parfor l = 1:nsamples
        mu = mu + wj_norm(l)*thetaj(l,:);
    end
    
    cov_gauss = zeros(Dimensions);
    parfor k = 1:nsamples
        tk_mu = thetaj(k,:) - mu;
        cov_gauss = cov_gauss + wj_norm(k)*(tk_mu'*tk_mu);
    end
    cov_gauss = beta^2 * cov_gauss + 1e-10*eye(size(cov_gauss));
    
    proppdf = @(x,y) prop_pdf(x, y, cov_gauss, priorpdf); %q(x,y) = q(x|y).
    proprnd = @(x)   prop_rnd(x, cov_gauss, priorpdf);   
    
    if pj1 == 1
        burnin = lastBurnin;
    end

    % Parallel Metropolis-Hastings
    parfor i = 1:nsamples
        [thetaj1(i,:), acceptance_rate(i)] = mhsample(thetaj(idx(i), :), 1, ...
            'logpdf',  log_posterior, ...
            'proppdf', proppdf, ...
            'proprnd', proprnd, ...
            'thin',    3,       ...
            'burnin',  burnin);
    end

    acceptance(count) = mean(acceptance_rate);
    
    %% Adjust acceptance rate and prepare for next iteration
    c_a = (acceptance(count) - ((0.21./Dimensions) + 0.23))./sqrt(j);
    beta = beta .* exp(c_a);
    
    count = count+1;
    scale(count) = beta;
    samps(:,:,count) = thetaj1;
    thetaj = thetaj1;
    pj     = pj1;
    beta_j(count) = pj;
end

%% Outputs
log_fD = sum(log(S(1:j)));
output.allsamples = samps;
output.samples = samps(:,:,end);
output.log_evidence = log_fD;
output.acceptance = acceptance;
output.beta = beta_j;
output.scale = scale(1:end-1);

return;

%% Calculate the tempering parameter p(j+1)
function pj1 = calculate_pj1(log_fD_T_thetaj, pj)
% find pj1 such that COV <= threshold, that is
%
%  std(wj)
% --------- <= threshold
%  mean(wj)
%
% here
% size(thetaj) = N x D,
% wj = fD_T(thetaj).^(pj1 - pj)
% e = pj1 - pj

threshold = 1; % 100% = threshold on the COV

% wj = @(e) fD_T_thetaj^e; % N x 1
% Note the following trick in order to calculate e:
% Take into account that e>=0
wj = @(e) exp(abs(e)*log_fD_T_thetaj); % N x 1
fmin = @(e) std(wj(e)) - threshold*mean(wj(e)) + realmin;
e = abs(fzero(fmin, 0)); % e is >= 0, and fmin is an even function
if isnan(e)
    error('There is an error finding e');
end

pj1 = min(1, pj + e);

return; % End

function proppdf = prop_pdf(x, mu, covmat, box)
% This is the Proposal PDF for the Markov Chain.

% Box function is the Prior PDF in the feasible region. 
% So if a point is out of bounds, this function will
% return 0.

proppdf = mvnpdf(x, mu, covmat).*box(x); %q(x,y) = q(x|y).

return;


function proprnd = prop_rnd(mu, covmat, box)
% Sampling from the proposal PDF for the Markov Chain.

while true
    proprnd = mvnrnd(mu, covmat, 1);
    if box(proprnd)
        % The box function is the Prior PDF in the feasible region.
        % If a point is out of bounds, this function will return 0 = false.
        break;
    end
end

return
