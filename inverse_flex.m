clc; clear; close all;

num_layers = 6;  
sigma_res = 0.01;
sigma_phase = 0.01;
NumFreq = 50;    
frequency = logspace(-1, 2, NumFreq); 

res_true = [500, 1000, 100, 500, 1000, 200]; 
thickness_true = [210, 1624, 1346, 1435, 1800];  

[AppRes_obs, Phase_obs] = deal(zeros(NumFreq, 1));
for j = 1:NumFreq
    [tmp1, tmp2, ~] = MTmodeling1D(res_true, thickness_true, frequency(j));
    AppRes_obs(j) = tmp1 * (1 + randn * sigma_res);  
    Phase_obs(j) = tmp2 * (1 + randn * sigma_phase); 
end

disp('AppRes_obs:');
disp(AppRes_obs);
disp('Phase_obs:');
disp(Phase_obs);
disp('Frequency:');
disp(frequency);

sigma_res1 = 5;
sigma_phase1 = 5;

priorpdf = @(theta) double(... 
    (theta(1) > 450 && theta(1) < 550) && ...
    (theta(2) > 950 && theta(2) < 1050) && ...
    (theta(3) > 50 && theta(3) < 150) && ...
    (theta(4) > 450 && theta(4) < 550) && ...
    (theta(5) > 950 && theta(5) < 1050) && ...
    (theta(6) > 150 && theta(6) < 250) && ...
    (theta(7) > 200 && theta(7) < 300) && ...
    (theta(8) > 1600 && theta(8) < 1700) && ...
    (theta(9) > 1300 && theta(9) < 1400) && ...
    (theta(10) > 1400 && theta(10) < 1500) && ...
    (theta(11) > 1800 && theta(11) < 1900));

priorrnd = @(n) [...
    rand(n, 1) * (100) + 450, ...  % res_1
    rand(n, 1) * (100) + 950, ...  % res_2
    rand(n, 1) * (100) + 50,  ...  % res_3
    rand(n, 1) * (100) + 450, ...  % res_4
    rand(n, 1) * (100) + 950, ...  % res_5
    rand(n, 1) * (100) + 150, ...  % res_6
    rand(n, 1) * (100) + 200,  ... % thickness_1
    rand(n, 1) * (100) + 1600, ... % thickness_2
    rand(n, 1) * (100) + 1300, ... % thickness_3
    rand(n, 1) * (100) + 1400, ... % thickness_4
    rand(n, 1) * (100) + 1800];    % thickness_5


loglikelihood_handle = @(theta) loglikelihood(theta, AppRes_obs, Phase_obs, frequency, sigma_res1, sigma_phase1, num_layers, NumFreq);

fprintf('Running TMCMC...\n');
tic;
output_tmcmc = tmcmc_par('nsamples', 50000, ...
                         'loglikelihood', loglikelihood_handle, ...
                         'priorpdf', priorpdf, ...
                         'priorrnd', priorrnd, ...
                         'burnin', 10, ...
                         'lastburnin', 50);
time_tmcmc = toc;

samples_tmcmc = output_tmcmc.samples;
save('output_tmcmc5', 'output_tmcmc');
save('time_tmcmc5', 'time_tmcmc');

fprintf('Running Langevin TMCMC...\n');
tic;
output_langevin = ltmcmc_par('nsamples', 5000, ...
                             'loglikelihood', loglikelihood_handle, ...
                             'priorpdf', priorpdf, ...
                             'priorrnd', priorrnd, ...
                             'burnin', 20, ...
                             'epsilon', 0.5);
time_ltmcmc = toc;

samples_langevin = output_langevin.samples;
save('output_langevin5', 'output_langevin');
save('time_ltmcmc5', 'time_ltmcmc');


function loglike = loglikelihood(theta, AppRes_obs, Phase_obs, frequency, sigma_res1, sigma_phase1, num_layers, NumFreq)
    misfit_res = zeros(NumFreq, 1);
    misfit_phase = zeros(NumFreq, 1);

    for j = 1:NumFreq
        [tmp_res, tmp_phase, ~] = MTmodeling1D(theta(1:num_layers), theta(num_layers+1:end), frequency(j));
        misfit_res(j) = (AppRes_obs(j) - tmp_res) / sigma_res1;
        misfit_phase(j) = (Phase_obs(j) - tmp_phase) / sigma_phase1;
    end

    loglike = -0.5 * (norm(misfit_res)^2 + norm(misfit_phase)^2);
end
