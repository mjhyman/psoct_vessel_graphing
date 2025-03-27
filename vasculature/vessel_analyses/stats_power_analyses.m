%% Power analysis to determine sample size for LME model
% The purpose is to determine the minimum sample size to achieve an alpha
% of 0.05 (p value) and a beta of 0.80 (power) for the linear mixed effects
% model (LMM). The literature states that the optimal method of calculating
% the minimum sample size is through simulations.
%
% This script examines the vascular metrics for AD vs. control


%% Clear workspace & add top-level directory
clear; clc; close all;
% Start in current directory
mydir  = pwd;
% Find indices of slashes separating directories
if ispc
    idcs = strfind(mydir,'\');
elseif isunix
    idcs = strfind(mydir,'/');
end
% Truncate path to reach top-level directory (psoct_vessel_graphing)
topdir = mydir(1:idcs(end-1));
addpath(genpath(topdir));

%% Initialize the subject ID lists
% IDs of each subject
subids = {'AD_10382', 'AD_20832', 'AD_20969', 'AD_21354', 'AD_21424',...
         'CTE_6489', 'CTE_6912', 'CTE_7019', 'CTE_7126',...
         'NC_8095','NC_6839','NC_6974','NC_8653','NC_21499'};
ad_subs = {'AD_10382', 'AD_20832', 'AD_20969', 'AD_21354', 'AD_21424'};
cte_subs = {'CTE_6489', 'CTE_6912', 'CTE_7019', 'CTE_7126'};
nc_subs = {'NC_8095','NC_6839','NC_6974','NC_8653','NC_21499'};

%% Import the vascular metrics heatmap ROIs
% Path to top-level directory
dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/';
% Path to metrics
mpath = ['/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/metrics/' ...
    'gsigma_1-3-5_2-3-4_3-5-7_5-7-9_7-9-11/p18/'];
% Heatmaps struct - distribution of values for each region
roi_size = 1000;
hm_filename = append('heatmap_distro_',num2str(roi_size),'.mat');
hm = load(fullfile(mpath,hm_filename));
hm = hm.hm_distro;
% Regions to run LME model (just use the entire tissue)
regions = {'tiss','wm','gm','sulci','gyri','gm_sulci','wm_sulci',...
            'gm_gyri','wm_gyri'};
% Parameters
hm_params = {'vf','ld','bd','tort'};

%%% Reorganize heatmap data
% This script will reorganize the metrics struct by taking the region from
% each subject and combining them into a single substruct (e.g. taking the
% "tiss" region from each subject and creating a substruct named
% "metrics.tiss". This section can be skipped if it was already run.

% Iterate over each tissue region
for ii = 1:length(regions)
    % Iterate over heatmap parameters
    for j=1:length(hm_params)
        % Retrieve the heatmap ROI values for this region/parameter from
        % all AD, CTE, and NC subjects. Concatenate into arrays.
        [ad, cte, nc] = organize_metrics(hm,subids,...
                                         regions{ii},hm_params{j});
        % Save parameter by the group to fascilitate statistical analyses
        % in a later step.
        hm.(regions{ii}).(hm_params{j}).ad = ad;
        hm.(regions{ii}).(hm_params{j}).cte = cte;
        hm.(regions{ii}).(hm_params{j}).nc = nc;
    end
end

%% Calculate variances (rand. int & resid err) + group means for simulation

%%% Define combination of regions and vascular metrics
combos = struct();
% Length density regions
ld_reg = {'tiss','wm','gyri','sulci','wm_gyri','wm_sulci'};
% Branch density regions
bd_reg = {'wm_sulci'};
for ii = 1:length(ld_reg)
    combos.ld.(ld_reg{ii}) = 1;
end
for ii = 1:length(bd_reg)
    combos.bd.(bd_reg{ii}) = 1;
end

%%% Run the LME for each metric and region and each group
lme_var = calc_lme_model_variances(hm, regions, hm_params, subids);

%%% Define average number of observations per region
nobs = struct();
nobs.tiss = 2000;
nobs.wm = floor((2921/5 + 3317/4 + 3131/5)/3);
nobs.gyri = floor((1436/5 + 1368/4 + 922/5)/3);
nobs.sulci = floor((1240/5 + 1200/4 + 858/5)/3);
nobs.wm_sulci = floor((694/5 + 769/4 + 444/5)/3);
nobs.wm_gyri = floor((320/5 + 525/4 + 118/5)/3);

%% Run the power test & determine # tissue volumes

%%% Struct to store results
power_struct = struct();

%%% Paremters for simulation
alpha = 0.05;       % Significance level (alpha)
desired_power = 0.8;% Desired power of study
n_groups = 2;       % Number of groups (AD, CTE, HC)
max_sim = 100;     % Max # simulations for power estimation

%%% Run power test on significant results
met = fields(combos);

for ii = 1:length(met)
    % Extract regions
    reg = fields(combos.(met{ii}));
    % Iterate over each region
    for j = 1:length(reg)
        % Print region and metric
        fprintf('Starting metric %s & region %s\n',met{ii},reg{j});

        % Extract variance of random effect and residual error
        var_b = lme_var.ad.(reg{j}).(met{ii}).var_rint;
        var_e = lme_var.ad.(reg{j}).(met{ii}).var_err;
        % Extract group means
        ad_mean = lme_var.ad.(reg{j}).(met{ii}).mean;
        nc_mean = lme_var.nc.(reg{j}).(met{ii}).mean;
        
        %%% Iterate over reductions in group mean differences
        % Divisor to reduce the difference of group mean
        divisor = [1,2,4];
        % Structure name to track divisor
        div_name = {'div1','div2','div4'};
        for k = 1:3            
            %%% Reduce the group mean difference by divisor
            fprintf('Starting Divisor %s\n',div_name{k});
            d = abs(ad_mean - nc_mean);
            d = d ./ divisor(k);
            % Reduce the difference for all except first divisor
            if k>1
                if ad_mean < nc_mean
                    ad_mean = ad_mean + d;
                else
                    ad_mean = ad_mean - d;
                end
            end
            group_means = [ad_mean,nc_mean];
            % Extract number of observations
            obs_sub = nobs.(reg{j});
    
            %%% Function to find minimum number of tissue volumes
            [n_subjects,power,iteration] =...
                ptest(alpha,desired_power,n_groups,obs_sub,...
                      var_b,var_e,max_sim,group_means);
            fprintf('\n')
            % Store results
            power_struct.ad.(div_name{k}).(reg{j}).(met{ii}).n = n_subjects;
            power_struct.ad.(div_name{k}).(reg{j}).(met{ii}).power = power;
            power_struct.ad.(div_name{k}).(reg{j}).(met{ii}).iteration = iteration;
        end
    end

end

% for ii = 1:length(hm_params)
%     % Print to console
%     fprintf('Starting on parameter %s\n',hm_params{ii});
%     
%     %%% Extract the variances of random intercept and residual error
%     % TODO: extract from lme_var
%     % These are used for the sampling of the random effect (subject-level
%     % variability) and the error.
%     var_b = sig_params.(hm_params{ii}).var_b;
%     var_e = sig_params.(hm_params{ii}).var_e;
% 
%     %%% Extract group means (fixed effect)
%     % TODO: extract from lme_var
%     % These are are used to set the group-level mean when simulating Y
%     % e.g. [mean_AD, mean_HC] (note the order doesn't matter)
%     group_means = sig_params.(hm_params{ii}).group_means;
% 
%     %%% Average number of observations per region per subject
%     % TODO: change this to call the structure based on region
%     obs_sub = sig_params.(hm_params{ii}).obs_sub;
% 
%     %%% Function to find minimum number of tissue volumes
%     [n_subjects,power,iteration] =...
%         ptest(alpha,desired_power,n_groups,obs_sub,...
%               var_b,var_e,max_sim,group_means);
%     
%     % Store results
%     power_struct.(hm_params{ii}).n = n_subjects;
%     power_struct.(hm_params{ii}).power = power;
%     power_struct.(hm_params{ii}).iteration = iteration;
% end


% Save the results
save(fullfile(mpath,'power_analysis_reduced_difference.mat'),...
    'power_struct','-v7.3');

%% Power Test

function [n_subjects,power,iteration] =...
    ptest(alpha,desired_power,n_groups,obs_sub,var_b,var_e,max_sim,group_means)
% Simulate data to calculate minimum number of total observations across
% all groups to achieve the desired power and alpha.
%
%%% INPUTS
%   alpha (double): Significance level (alpha)
%   desired_power (double):  Desired power
%   n_groups (double): Number of groups (fixed effect) (AD, CTE, HC)
%   obs_sub (array): range of observations per subject [min, max] 
%   var_b (double): Variance of the random intercept (group-level effect).
%       This measures how much subjects differ from each other in their
%       baseline.
%   var_e (double): Residual error variance
%       This measures how much subjects differ from each other in their
%       baseline.
%   max_sim (double):  Max # simulations for power estimation
%   group_means (array): Mean values for each group
%                   (e.g., [mean_AD, mean_CTE, mean_HC])
%
%%% OUTPUTS:
%   n (double): number of observations
%   power (double): power of study
%   iteration (double): number of interations simulation ran through
%
% The simulation defines an LMM where:
%   Y: outcome (vascular metric observations)
%   X: disease group (AD, CTE, HC)
% The output will be the minimum number of observations required. These
% observations are the discrete regions of interest (ROIs), which are
% created when dividing the tissue volumes into subsamples.

% Set the parameters for the power calculation
max_n_subjects = 60;  % Maximum number of subjects
max_iterations = 20;  % Maximum iterations to avoid infinite loop
iteration = 0;  % Iteration counter

%%% Start with a small number of subjects
n_subjects_start = 10;  % Initial number of subjects
increment = 2;  % Increment to increase n_subjects after each iteration
n_subjects = n_subjects_start;  % Start with n_subjects_start subjects
power = 0;  % Initial power

%%% Loop to find the minimum n_subjects to achieve the desired power
while power < desired_power && n_subjects < max_n_subjects && iteration < max_iterations
    % Store results of p-values
    p_values = zeros(max_sim, 1);
    
    % Run simulations for the current number of subjects
    for i = 1:max_sim
        %%% Simulate random group assignments (1 to n_groups)
        % Ensure that the number of subjects is divisible by n groups
        if mod(n_subjects,n_groups) ~= 0
            error('Number subjects must be divisible by N groups')
        end
        % Calculate number of occurrences of each group assignment
        n_occur = n_subjects / n_groups;
        % Create vector of assignments
        group = repmat(1:n_groups,1,n_occur);
        % Shuffle assignments
        group = group(randperm(length(group)));
                
        % Create vectors for group and subjectIDs (same total number of rows as Y)
        group_all = [];  % All group values
        subjectIDs_all = [];  % All subject IDs
        Y_all = [];  % All outcome values (Y)

        % Simulate data for each subject
        for j = 1:n_subjects
            %%% Assign group and subjectID to each observations
            % Group repeated for each observation
            group_all = [group_all; repmat(group(j), obs_sub, 1)];
            % Subject ID repeated for each observation
            subjectIDs_all = [subjectIDs_all; repmat(j, obs_sub, 1)];
            
            %%% Fixed effect, random intercept, error
            % Get the group mean (fixed effect) for this group
            group_effect = group_means(group(j));
            % Simulate random intercept (u_j) for this subject
            u_j = randn(1) * sqrt(var_b);
            % Simulate residual errors (epsilon_ij) for each observation
            epsilon = randn(obs_sub, 1) * sqrt(var_e);
            
            %%% Calculate outcome variable (Y) for observation
            Y_all = [Y_all; group_effect + u_j + epsilon];  
        end
        
        % Create a table for the mixed model
        T = table(Y_all, group_all, subjectIDs_all);
        
        % Fit the linear mixed-effects model (using fitlme)
        try
            lme = fitlme(T, 'Y_all ~ group_all + (1|subjectIDs_all)');
        catch
            pause(0.1)
        end
        
        % Extract the p-value for the fixed effect of group
        p_values(i) = lme.Coefficients.pValue(2);
    end
    
    % Calculate the power as the proportion of significant p-values
    power = mean(p_values < alpha);
    
    % If the power is less than desired, increase the number of subjects
    if power < desired_power
        n_subjects = n_subjects + increment;  % Increase number of subjects
    end
    
    iteration = iteration + 1;  % Increment iteration count
end

% Display the result
fprintf('Minimum number of subjects required: %d\n', n_subjects);
fprintf('Achieved power: %.4f\n', power);
fprintf('Total iterations: %d\n', iteration);

end