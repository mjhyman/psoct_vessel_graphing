%% Power analysis to determine sample size for LME model
% The purpose is to determine the minimum sample size to achieve an alpha
% of 0.05 (p value) and a beta of 0.80 (power) for the linear mixed effects
% model (LMM). The literature states that the optimal method of calculating
% the minimum sample size is through simulations.
% 
% This script examines the myelin defects for AD vs. control. This model
% incorporates the post-mortem interval (PMI).

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
topdir = mydir(1:idcs(end));
addpath(genpath(topdir));

%% Import the vascular metrics heatmap ROIs
% Import spreadsheet of mus values
data_info = 'table 1.xlsx';
data = 'table 3.xlsx';
Info = readtable(data_info);
T = readtable(data);

% Separate table into subID, group, mus
subid = table2array(T(:,1));
grp_table = table2array(Info(:,2));
id_table = table2array(Info(:,1));
pmi_table = table2array(Info(:,3));
mus = table2array(T(:,3));

%%% Initialize group indices
% grp1 = control
% grp2 = AD
% grp3 = CTE
grp = zeros(size(subid));
pmi = zeros(size(subid));
for i=1:length(subid)
    % Find group index
    group_idx = find(id_table == subid(i));
    pmi(i) = pmi_table(group_idx);
    % Compare grp index
    if strcmp(grp_table(group_idx),'Control')
        grp(i) = 1;
    elseif strcmp(grp_table(group_idx),'AD')
        grp(i) = 2;
    elseif strcmp(grp_table(group_idx),'CTE')
        grp(i) = 3;
    end
end

% Extracat row identifier for each group by "Category" entry in table
grp1 = find(grp==1);
grp2 = find(grp==2);
grp3 = find(grp==3);

%%% Extract subID, mus for each group
% Extract subject ID for each group
grp1_subid = subid(grp1,:);
grp2_subid = subid(grp2,:);
grp3_subid = subid(grp3,:);
% Extract subject ID for each group
grp1_mus = mus(grp1,:);
grp2_mus = mus(grp2,:);
grp3_mus = mus(grp3,:);
% Extract PMI for each group
grp1_pmi = pmi(grp1,:);
grp2_pmi = pmi(grp2,:);
grp3_pmi = pmi(grp3,:);

%% Create LME model
%%% Create a table for fitting the LME model
% Column 1 = group label ('experimental' or 'control')
% Column 2 = value of each observation
% Column 3 = post-mortem interval
% grp1 = control
% grp2 = AD
% Create the group labels
g_cntrl = repmat({'CNTRL'},[length(grp1_mus),1]);
g_exp = repmat({'EXP'},[length(grp2_mus),1]);
grp_lbl = vertcat(g_exp, g_cntrl);
% Combine PMI into column vector
pmi = [grp2_pmi; grp1_pmi];
% Combine observations into column vector
obs = [grp2_mus; grp1_mus];
% Combine the subject IDs into column vector
subids = vertcat(grp2_subid, grp1_subid);
% Table (group labels, subjectID, vascular metric values)
tbl_ad_nc = table(grp_lbl,subids, pmi,obs,...
                  'VariableNames',{'Groups','subID','PMI','Observation'});

%% Calculate LMM variances
% Variances of Random Intercept & Resid Error
% Group Means of AD and normal control

%%% Fit linear mixed-effects model (LME) (table)
% Define the model:
%   response = measurement (Observation)
%   random effect (intercept/subject): subject identifier (subID)
%   fixed effect: group (experimental vs. control) (Groups)
%   fixed effect: post-mortem interval
fml = 'Observation ~ Groups + PMI + (1 | subID)';
% Fit the model for AD vs. HC
lme = fitlme(tbl_ad_nc,fml);
% Extract residual error variance
var_e = lme.MSE;
% Extract random intercept variance
[re, ~, ~] = randomEffects(lme);
% Take square of the standard error of the prediction of error
var_b = var(re);

%% Calculate Group means and Poisson distribution Lambda
% The group means are used for offseting the groups during random sampling.
% The Poisson Lambda is used for determining the shape of the Poisson
% distribution during random sampling.

%%% Group means of AD and NC
% NC group mean
nc_mean = mean(grp1_mus);
% AD group mean
ad_mean = mean(grp2_mus);
% array of group means
group_means = [ad_mean,nc_mean];

%%% Mean PMI of AD and NC combined
% This is used for random sampling the PMI in the simulation
mean_pmi = mean([grp1_pmi;grp2_pmi]);

%%% Poisson Lamba of AD and NC
% NC Lambda
nc_lambda = poissfit(grp1_mus);
% AD Lambda
ad_lambda = poissfit(grp2_mus);
% vector of lambda
group_lambda = [ad_lambda,nc_lambda];

%% Define average number of observations per region
% Maximum number of tissue slices per subject
nslice = 30;
% Maximum number of measurements per tissue slice
nobs = 25;
% Total number of observations
nobs = nslice .* nobs;
nobs = 8;

%% Run the power test & determine # tissue volumes

%%% Struct to store results
power_struct = struct();

%%% Paremters for simulation
alpha = 0.05;       % Significance level (alpha)
desired_power = 0.8;% Desired power of study
n_groups = 2;       % Number of groups (AD, CTE, HC)
max_sim = 100;     % Max # simulations for power estimation

%%% Run power test for AD vs. Control
fprintf('Starting power analysis of AD vs. Control.\n');
% Function to find minimum number of tissue volumes
[n_subjects,power,iteration] =...
    ptest(alpha,desired_power,n_groups,nobs,...
          var_b,var_e,max_sim,group_means,group_lambda,mean_pmi);
fprintf('Finished Power Analyses\n')
% Store results
power_struct.n = n_subjects;
power_struct.power = power;
power_struct.iteration = iteration;

% Save the results
save(fullfile(mydir,'power_analysis_myelin.mat'),...
    'power_struct','-v7.3');

%% Power Test

function [n_subjects,power,iteration] =...
    ptest(alpha,desired_power,n_groups,obs_sub,var_b,var_e,max_sim,...
          group_means,group_lambda,mean_pmi)
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
%                   (e.g., [mean_AD, mean_HC])
%   group_lambda (vector): Poisson distribution lambda
%                   (e.g., [AD_lambda, HC_lambda])
%   mean_pmi (double): the average PMI of the groups. This is used for
%           simulating the PMI for each subject
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
    fprintf('Starting iteration %d\n',iteration)
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
                
        %%% Group, PMI, subjectIDs (same total number of rows as Y)
        group_all = [];       % All group values
        pmi_all = [];         % All PMI values
        subjectIDs_all = [];  % All subject IDs
        Y_all = [];  % All outcome values (Y)

        %%% Simulate the PMI for each subject
        pmi = randn(n_subjects,1) + mean_pmi;

        % Simulate data for each subject
        for j = 1:n_subjects
            %%% Assign group and subjectID to each observations
            % Group repeated for each observation
            group_all = [group_all; repmat(group(j), obs_sub, 1)];
            % PMI for each observation
            pmi_all = [pmi_all; repmat(pmi(j), obs_sub, 1)];
            % Subject ID repeated for each observation
            subjectIDs_all = [subjectIDs_all; repmat(j, obs_sub, 1)];
            
            %%% Fixed effect, random intercept, error
            % Get the group mean (fixed effect) for this group
            group_effect = group_means(group(j));
            % Simulate random intercept (u_j) for this subject (Poisson)
            lambda = group_lambda(group(j));
            u_j = poissrnd(lambda,1,1) * sqrt(var_b);
            % Simulate residual errors (epsilon_ij) for each observation
            epsilon = randn(obs_sub, 1) * sqrt(var_e);
            
            %%% Calculate outcome variable (Y) for observation
            Y_all = [Y_all; group_effect + u_j + epsilon];  
        end
        
        % Create a table for the mixed model
        T = table(Y_all, group_all, pmi_all, subjectIDs_all);
        
        % Fit the linear mixed-effects model (using fitlme)
        try
            lme = fitlme(T, 'Y_all ~ group_all + pmi_all + (1|subjectIDs_all)');
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