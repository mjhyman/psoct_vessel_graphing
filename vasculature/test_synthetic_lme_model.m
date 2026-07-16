% Load pilot data (assuming a dataset with columns for the outcome, group, and subject ID)
% Example data format: 
% Outcome variable 'Y', Group variable 'Group' (categorical), Subject ID 'SubjectID'

% Create a sample pilot dataset (replace with your actual data)
% Assuming we have three groups and several subjects
subjects = repmat(1:10, 1, 3);  % 10 subjects with 3 observations each
group = repmat([1, 2, 3], 10, 1);  % 3 groups per subject
group = group(:);  % Reshape to make it a column vector
Y = randn(30, 1) + group;  % Outcome variable with random noise + group effect

% Subject IDs (each subject has multiple observations)
subjectIDs = repmat(1:10, 1, 3);  % Subject IDs for each observation

% Create a table for fitting the mixed model
T = table(Y, group, subjectIDs','VariableNames',{'Y','group','subjectIDs'});

% Fit a linear mixed-effects model
% Random intercept for SubjectID, fixed effect for Group
lme = fitlme(T, 'Y ~ group + (1|subjectIDs)');

% Display the model results
disp(lme);

% Extract the residual error variance (MSE) from the model fit
residual_error_var = lme.MSE;  % Residual error variance (MSE)

% Extract random effects and standard errors using randomEffects function
[~, ~, stats] = randomEffects(lme);

% Calculate the variance of the random intercept
random_intercept_var = (stats.SEPred) .^ 2;  % Squaring the standard error

% Display the estimated variances
fprintf('Estimated random intercept variance: %.4f\n', random_intercept_var);
fprintf('Estimated residual error variance: %.4f\n', residual_error_var);
