function [lme_var] = calc_lme_model_variances(hm, regions, params, subids)
%CALC_HEATMAP_STATS Linear Mixed Effects (LME) model
% This function creates an LME model and then tests the fixed effects
% (diseased vs. control) while accounting for random effects within groups.
% This method accounts for multiple observations (samples) from the same
% subject.
% 
% In addition, this script measures the random intercept variance and the
% residual error variance.
%
%   INPUTS:
%       hm (struct): contains substructures of vascular parameters
%                       metrics.[region].[parameter].[group]
%       regions (cell array): tissue regions (tiss, gyri, sulci, gm, wm)
%       params (cell array): vascular parameters (length_density,
%                            branch_density, fraction_volume, tortuosity
%                            diameter)
%       subids (cell array): subject IDs
%   OUTPUTS:
%       lme_var (struct): contains random intercept variance and variance
%       of random intercept.
%
%% Initialize workspace
% Struct for storing variances of random intercept and residual error
lme_var = struct();

% Iterate over tissue regions
for ii = 1:length(regions)
    % Iterate over parameters
    for j = 1:length(params) 
        %%% Load in the arrays from each group (ad, cte, nc)
        ad = hm.(regions{ii}).(params{j}).ad;
        cte = hm.(regions{ii}).(params{j}).cte;
        nc = hm.(regions{ii}).(params{j}).nc;
        % Calculate mean values of each (used to set fixed effect)
        lme_var.ad.(regions{ii}).(params{j}).mean = mean(ad);
        lme_var.cte.(regions{ii}).(params{j}).mean = mean(cte);
        lme_var.nc.(regions{ii}).(params{j}).mean = mean(nc);

        %% Generate a linear mixed-effects model
        % The purpose is to incorporate a random effect into the analysis
        % to account for a potential correlation between outcomes on the
        % same person.

        %%% Create column vector of subject IDs for each group
        % Initialize counter for N values per group
        nval_ad = 0;
        nval_cte = 0;
        nval_nc = 0;
        % Initialize counter for number of values per subject
        nval = zeros(length(subids),1);
        % Iterate over subjects
        for s = 1:length(subids)
            % Count number of values in array
            tmp = hm.(subids{s}).(regions{ii}).(params{j});
            nval(s) = length(tmp);
            % Add N values to the respective group list
            if contains(subids(s), 'AD_')
                nval_ad = nval_ad + nval(s);
            elseif contains(subids(s), 'CTE_')
                nval_cte = nval_cte + nval(s);
            elseif contains(subids(s), 'NC_')
                nval_nc = nval_nc + nval(s);
            end
        end
        % Create column vector with subject identifier for each value
        subID_vec = zeros(sum(nval),1);
        % Start index
        s0 = 1;
        for s = 1:length(subids)
            % Update end index
            sf = s0 + nval(s) - 1;
            % Add subject ID for each value
            subID_vec(s0:sf) = ones(nval(s),1) .* s;
            % Update start index
            s0 = sf + 1;
        end
        % Separate the groups (AD, CTE, NC)
        ad_vec = subID_vec(1 : nval_ad);
        cte_vec = subID_vec(nval_ad+1 : nval_ad+nval_cte);
        nc_vec = subID_vec(nval_ad+nval_cte+1 : nval_ad+nval_cte+nval_nc);
        
        %%% LME model of AD vs. HC
        % Column 1 = group label
        % Column 2 = vascular metric value
        % Column 3 = vascular metric
        g_ad = repmat({'AD'},[length(ad),1]);
        g_cte = repmat({'CTE'},[length(cte),1]);
        g_nc = repmat({'HC'},[length(nc),1]);
        groups_ad_nc = vertcat(g_ad,g_nc);
        subids_ad_nc = vertcat(ad_vec, nc_vec);
        vmets_ad_nc = [ad; nc];
        tbl_ad_nc = table(groups_ad_nc,subids_ad_nc,vmets_ad_nc,...
                  'VariableNames',{'Groups','subID','VascularMetric'});
        
        %%% LME model of CTE vs. HC
        groups_cte_nc = vertcat(g_cte, g_nc);
        subids_cte_nc = vertcat(cte_vec, nc_vec);
        vmets_cte_nc = [cte; nc];
        tbl_cte_nc = table(groups_cte_nc,subids_cte_nc,vmets_cte_nc,...
                  'VariableNames',{'Groups','subID','VascularMetric'});

        %%% Return CTE vs HC and AD vs. HC
        [ad_var_err, ad_var_rint] = lme_vars(tbl_ad_nc);
        [cte_var_err, cte_var_rint] = lme_vars(tbl_cte_nc);
        % Add to struct
        lme_var.ad.(regions{ii}).(params{j}).var_err = ad_var_err;
        lme_var.ad.(regions{ii}).(params{j}).var_rint = ad_var_rint;
        lme_var.cte.(regions{ii}).(params{j}).var_err = cte_var_err;
        lme_var.cte.(regions{ii}).(params{j}).var_rint = cte_var_rint;            
    end
end

    function [var_err, var_rint] = lme_vars(tbl)
        %%% Fit linear mixed-effects model (LME) (table)
        % Define the model:
        %   response = vascular metric array
        %   random effect (intercept/subject): subID
        %   fixed effect: Groups
        fml = 'VascularMetric ~ Groups + (1 | subID)';
        % Fit the model for AD vs. HC
        lme = fitlme(tbl,fml);

        %%% Extract residual error variance
        var_err = lme.MSE;

        %%% Extract random intercept variance
        [re, ~, ~] = randomEffects(lme);
        % Take square of the standard error of the prediction of error
        var_rint = var(re);

    end

end