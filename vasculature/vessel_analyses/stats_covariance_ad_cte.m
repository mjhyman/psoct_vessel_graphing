%% Measure covariance between age and average vascular metrics
% The purpose is to determine whether there is a covariance between each of
% the vascular metrics and the age of the tissue donors. This is to
% determine whether the age range of each group impacts the statistical
% outcomes.

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
subid = {'AD_10382', 'AD_20832', 'AD_20969', 'AD_21354', 'AD_21424',...
         'CTE_6489', 'CTE_6912', 'CTE_7019', 'CTE_7126',...
         'NC_6839',  'NC_6974',  'NC_8653',  'NC_21499', 'NC_301181'};
ad_subs = {'AD_10382', 'AD_20832', 'AD_20969', 'AD_21354', 'AD_21424'};
cte_subs = {'CTE_6489', 'CTE_6912', 'CTE_7019', 'CTE_7126'};
nc_subs = {'NC_8095','NC_6839','NC_6974','NC_8653','NC_21499'};


%% Import the vascular metrics
% Path to top-level directory
dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/';
% Path to metrics
mpath = ['/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/metrics/' ...
    'gsigma_1-3-5_2-3-4_3-5-7_5-7-9_7-9-11/p18/'];
% Load the metrics struct
metrics = load(fullfile(mpath,'metrics.mat'));
metrics = metrics.metrics;

%% Calculate covariance for each metric
% parameters to iterate over
params = {'length_density', 'branch_density', 'fraction_volume'};
% struct to store covariance
covariance = struct();
% struct to store correlation
correlation = struct();
% Scatter plot figure title 
fig_titles = {'Length Density','Branch Density', 'Fraction Volume'};
% y-axis units
ylabels = {'Length (\mum) / Volume (\mum^3)',...
        'Branches / Volume (\mum^3)','(a.u.)'};
% Marker Size
marker_size = 100;
% Font Size 
f_size = 18;

%%% Iterate over vascular metrics
for ii = 1:length(params)
    % Select specific vascular metric
    met = params{ii};
    
    % Create x/y arrays for the metrics
    x_ad = zeros(1,size(ad_subs,2));
    x_cte = zeros(1,size(cte_subs,2));
    x_nc = zeros(1,size(nc_subs,2));
    y_ad = x_ad;
    y_cte = x_cte;
    y_nc = x_nc;
    
    % Iterate over AD subject IDs
    for j = 1:length(ad_subs)
        % Extract subject ID
        sub = ad_subs{j};
        % Add the age and metric to the respective arrays
        x_ad(j) = metrics.(sub).age;
        y_ad(j) = metrics.(sub).tiss.(met);
    end

    % Iterate over CTE subject IDs
    for j = 1:length(cte_subs)
        % Extract subject ID
        sub = cte_subs{j};
        % Add the age and metric to the respective arrays
        x_cte(j) = metrics.(sub).age;
        y_cte(j) = metrics.(sub).tiss.(met);
    end

    % Iterate over NC subject IDs
    for j = 1:length(nc_subs)
        % Extract subject ID
        sub = nc_subs{j};
        % Add the age and metric to the respective arrays
        x_nc(j) = metrics.(sub).age;
        y_nc(j) = metrics.(sub).tiss.(met);
    end

    %%% Calculate covariance between age/metric for each group
    % AD covariance
    c = cov(x_ad, y_ad);
    covariance.(met).ad = c(1,2);
    % CTE covariance
    c = cov(x_cte, y_cte);
    covariance.(met).cte = c(1,2);
    % NC covariance
    c = cov(x_nc, y_nc);
    covariance.(met).nc = c(1,2);

    %%% Calculate correlation between age/metric for each group
    % AD
    [ad_corr, ad_p] = corrcoef(x_ad, y_ad);
    correlation.(met).ad = ad_corr(1,2);
    % CTE
    [cte_corr, cte_p] = corrcoef(x_cte, y_cte);
    correlation.(met).cte = cte_corr(1,2);
    % HC
    [nc_corr, nc_p] = corrcoef(x_nc, y_nc);
    correlation.(met).nc = nc_corr(1,2);

    %%% Create text string for reporting covariance
    ad_cov_str = append('AD = ', num2str(covariance.(met).ad));
    cte_cov_str = append('CTE = ', num2str(covariance.(met).cte));
    nc_cov_str = append('NC = ', num2str(covariance.(met).nc));
    cov_str = {'Covariances'; ad_cov_str; cte_cov_str; nc_cov_str};

    %%% Create text string for reporting correlation
    ad_corr_str = append('AD = ', num2str(correlation.(met).ad),...
        '. p = ', num2str(ad_p(1,2)));
    cte_corr_str = append('CTE = ', num2str(correlation.(met).cte),...
        '. p = ', num2str(cte_p(1,2)));
    nc_corr_str = append('NC = ', num2str(correlation.(met).nc),...
        '. p = ', num2str(nc_p(1,2)));
    corr_str = {'Correlations'; ad_corr_str; cte_corr_str; nc_corr_str};

    %% Scatter Plots of Covariance and Correlation

    %%% Overlay of scatter plots
    % Initialize the scatter plot
    fh = figure();
    fh.WindowState = 'maximized';
    % Overlay scatter plots for AD, CTE, NC
    scatter(x_ad, y_ad, marker_size,'square','b','filled'); hold on;
    scatter(x_cte, y_cte, marker_size,'o','r','filled');
    scatter(x_nc, y_nc, marker_size,'diamond','k','filled');
    % Textbox with covariance and correlation
    dim = [.2 .5 .3 .3];
    annotation('textbox',dim, 'String',cov_str,'FitBoxToText','on',...
        'FontSize',f_size);
    dim = [.35 .5 .3 .3];
    annotation('textbox',dim, 'String',corr_str,'FitBoxToText','on',...
        'FontSize',f_size);
    % Format the figure
    xlabel('Age (years)')
    ylabel(ylabels{ii})
    legend({'AD', 'CTE', 'NC'})
    title(fig_titles{ii})
    set(gca, 'FontSize', f_size)
    % Save the figure
    fout = append('cov_corr_',met);
    fout = fullfile(mpath, fout);
    saveas(gca,fout,'png');

    %%% Create subplots and overlay with least squares fitted line
    fh = figure();
    fh.WindowState = 'maximized';
    
    % Scatter plots of values for AD, CTE, HC
    ax1 = subplot(3,1,1);
    s1 = scatter(ax1, x_ad, y_ad, marker_size,'square','b','filled');
    title('AD'); set(gca, 'FontSize', f_size)
    ax2 = subplot(3,1,2);
    s2 = scatter(ax2, x_cte, y_cte, marker_size,'o','r','filled');
    ylabel(ylabels{ii})
    title('CTE'); set(gca, 'FontSize', f_size)
    ax3 = subplot(3,1,3);
    s3 = scatter(ax3, x_nc, y_nc, marker_size,'diamond','k','filled');
    title('HC');
    
    % Overlay LS line for each group
    h = lsline(ax1); h.LineWidth = 3; h.Color = s1.CData;
    h = lsline(ax2); h.LineWidth = 3; h.Color = s2.CData;
    h = lsline(ax3); h.LineWidth = 3; h.Color = s3.CData;
    
    % Format the figure
    xlabel('Age (years)')
    set(gca, 'FontSize', f_size)
    % Save the figure
    fout = append('cov_corr_lsline_',met);
    fout = fullfile(mpath, fout);
    saveas(gca,fout,'png');
end