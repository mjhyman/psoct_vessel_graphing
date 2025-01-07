%% Correlation b/w age and heatmap vascular metrics
% The purpose is to determine whether there is a covariance between each of
% the vascular metrics and the age of the tissue donors. This is to
% determine whether the age range of each group impacts the statistical
% outcomes. This is performed on the heatmap regions of interest (ROIs).

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

%% Import the vascular metrics (heatmap ROIs)
% Path to top-level directory
dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/';
% Path to metrics
mpath = ['/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/metrics/' ...
    'gsigma_1-3-5_2-3-4_3-5-7_5-7-9_7-9-11/p18/'];
% Load the metrics
metrics = load(fullfile(mpath,'metrics.mat'));
metrics = metrics.metrics;
% Load the Heatmap ROI struct
hm = load(fullfile(mpath,'heatmap_distro_1000.mat'));
hm = hm.hm_distro;

%% Calculate covariance for each metric
% vascular metrics to iterate over
mets = {'ld','bd','vf','tort'};
% struct to store covariance
covariance = struct();
% struct to store correlation
hm_corr_cov = struct();
% Scatter plot figure title 
fig_titles = {'Length Density','Branch Density',...
               'Volume Fraction','Tortuosity'};
% y-axis units
ylabels = {'Length (\mum) / Volume (\mum^3)',...
           'Branches / Volume (\mum^3)','(unitless)','Tortuosity'};
% Marker Size
marker_size = 100;
% Font Size 
f_size = 18;

%%% Iterate over vascular metrics
for ii = 1:length(mets)
    % Select specific vascular metric
    met = mets{ii};

    %%% Initialize arrays to store x,y values
    % Count number of values in each vascular metric array. This is used to
    % initialize the age arrays.
    ad_roi = 0; cte_roi = 0; hc_roi = 0;
    for j = 1:length(ad_subs)
        % Add number of ROIs to counter
        ad_roi = ad_roi + length(hm.(ad_subs{j}).tiss.(met));
        hc_roi = hc_roi + length(hm.(nc_subs{j}).tiss.(met));
        try
            cte_roi = cte_roi + length(hm.(cte_subs{j}).tiss.(met));
        catch
            continue;
        end
    end

    x_ad = zeros(ad_roi,1);     y_ad = zeros(ad_roi,1);
    x_cte = zeros(cte_roi,1);   y_cte = zeros(cte_roi,1);
    x_hc = zeros(hc_roi,1);     y_hc = zeros(hc_roi,1);
    
    %%% Import AD, CTE, HC subjects
    [x_ad, y_ad] = consolidate_group(hm,metrics,ad_subs,'tiss',met,x_ad,y_ad);
    [x_cte, y_cte] = consolidate_group(hm,metrics,cte_subs,'tiss',met,x_cte,y_cte);
    [x_hc, y_hc] = consolidate_group(hm,metrics,nc_subs,'tiss',met,x_hc,y_hc);
    
    %%% Visualize (to determine regression type)
    % Plot AD
    figure; tstr = strcat(fig_titles{ii}, ' AD');
    scatter(x_ad,y_ad,'filled','o');
    title(tstr); ylabel(ylabels{ii}); xlabel('Age (years)');
    % Plot CTE
    figure; tstr = strcat(fig_titles{ii}, ' CTE');
    scatter(x_cte,y_cte,'filled','o');
    title(tstr); ylabel(ylabels{ii}); xlabel('Age (years)');
    % Plot HC
    figure; tstr = strcat(fig_titles{ii}, ' HC');
    scatter(x_hc,y_hc,'filled','o');
    title(tstr); ylabel(ylabels{ii}); xlabel('Age (years)');
        
    %%% Correlation between age & vascular metric
    % AD
    [ad_corr, ad_p] = corrcoef(x_ad, y_ad);
    hm_corr_cov.(met).ad.corr = ad_corr(1,2);
    hm_corr_cov.(met).ad.p = ad_p(1,2);
    % CTE
    [cte_corr, cte_p] = corrcoef(x_cte, y_cte);
    hm_corr_cov.(met).cte.corr = cte_corr(1,2);
    hm_corr_cov.(met).cte.p = cte_p(1,2);
    % HC
    [nc_corr, nc_p] = corrcoef(x_hc, y_hc);
    hm_corr_cov.(met).hc.corr = nc_corr(1,2);
    hm_corr_cov.(met).hc.p = nc_p(1,2);
    % AD, CTE, NC combined correlation
    x = [x_ad;x_cte;x_hc];
    y = [y_ad;y_cte;y_hc];
    [all_corr, all_p] = corrcoef(x, y);
    hm_corr_cov.(met).all.corr = all_corr(1,2);
    hm_corr_cov.(met).all.p = all_p(1,2);
    
    %%% Covariance between age & vascular metric
    % AD covariance
    c = cov(x_ad, y_ad);
    hm_corr_cov.(met).ad.cov = c(1,2);
    % CTE covariance
    c = cov(x_cte, y_cte);
    hm_corr_cov.(met).cte.cov = c(1,2);
    % NC covariance
    c = cov(x_hc, y_hc);
    hm_corr_cov.(met).hc.cov = c(1,2);
    % AD, CTE, NC combined covariance
    c = cov(x,y);
    hm_corr_cov.(met).all.cov = c(1,2);

    %%% Regression
    beta0 = [1,1,1,1,1];
    b = nlinfit(x_ad',y_ad',@hougen,beta0');

    %%% Univariate Analysis (sex for AD)
    

    %%% Multivariate Linear Regression (sex + age for AD)    
    

    %%% Calculate correlation between age/metric for each group
    

    %%% Create text string for reporting covariance
    ad_cov_str = append('AD = ', num2str(covariance.(met).ad));
    cte_cov_str = append('CTE = ', num2str(covariance.(met).cte));
    nc_cov_str = append('NC = ', num2str(covariance.(met).hc));
    cov_str = {'Covariances'; ad_cov_str; cte_cov_str; nc_cov_str};

    %%% Create text string for reporting correlation
    ad_corr_str = append('AD = ', num2str(hm_corr_cov.(met).ad),...
        '. p = ', num2str(ad_p(1,2)));
    cte_corr_str = append('CTE = ', num2str(hm_corr_cov.(met).cte),...
        '. p = ', num2str(cte_p(1,2)));
    nc_corr_str = append('NC = ', num2str(hm_corr_cov.(met).hc),...
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
    scatter(x_hc, y_hc, marker_size,'diamond','k','filled');
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
    s3 = scatter(ax3, x_hc, y_hc, marker_size,'diamond','k','filled');
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

%% Extract age and ROIs from specific group
function [x, y] = consolidate_group(hm,metrics,subjects,region,metric,x,y)
% CONSOLIDATE_GROUP combine ROIs from all subjects into single array
% Iterate over the heatmap to extract all ROIs. Iterate over the "metrics"
% structure to retrieve the age of the subject.
% INPUTS:
%   hm (struct): heatmap ROIs
%   metrics (struct): metrics structure
%   subjects (cell array): subject IDs
%   region (string): region of sample (i.e. tiss, sulci, gyri)
%   metric (string): vascular metric (vf, ld, bd, tort)
%   x (double array): zeros array to store age at death of each ROI
%   y (double array): zeros array to store the ROI value
% OUTPUTS:
%   x (double array): age of each ROI
%   y (double array): value of each ROI

% Initialize variable for storing index
idx = 1;
for j = 1:length(subjects)
    % Extract all ROI values
    roi = hm.(subjects{j}).(region).(metric);
    y(idx:idx+length(roi)-1) = roi;
    % Add age to x-axis for each entry
    age = metrics.(subjects{j}).age;
    x(idx:idx+length(roi)-1) = age;
    % Increment counter
    idx = idx + length(roi);
end

end