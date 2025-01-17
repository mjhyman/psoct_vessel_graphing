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
% average value parameters
params = {'length_density','branch_density',...
    'fraction_volume','tortuosity'};
% Load the Heatmap ROI struct
hm = load(fullfile(mpath,'heatmap_distro_1000.mat'));
hm = hm.hm_distro;
% heatmap vascular metrics to iterate over
mets = {'ld','bd','vf','tort'};
% Regions
regions = {'tiss','wm','gm','sulci','gyri','gm_sulci','wm_sulci',...
    'gm_gyri','wm_gyri'};

%% Correlation average vascular metric vs. age
% struct to store covariance
covariance = struct();
% struct to store correlation
correlation = struct();

%%% Iterate over vascular metrics
for ii = 1:length(params)
    % Select specific vascular metric
    met = params{ii};
    for k = 1:length(regions)       
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
            y_ad(j) = mean(metrics.(sub).(regions{k}).(met));
        end
    
        % Iterate over CTE subject IDs
        for j = 1:length(cte_subs)
            % Extract subject ID
            sub = cte_subs{j};
            % Add the age and metric to the respective arrays
            x_cte(j) = metrics.(sub).age;
            y_cte(j) = mean(metrics.(sub).(regions{k}).(met));
        end
    
        % Iterate over NC subject IDs
        for j = 1:length(nc_subs)
            % Extract subject ID
            sub = nc_subs{j};
            % Add the age and metric to the respective arrays
            x_nc(j) = metrics.(sub).age;
            y_nc(j) = mean(metrics.(sub).(regions{k}).(met));
        end
    
        %%% Calculate covariance between age/metric for each group
        % AD covariance
        c = cov(x_ad, y_ad);
        covariance.(met).(regions{k}).ad = c(1,2);
        % CTE covariance
        c = cov(x_cte, y_cte);
        covariance.(met).(regions{k}).cte = c(1,2);
        % NC covariance
        c = cov(x_nc, y_nc);
        covariance.(met).(regions{k}).nc = c(1,2);
    
        %%% Calculate correlation between age/metric for each group
        % AD
        [ad_corr, ad_p] = corrcoef(x_ad, y_ad);
        correlation.(met).(regions{k}).ad = ad_corr(1,2);
        % CTE
        [cte_corr, cte_p] = corrcoef(x_cte, y_cte);
        correlation.(met).(regions{k}).cte = cte_corr(1,2);
        % HC
        [nc_corr, nc_p] = corrcoef(x_nc, y_nc);
        correlation.(met).(regions{k}).nc = nc_corr(1,2);
        % Print if correlation at least moderate
        if ad_corr(1,2) >= 0.7
            sprintf('AD correlation = %f for %s in %s',ad_corr(1,2),...
                met,regions{k})
            figure; scatter(x_ad,y_ad); title({'AD',met,regions{k}});
        elseif cte_corr(1,2) >= 0.7
            sprintf('CTE correlation = %f for %s in %s',cte_corr(1,2),...
                met,regions{k})
            figure; scatter(x_cte,y_cte); title({'CTE',met,regions{k}});
        elseif nc_corr(1,2) >= 0.7
            sprintf('NC correlation = %f for %s in %s',nc_corr(1,2),...
                met,regions{k})
            figure; scatter(x_nc,y_nc); title({'NC',met,regions{k}});
        end
    end
end

%% Correlation b/w vascular metric heatmap
% X = ROIs from vascular metric A
% Y = ROIs from vascular metric B
% measure correlation between X,Y
% repeat for each pair of metrics (24 combinations)

%% Correlation heatmap vascular metric vs. age
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
    % Iterate over regions
    for j = 1:length(regions)
        %%% Initialize arrays to store x,y values
        % Count number of values in each vascular metric array. This is
        % used to initialize the age arrays.
        ad_roi = 0; cte_roi = 0; hc_roi = 0;
        for k = 1:length(ad_subs)
            % Add number of ROIs to counter
            ad_roi = ad_roi + length(hm.(ad_subs{k}).(regions{j}).(met));
            hc_roi = hc_roi + length(hm.(nc_subs{k}).(regions{j}).(met));
            try
                cte_roi = cte_roi + length(hm.(cte_subs{k}).(regions{j}).(met));
            catch
                continue;
            end
        end
        % Initialize x and y vectors for each group
        x_ad = zeros(ad_roi,1);     y_ad = zeros(ad_roi,1);
        x_cte = zeros(cte_roi,1);   y_cte = zeros(cte_roi,1);
        x_hc = zeros(hc_roi,1);     y_hc = zeros(hc_roi,1);
        
        %%% Import AD, CTE, HC subjects
        [x_ad, y_ad] = consolidate_group(hm,metrics,ad_subs,...
            (regions{j}),met,x_ad,y_ad);
        [x_cte, y_cte] = consolidate_group(hm,metrics,cte_subs,...
            (regions{j}),met,x_cte,y_cte);
        [x_hc, y_hc] = consolidate_group(hm,metrics,nc_subs,...
            (regions{j}),met,x_hc,y_hc);
            
        %%% Correlation between age & vascular metric
        % AD
        [ad_corr, ad_p] = corrcoef(x_ad, y_ad);
        hm_corr_cov.(met).(regions{j}).ad.corr = ad_corr(1,2);
        hm_corr_cov.(met).(regions{j}).ad.p = ad_p(1,2);
        % CTE
        [cte_corr, cte_p] = corrcoef(x_cte, y_cte);
        hm_corr_cov.(met).(regions{j}).cte.corr = cte_corr(1,2);
        hm_corr_cov.(met).(regions{j}).cte.p = cte_p(1,2);
        % HC
        [nc_corr, nc_p] = corrcoef(x_hc, y_hc);
        hm_corr_cov.(met).(regions{j}).hc.corr = nc_corr(1,2);
        hm_corr_cov.(met).(regions{j}).hc.p = nc_p(1,2);
        % AD, CTE, NC combined correlation
        x = [x_ad;x_cte;x_hc];
        y = [y_ad;y_cte;y_hc];
        [all_corr, all_p] = corrcoef(x, y);
        hm_corr_cov.(met).(regions{j}).all.corr = all_corr(1,2);
        hm_corr_cov.(met).(regions{j}).all.p = all_p(1,2);
        % Print if correlation at least moderate
        if ad_corr >= 0.3
            sprintf('AD correlation >= 0.5 for %s in %s',met, regions{j})
        elseif cte_corr >= 0.3
            sprintf('CTE correlation >= 0.5 for %s in %s',met, regions{j})
        elseif nc_corr >= 0.3
            sprintf('NC correlation >= 0.5 for %s in %s',met, regions{j})
        elseif all_corr >= 0.3
            sprintf('ALL correlation >= 0.5 for %s in %s',met, regions{j})
        end

        %%% Covariance between age & vascular metric
        % AD covariance
        c = cov(x_ad, y_ad);
        hm_corr_cov.(met).(regions{j}).ad.cov = c(1,2);
        % CTE covariance
        c = cov(x_cte, y_cte);
        hm_corr_cov.(met).(regions{j}).cte.cov = c(1,2);
        % NC covariance
        c = cov(x_hc, y_hc);
        hm_corr_cov.(met).(regions{j}).hc.cov = c(1,2);
        % AD, CTE, NC combined covariance
        c = cov(x,y);
        hm_corr_cov.(met).(regions{j}).all.cov = c(1,2);
  
    
        %%% Univariate Analysis (sex for AD)
        
    
        %%% Multivariate Linear Regression (sex + age for AD)    
        
    
        %% Scatter Plots of Covariance and Correlation
        %{
        %%% Create text string for reporting covariance
        ad_cov_str = append('AD = ', num2str(hm_corr_cov.(met).(regions{j}).ad.cov));
        cte_cov_str = append('CTE = ', num2str(hm_corr_cov.(met).(regions{j}).cte.cov));
        nc_cov_str = append('NC = ', num2str(hm_corr_cov.(met).(regions{j}).hc.cov));
        cov_str = {'Covariances'; ad_cov_str; cte_cov_str; nc_cov_str};
    
        %%% Create text string for reporting correlation
        ad_corr_str = append('AD = ', num2str(hm_corr_cov.(met).(regions{j}).ad.corr),...
            '. p = ', num2str(ad_p(1,2)));
        cte_corr_str = append('CTE = ', num2str(hm_corr_cov.(met).(regions{j}).cte.corr),...
            '. p = ', num2str(cte_p(1,2)));
        nc_corr_str = append('NC = ', num2str(hm_corr_cov.(met).(regions{j}).hc.corr),...
            '. p = ', num2str(nc_p(1,2)));
        corr_str = {'Correlations'; ad_corr_str; cte_corr_str; nc_corr_str};
    
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
        %}
    end
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