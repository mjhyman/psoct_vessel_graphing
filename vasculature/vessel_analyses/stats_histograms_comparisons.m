%% Histograms to investigate the discrepancy in p-values from the LME
% LME model outputs the following p-values for length density (AD vs HC):
% gyri = 6.2e-7
% wm gyri = 0.013
% gm > 0.05
%
% Since the GM is not significant, it is expected that the p-value of wm
% gyri should also be lower than gyri. However the inverse is true. This
% script will plot the histograms to determine why.
%
% TODO:
% - Tortuosity
%   - volume, wm, gm + gyri, WM gyri, GM gyri


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
% Import LME statistics
load(fullfile(mpath, 'heatmap_lme_stats.mat'));

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


%% Histograms of length density (Volume + Gyri) for AD vs. control
% x and y-axis limits
xlims = [0.25,20]; ylims1 = [0,450]; ylims2 = [0, 100];
% Bin width (units of vascular metric)
bw = 0.25;
% make histogram
generate_histograms(hm,hm_stats,mpath,'ld',xlims,ylims1,ylims2,bw,...
    'tiss','wm','gm','gyri','wm_gyri','gm_gyri',...
    {'Length Density - Entire Volume';'Length Density - White Matter';...
    'Length Density - Gray Matter';'Length Density - Gyri';...
    'Length Density - WM Gyri';'Length Density - GM Gyri'},...
    'mm^-^2','LMM_histograms_LengthDens_vol_wm_gm_gyri.png','ad');

%% Histograms of length density (Volume + Sulci) for AD vs. control
% x and y-axis limits
xlims = [0.25,20]; ylims1 = [0,450]; ylims2 = [0, 100];
% Bin width (units of vascular metric)
bw = 0.25;
% make histogram
generate_histograms(hm,hm_stats,mpath,'ld',xlims,ylims1,ylims2,bw,...
    'tiss','wm','gm','sulci','wm_sulci','gm_sulci',...
    {'Length Density - Entire Volume';'Length Density - White Matter';...
    'Length Density - Gray Matter';'Length Density - Sulci';...
    'Length Density - WM Sulci';'Length Density - GM Sulci'},...
    'mm^-^2','LMM_histograms_LengthDens_vol_wm_gm_sulci.png','ad');

%% Histogram of Branch Density (volume + sulci) for AD vs. control
% x and y-axis limits
xlims = [0.25,150]; ylims1 = [0,700]; ylims2 = [0, 150];
% Bin width (units of vascular metric)
bw = 5;
% make histogram
generate_histograms(hm,hm_stats,mpath,'bd',xlims,ylims1,ylims2,bw,...
    'tiss','wm','gm','sulci','wm_sulci','gm_sulci',...
    {'Branch Density - Entire Volume';'Branch Density - White Matter';...
    'Branch Density - Gray Matter';'Branch Density - Sulci';...
    'Branch Density - WM Sulci';'Branch Density - GM Sulci'},...
    'mm^-^3','LMM_histograms_BranchDens_vol_wm_gm_sulci.png','ad');

%% Histogram of Tortuosity (volume + gyri) for CTE vs. control
% TODO: change from AD to CTE in both function and call to function
% x and y-axis limits
xlims = [1,1.5]; ylims1 = [0,700]; ylims2 = [0, 150];
% Bin width (units of vascular metric)
bw = 0.01;
% make histogram
generate_histograms(hm,hm_stats,mpath,'tort',xlims,ylims1,ylims2,bw,...
    'tiss','wm','gm','gyri','wm_gyri','gm_gyri',...
    {'Tortuosity - Entire Volume';'Tortuosity - White Matter';...
    'Tortuosity - Gray Matter';'Tortuosity - Gyri';...
    'Tortuosity - WM Gyri';'Tortuosity - GM Gyri'},...
    '(unitless)','LMM_histograms_Tortuosity_vol_wm_gm_gyri.png','cte');


%% Function of histogram plots

function generate_histograms(hm,hm_stats,mpath,met,xlims,ylims1,ylims2,...
    bw,reg1, reg2, reg3, reg4, reg5, reg6, tstr, xlab, fname, group)
% Create 3x2 figure of histograms 
% INPUTS:
%   hm (struct): heatmap structure
%   hm_stats (struct): heatmap statistics structure
%   mpath (string): path to metrics output folder
%   met (string): vascular metric name
%   xlims (double vector): x-axis limits
%   ylims1 (double vector): y-axis limits for top row
%   ylims2 (double vector): y-axis limits for bottom row
%   bw (double): bin width
%   reg1-6 (string): regions 1-6 for histogram
%   tstrings (cell array): cell array of title strings
%   xlab (string): x-axis label
%   fname (string): name of output file
%   group (string): group for comparison with control ('ad' or 'cte')

%%% Textbox dimensions for storing p-value
dim1 = [0.207251301467108,0.686815568786649,...
    0.083481800597072,0.187179191038678];
dim2 = [0.48332809446773,0.679996455066376,...
    0.076890465758453,0.187179191038678];
dim3 = [0.7588,0.6661,0.0835,0.1872];
dim4 = [0.218,0.19,0.083,0.1872];
dim5 = [0.482536162411161,0.193,0.0835,0.1872];
dim6 = [0.7598,0.2067,0.0835,0.1872];

%%% Select the statistical comparisons
if strcmp(group,'ad')
    grp = 'ad_nc';
else
    grp = 'cte_nc';
end

%%% Region 1
figure;
subplot(2,3,1)
ad = hm.(reg1).(met).(group);
hc = hm.(reg1).(met).nc;
histogram(ad,'FaceColor','r','BinWidth',bw); hold on;
histogram(hc,'FaceColor','g','BinWidth',bw);
xlim([xlims(1),xlims(2)]); ylim([ylims1(1),ylims1(2)])
title(tstr{1},'FontSize',20);
xlabel(xlab); ylabel('Counts',"Rotation",90)
% Annotations
bhat = hm_stats.(reg1).(met).coef.(grp);
se = hm_stats.(reg1).(met).se.(grp);
t = hm_stats.(reg1).(met).t.(grp);
df = hm_stats.(reg1).(met).df.(grp);
pval = hm_stats.(reg1).(met).p.(grp);
txt = {[strcat('$\hat{\beta}$ = ',num2str(bhat,2))],...
        [strcat('se = ',num2str(se,2))],...
        [strcat('t = ',num2str(t,2))],...
        [strcat('df = ',num2str(df,2))],...
        [strcat('p = ',num2str(pval,2))]};
t=annotation('textbox',dim1,'String',txt,'FitBoxToText','on',...
    'Interpreter','latex');
t.FontSize = 20;
xline(mean(ad),'Color','r','LineStyle',':','LineWidth',4);
xline(mean(hc),'Color','g','LineStyle',':','LineWidth',4);
if strcmp(group,'ad')
    legend({'AD','HC','',''});
else
    legend({'CTE','HC','',''});
end
set(gca,'FontSize',20)

%%% Region 2
subplot(2,3,2)
ad = hm.(reg2).(met).(group);
hc = hm.(reg2).(met).nc;
histogram(ad,'FaceColor','r','BinWidth',bw); hold on;
histogram(hc,'FaceColor','g','BinWidth',bw);
xlim([xlims(1),xlims(2)]); ylim([ylims1(1),ylims1(2)])
title(tstr{2},'FontSize',20);
xlabel(xlab); ylabel('Counts',"Rotation",90)
% Annotations
bhat = hm_stats.(reg2).(met).coef.(grp);
se = hm_stats.(reg2).(met).se.(grp);
t = hm_stats.(reg2).(met).t.(grp);
df = hm_stats.(reg2).(met).df.(grp);
pval = hm_stats.(reg2).(met).p.(grp);
txt = {[strcat('$\hat{\beta}$ = ',num2str(bhat,2))],...
        [strcat('se = ',num2str(se,2))],...
        [strcat('t = ',num2str(t,2))],...
        [strcat('df = ',num2str(df,2))],...
        [strcat('p = ',num2str(pval,2))]};
t=annotation('textbox',dim2,'String',txt,'FitBoxToText','on',...
    'Interpreter','latex');
t.FontSize = 20;
xline(mean(ad),'Color','r','LineStyle',':','LineWidth',4);
xline(mean(hc),'Color','g','LineStyle',':','LineWidth',4);
if strcmp(group,'ad')
    legend({'AD','HC','',''});
else
    legend({'CTE','HC','',''});
end
set(gca,'FontSize',20)

%%% Region 3
subplot(2,3,3)
ad = hm.(reg3).(met).(group);
hc = hm.(reg3).(met).nc;
histogram(ad,'FaceColor','r','BinWidth',bw); hold on;
histogram(hc,'FaceColor','g','BinWidth',bw);
xlim([xlims(1),xlims(2)]); ylim([ylims1(1),ylims1(2)])
title(tstr{3},'FontSize',20);
xlabel(xlab); ylabel('Counts',"Rotation",90)
% Annotations
bhat = hm_stats.(reg3).(met).coef.(grp);
se = hm_stats.(reg3).(met).se.(grp);
t = hm_stats.(reg3).(met).t.(grp);
df = hm_stats.(reg3).(met).df.(grp);
pval = hm_stats.(reg3).(met).p.(grp);
txt = {[strcat('$\hat{\beta}$ = ',num2str(bhat,2))],...
        [strcat('se = ',num2str(se,2))],...
        [strcat('t = ',num2str(t,2))],...
        [strcat('df = ',num2str(df,2))],...
        [strcat('p = ',num2str(pval,2))]};
t=annotation('textbox',dim3,'String',txt,'FitBoxToText','on',...
    'Interpreter','latex');
t.FontSize = 20;
xline(mean(ad),'Color','r','LineStyle',':','LineWidth',4);
xline(mean(hc),'Color','g','LineStyle',':','LineWidth',4);
if strcmp(group,'ad')
    legend({'AD','HC','',''});
else
    legend({'CTE','HC','',''});
end
set(gca,'FontSize',20)

%%% Region 4
subplot(2,3,4)
ad = hm.(reg4).(met).(group);
hc = hm.(reg4).(met).nc;
histogram(ad,'FaceColor','r','BinWidth',bw); hold on;
histogram(hc,'FaceColor','g','BinWidth',bw);
xlim([xlims(1),xlims(2)]); ylim([ylims2(1),ylims2(2)])
title(tstr{4},'FontSize',20);
xlabel(xlab); ylabel('Counts',"Rotation",90)
% Annotations
bhat = hm_stats.(reg4).(met).coef.(grp);
se = hm_stats.(reg4).(met).se.(grp);
t = hm_stats.(reg4).(met).t.(grp);
df = hm_stats.(reg4).(met).df.(grp);
pval = hm_stats.(reg4).(met).p.(grp);
txt = {[strcat('$\hat{\beta}$ = ',num2str(bhat,2))],...
        [strcat('se = ',num2str(se,2))],...
        [strcat('t = ',num2str(t,2))],...
        [strcat('df = ',num2str(df,2))],...
        [strcat('p = ',num2str(pval,2))]};
t=annotation('textbox',dim4,'String',txt,'FitBoxToText','on',...
    'Interpreter','latex');
t.FontSize = 20;
xline(mean(ad),'Color','r','LineStyle',':','LineWidth',4);
xline(mean(hc),'Color','g','LineStyle',':','LineWidth',4);
if strcmp(group,'ad')
    legend({'AD','HC','',''});
else
    legend({'CTE','HC','',''});
end
set(gca,'FontSize',20)

%%% Region 5
subplot(2,3,5)
ad = hm.(reg5).(met).(group);
hc = hm.(reg5).(met).nc;
histogram(ad,'FaceColor','r','BinWidth',bw); hold on;
histogram(hc,'FaceColor','g','BinWidth',bw);
xlim([xlims(1),xlims(2)]); ylim([ylims2(1),ylims2(2)])
title(tstr{5},'FontSize',20);
xlabel(xlab); ylabel('Counts',"Rotation",90)
% Annotations
bhat = hm_stats.(reg5).(met).coef.(grp);
se = hm_stats.(reg5).(met).se.(grp);
t = hm_stats.(reg5).(met).t.(grp);
df = hm_stats.(reg5).(met).df.(grp);
pval = hm_stats.(reg5).(met).p.(grp);
txt = {[strcat('$\hat{\beta}$ = ',num2str(bhat,2))],...
        [strcat('se = ',num2str(se,2))],...
        [strcat('t = ',num2str(t,2))],...
        [strcat('df = ',num2str(df,2))],...
        [strcat('p = ',num2str(pval,2))]};
t=annotation('textbox',dim5,'String',txt,'FitBoxToText','on',...
    'Interpreter','latex');
t.FontSize = 20;
xline(mean(ad),'Color','r','LineStyle',':','LineWidth',4);
xline(mean(hc),'Color','g','LineStyle',':','LineWidth',4);
if strcmp(group,'ad')
    legend({'AD','HC','',''});
else
    legend({'CTE','HC','',''});
end
set(gca,'FontSize',20)

%%% Region 6
subplot(2,3,6)
ad = hm.(reg6).(met).(group);
hc = hm.(reg6).(met).nc;
histogram(ad,'FaceColor','r','BinWidth',bw); hold on;
histogram(hc,'FaceColor','g','BinWidth',bw);
xlim([xlims(1),xlims(2)]); ylim([ylims2(1),ylims2(2)])
title(tstr{6},'FontSize',20);
xlabel(xlab); ylabel('Counts',"Rotation",90)
% Annotations
bhat = hm_stats.(reg6).(met).coef.(grp);
se = hm_stats.(reg6).(met).se.(grp);
t = hm_stats.(reg6).(met).t.(grp);
df = hm_stats.(reg6).(met).df.(grp);
pval = hm_stats.(reg6).(met).p.(grp);
txt = {[strcat('$\hat{\beta}$ = ',num2str(bhat,2))],...
        [strcat('se = ',num2str(se,2))],...
        [strcat('t = ',num2str(t,2))],...
        [strcat('df = ',num2str(df,2))],...
        [strcat('p = ',num2str(pval,2))]};
t=annotation('textbox',dim6,'String',txt,'FitBoxToText','on',...
    'Interpreter','latex');
t.FontSize = 20;
xline(mean(ad),'Color','r','LineStyle',':','LineWidth',4);
xline(mean(hc),'Color','g','LineStyle',':','LineWidth',4);
if strcmp(group,'ad')
    legend({'AD','HC','',''});
else
    legend({'CTE','HC','',''});
end
set(gca,'FontSize',20)

%%% Save figure
fout = fullfile(mpath,fname);
saveas(gca,fout);
end

