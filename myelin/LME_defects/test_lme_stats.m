%% Test the function "lme_stats" using data from Anna
% This function was writting for testing Anna's experimental data. This
% test script was written to test the function.
clear; clc; close all;


%% Import the spreadsheet of data

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
var = table2array(T(:,3));

grp = zeros(size(subid));
pmi = zeros(size(subid));

for i=1:length(subid)
    group_idx = find(id_table == subid(i));
    pmi(i) = pmi_table(group_idx);

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
grp1_var = var(grp1,:);
grp2_var = var(grp2,:);
grp3_var = var(grp3,:);

% Extract PMI for each group
grp1_pmi = pmi(grp1,:);
grp2_pmi = pmi(grp2,:);
grp3_pmi = pmi(grp3,:);

%% Call lme_stats function on a subset of data
%function [p, coef, lower, upper] =...
    %lme_stats(exp_values, exp_sub_id, cntrl_values, cntrl_sub_id)
% Compare group 2 and 3
[p, coef, lower, upper] =...
    lme_stats(grp3_var,grp3_subid, grp3_pmi, grp2_var,grp2_subid, grp2_pmi);

