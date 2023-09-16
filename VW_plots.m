%Steffen Docken (Lewis Lab). All rights reserved.
%  Published in the Journal of Theoretical Biology under the title "Rate-dependent effects of state-specific sodium channel blockers in cardiac tissue: Insights from idealized models"
%9-30-17
%Code to create plots of VW

clear
close all

oldFolder = cd('1D_modified_Ten_Tusscher_with_drug_output/VW_results_12_22_17');


%no drug
no_drug_400ms_VW_times = load('Simple_nodrug_VW_find_VW_BCL_400ms_dt1eneg3_dx1eneg2_VW_times.txt');
no_drug_675ms_VW_times = load('Simple_nodrug_VW_find_VW_BCL_675ms_dt1eneg3_dx1eneg2_VW_times.txt');
no_drug_1000ms_VW_times = load('Simple_nodrug_VW_find_VW_BCL_1000ms_dt1eneg3_dx1eneg2_VW_times.txt');

% GRI
GRI_675ms_VW_times = load('Simple_HH_guarded_receptor_inactive_VW_find_VW_BCL_675ms_dt1eneg3_dx1eneg2_kon_1e1_kd10e-6_Drug20e-6_VW_times.txt');
GRI_1000ms_VW_times = load('Simple_HH_guarded_receptor_inactive_VW_find_VW_BCL_1000ms_dt1eneg3_dx1eneg2_kon_1e1_kd10e-6_Drug20e-6_VW_times.txt');

% GRN
GRN_400ms_VW_times = load('Simple_HH_guarded_receptor_noninactive_VW_find_VW_BCL_400ms_dt1eneg3_dx1eneg2_kon_1e1_kd10e-6_Drug20e-6_VW_times.txt');
GRN_675ms_VW_times = load('Simple_HH_guarded_receptor_noninactive_VW_find_VW_BCL_675ms_dt1eneg3_dx1eneg2_kon_1e1_kd10e-6_Drug20e-6_VW_times.txt');
GRN_1000ms_VW_times = load('Simple_HH_guarded_receptor_noninactive_VW_find_VW_BCL_1000ms_dt1eneg3_dx1eneg2_kon_1e1_kd10e-6_Drug20e-6_VW_times.txt');

%GII
GII_400ms_VW_times = load('Simple_HH_gate_immobilization_inactive_VW_find_VW_BCL_400ms_dt1eneg3_dx1eneg2_kon_1e1_kd10e-6_Drug20e-6_VW_times.txt');
GII_675ms_VW_times = load('Simple_HH_gate_immobilization_inactive_VW_find_VW_BCL_675ms_dt1eneg3_dx1eneg2_kon_1e1_kd10e-6_Drug20e-6_VW_times.txt');
GII_1000ms_VW_times = load('Simple_HH_gate_immobilization_inactive_VW_find_VW_BCL_1000ms_dt1eneg3_dx1eneg2_kon_1e1_kd10e-6_Drug20e-6_VW_times.txt');

%GIN
GIN_400ms_VW_times = load('Simple_HH_gate_immobilization_noninactive_VW_find_VW_BCL_400ms_dt1eneg3_dx1eneg2_kon_1e1_kd10e-6_Drug20e-6_VW_times.txt');
GIN_675ms_VW_times = load('Simple_HH_gate_immobilization_noninactive_VW_find_VW_BCL_675ms_dt1eneg3_dx1eneg2_kon_1e1_kd10e-6_Drug20e-6_VW_times.txt');
GIN_1000ms_VW_times = load('Simple_HH_gate_immobilization_noninactive_VW_find_VW_BCL_1000ms_dt1eneg3_dx1eneg2_kon_1e1_kd10e-6_Drug20e-6_VW_times.txt');

%no drug movie output
no_drug_1000ms_1way_prop_previous_Vs = load('Simple_nodrug_VW_find_VW_BCL_1000ms_dt1eneg3_dx1eneg2_1way_prop_previous_Vs.txt');
no_drug_1000ms_2way_prop_previous_Vs = load('Simple_nodrug_VW_find_VW_BCL_1000ms_dt1eneg3_dx1eneg2_2way_prop_previous_Vs.txt');

no_drug_1000ms_1way_prop_Vs = load('Simple_nodrug_VW_find_VW_BCL_1000ms_dt1eneg3_dx1eneg2_1way_prop_Vs.txt');
no_drug_1000ms_2way_prop_Vs = load('Simple_nodrug_VW_find_VW_BCL_1000ms_dt1eneg3_dx1eneg2_2way_prop_Vs.txt');

cd(oldFolder);

oldFolder = cd('1D_modified_Ten_Tusscher_with_drug_output/VW_initialization_12_15_17');

no_drug_1000ms_last2AP = load('Simple_nodrug_VW_initialization_BCL_1000ms_dt1eneg3_dx1eneg2_params_cell260.txt');
GRI_1000ms_last2AP = load('Simple_HH_guarded_receptor_inactive_VW_initialization_BCL_1000ms_dt1eneg3_dx1eneg2_kon_1e1_kd10e-6_Drug20e-6_params_cell260.txt');
GRN_1000ms_last2AP = load('Simple_HH_guarded_receptor_noninactive_VW_initialization_BCL_1000ms_dt1eneg3_dx1eneg2_kon_1e1_kd10e-6_Drug20e-6_params_cell260.txt');
GII_1000ms_last2AP = load('Simple_HH_gate_immobilization_inactive_VW_initialization_BCL_1000ms_dt1eneg3_dx1eneg2_kon_1e1_kd10e-6_Drug20e-6_params_cell260.txt');
GIN_1000ms_last2AP = load('Simple_HH_gate_immobilization_noninactive_VW_initialization_BCL_1000ms_dt1eneg3_dx1eneg2_kon_1e1_kd10e-6_Drug20e-6_params_cell260.txt');


no_drug_675ms_last2AP = load('Simple_nodrug_VW_initialization_BCL_675ms_dt1eneg3_dx1eneg2_params_cell260.txt');
GRI_675ms_last2AP = load('Simple_HH_guarded_receptor_inactive_VW_initialization_BCL_675ms_dt1eneg3_dx1eneg2_kon_1e1_kd10e-6_Drug20e-6_params_cell260.txt');
GRN_675ms_last2AP = load('Simple_HH_guarded_receptor_noninactive_VW_initialization_BCL_675ms_dt1eneg3_dx1eneg2_kon_1e1_kd10e-6_Drug20e-6_params_cell260.txt');
GII_675ms_last2AP = load('Simple_HH_gate_immobilization_inactive_VW_initialization_BCL_675ms_dt1eneg3_dx1eneg2_kon_1e1_kd10e-6_Drug20e-6_params_cell260.txt');
GIN_675ms_last2AP = load('Simple_HH_gate_immobilization_noninactive_VW_initialization_BCL_675ms_dt1eneg3_dx1eneg2_kon_1e1_kd10e-6_Drug20e-6_params_cell260.txt');

GRI_675ms_500 = load('Simple_HH_guarded_receptor_inactive_VW_initialization_BCL_675ms_dt1eneg3_dx1eneg2_kon_1e1_kd10e-6_Drug20e-6_params_500.txt');


no_drug_400ms_last2AP = load('Simple_nodrug_VW_initialization_BCL_400ms_dt1eneg3_dx1eneg2_params_cell260.txt');
%GRI_400ms_last2AP =
%load('Simple_HH_guarded_receptor_inactive_VW_initialization_BCL_400ms_dt1eneg3_dx1eneg2_kon_1e1_kd10e-6_Drug20e-6_params_cell260.txt');
%This file doesn't exist, because a 1:1 ratio of stimulus to wave does not
%exist in this model at BCL = 400ms.
GRN_400ms_last2AP = load('Simple_HH_guarded_receptor_noninactive_VW_initialization_BCL_400ms_dt1eneg3_dx1eneg2_kon_1e1_kd10e-6_Drug20e-6_params_cell260.txt');
GII_400ms_last2AP = load('Simple_HH_gate_immobilization_inactive_VW_initialization_BCL_400ms_dt1eneg3_dx1eneg2_kon_1e1_kd10e-6_Drug20e-6_params_cell260.txt');
GIN_400ms_last2AP = load('Simple_HH_gate_immobilization_noninactive_VW_initialization_BCL_400ms_dt1eneg3_dx1eneg2_kon_1e1_kd10e-6_Drug20e-6_params_cell260.txt');

cd(oldFolder);

h = plot(1:10,1:10,1:10,2:11,1:10,3:12,1:10,4:13,1:10,5:14,1:10,6:15);
c = get(h,'Color'); %This graph is just so I can reference the colors in 
%future plots

BCL_vec = [400, 675, 1000];

no_drug_VW = [no_drug_400ms_VW_times, no_drug_675ms_VW_times, ...
    no_drug_1000ms_VW_times];

BCL_vec_GRI = [675, 1000];
GRI_VW = [GRI_675ms_VW_times, GRI_1000ms_VW_times];

GRN_VW = [GRN_400ms_VW_times, GRN_675ms_VW_times, GRN_1000ms_VW_times];

GII_VW = [GII_400ms_VW_times, GII_675ms_VW_times, GII_1000ms_VW_times];

GIN_VW = [GIN_400ms_VW_times, GIN_675ms_VW_times, GIN_1000ms_VW_times];
%combining results into arrays

figure(1)
x_tick_vals = [400, 600, 800, 1000];
y_tick_vals = [0, .04, .08, .12];
y_tick_labels = {'0%', '4%', '8%', '12%'};
p1 = plot(BCL_vec_GRI, GRI_VW(2,:) - GRI_VW(1,:), '-o', BCL_vec, GRN_VW(2,:) - ...
    GRN_VW(1,:), '--o', BCL_vec, GII_VW(2,:) - GII_VW(1,:), '-o', BCL_vec, ...
    GIN_VW(2,:) - GIN_VW(1,:), '--o', BCL_vec, no_drug_VW(2,:) - no_drug_VW(1,:), '-ok');
title('Length of VW vs. BCL');
ylabel('Length of VW (ms)');
xlabel('BCL (ms)');
legend('GRI', 'GRN', 'GII', 'GIN', 'no drug');
p1(1).Color = c{1};
p1(2).Color = c{1};
p1(3).Color = c{2};
p1(4).Color = c{2};
%box off
set(gca,'XTick',x_tick_vals);

figure(2)
p1 = plot(BCL_vec_GRI, (GRI_VW(2,:) - GRI_VW(1,:))./BCL_vec_GRI, '-o' ,BCL_vec, (GRN_VW(2,:) - ...
    GRN_VW(1,:))./BCL_vec, '--o', BCL_vec, (GII_VW(2,:) - GII_VW(1,:))./BCL_vec, '-o', BCL_vec, ...
    (GIN_VW(2,:) - GIN_VW(1,:))./BCL_vec, '--o', BCL_vec,...
    (no_drug_VW(2,:) - no_drug_VW(1,:))./BCL_vec, '-ok');
%title('fraction of BCL that is VW vs. BCL');
% ylabel('fraction of BCL');
% xlabel('BCL (ms)');
% legend('GRI', 'GRN', 'GII', 'GIN', 'no drug');
box off
xlim([350, 1000]);
p1(1).Color = c{1};
p1(2).Color = c{1};
p1(3).Color = c{2};
p1(4).Color = c{2};
set(gca,'XTick',x_tick_vals);
set(gca,'YTick',y_tick_vals, 'YTickLabel', y_tick_labels); %fraction of the BCL that is composed of the
% vulnerable window

% figure(3)
% plot(BCL_vec_GRI, GRI_VW(1,:), BCL_vec, ...
%     GRN_VW(1,:), BCL_vec, GII_VW(1,:), BCL_vec, ...
%     GIN_VW(1,:), BCL_vec, no_drug_VW(1,:), 'k');
% ylabel('S2 (ms)');
% xlabel('BCL (ms)');
% legend('GRI', 'GRN', 'GII', 'GIN', 'no drug');
% hold on
% p = plot(BCL_vec_GRI, GRI_VW(2,:), BCL_vec, ...
%     GRN_VW(2,:), BCL_vec, GII_VW(2,:), BCL_vec, ...
%     GIN_VW(2,:), BCL_vec, no_drug_VW(2,:), 'k');
% hold off
% p(1).Color = c{1};
% p(2).Color = c{2};
% p(3).Color = c{3};
% p(4).Color = c{4}; %timing of smallest S2 of 1way prop vs. BCL and timing
% % of smallest S2 of 2way prop vs. BCL
% 
% figure(4)
% plot(BCL_vec_GRI, GRI_VW(1,:), BCL_vec, ...
%     GRN_VW(1,:), BCL_vec, GII_VW(1,:), BCL_vec, ...
%     GIN_VW(1,:), BCL_vec, no_drug_VW(1,:), 'k');
% ylabel('S2 (ms)');
% xlabel('BCL (ms)');
% legend('GRI', 'GRN', 'GII', 'GIN', 'no drug');
% title('smallest S2 of 1way prop vs. BCL');
% 
% figure(5)
% plot(BCL_vec_GRI, GRI_VW(2,:), BCL_vec, ...
%     GRN_VW(2,:), BCL_vec, GII_VW(2,:), BCL_vec, ...
%     GIN_VW(2,:), BCL_vec, no_drug_VW(2,:), 'k');
% ylabel('S2 (ms)');
% xlabel('BCL (ms)');
% legend('GRI', 'GRN', 'GII', 'GIN', 'no drug');
% title('smallest S2 of 2way prop vs. BCL');

end_prop_index = 351;

view_vec = [0, 70];
x_vec = 1:1:510;
figure(6)
AxesH = axes;
surf(no_drug_1000ms_1way_prop_previous_Vs(1:end_prop_index,1), x_vec,...
    no_drug_1000ms_1way_prop_previous_Vs(1:end_prop_index,end:-1:2)', 'Linestyle', ...
    'none');
colormap('gray')
view(view_vec)
axis off
InSet = get(AxesH, 'TightInset');
set(AxesH, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)])
%visualization of a propagating wave followed by a stimulus that elicits no
%propagating wave

figure(7)
AxesH = axes;
surf(no_drug_1000ms_2way_prop_previous_Vs(1:end_prop_index,1), x_vec,...
    no_drug_1000ms_2way_prop_previous_Vs(1:end_prop_index,end:-1:2)', 'Linestyle', ...
    'none');
colormap('gray')
view(view_vec)
axis off
InSet = get(AxesH, 'TightInset');
set(AxesH, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)])
%visualization of a propagating wave followed by a stimulus that elicits 
%a unidirectional propagating wave

figure(8)
AxesH = axes;
surf(no_drug_1000ms_2way_prop_Vs(1:end_prop_index,1), x_vec,...
    no_drug_1000ms_2way_prop_Vs(1:end_prop_index,end:-1:2)', 'Linestyle', ...
    'none');
colormap('gray')
view(view_vec)
axis off
InSet = get(AxesH, 'TightInset');
set(AxesH, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)])
%visualization of a propagating wave followed by a stimulus that elicits
%two new propagating waves

