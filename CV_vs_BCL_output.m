%Steffen Docken (Lewis Lab). All rights reserved.
%  Published in the Journal of Theoretical Biology under the title "Rate-dependent effects of state-specific sodium channel blockers in cardiac tissue: Insights from idealized models"
%8-21-17
%Code to visualize CV vs. BCL output

clear
close all

oldFolder = cd('1D_modified_Ten_Tusscher_with_drug_output/CV_results_12_22_17');

no_drug_dt001 = load('Simple_nodrug_CV_v_BCL_dt1eneg3_dx1eneg2_params_cell60.txt');

GRI_dt001 = load('Simple_HH_guarded_receptor_inactive_CV_v_BCL_dt1eneg3_dx1eneg2_kon_1e1_kd10e-6_Drug20e-6_params_cell60.txt');
GRN_dt001 = load('Simple_HH_guarded_receptor_noninactive_CV_v_BCL_dt1eneg3_dx1eneg2_kon_1e1_kd10e-6_Drug20e-6_params_cell60.txt');
GII_dt001 = load('Simple_HH_gate_immobilization_inactive_CV_v_BCL_dt1eneg3_dx1eneg2_kon_1e1_kd10e-6_Drug20e-6_params_cell60.txt');
GIN_dt001 = load('Simple_HH_gate_immobilization_noninactive_CV_v_BCL_dt1eneg3_dx1eneg2_kon_1e1_kd10e-6_Drug20e-6_params_cell60.txt');
%loading data for conduction velocity results

cd(oldFolder);

h = plot(1:10,1:10,1:10,2:11,1:10,3:12,1:10,4:13,1:10,5:14,1:10,6:15);
c = get(h,'Color'); %This graph is just so I can reference the colors in 
%future plots

UV_lim = [0, 140];
CV_lim = [0, 0.06];

x_tick_vals = [400, 600, 800, 1000];
y_tick_vals = [0, 0.02, 0.04, 0.06];
figure
p1 = plot(GRI_dt001(1:end-9,2), GRI_dt001(1:end-9,4), '-', GRN_dt001(:,2), GRN_dt001(:,4), '--',...
    GII_dt001(1:end-3,2), GII_dt001(1:end-3,4), '-',...
    GIN_dt001(:,2), GIN_dt001(:,4), '--', no_drug_dt001(:,2), no_drug_dt001(:,4), '-k');
ylim(CV_lim);
%legend('GRI', 'GRN', 'GII', 'GIN', 'no drug');
%title('CV 2cell vs. BCL (dt = 0.001, dx = 0.01)');
%ylabel('CV (cm/ms)');
%xlabel('BCL (ms)');

p1(1).Color = c{1};
p1(2).Color = c{1};
p1(3).Color = c{2};
p1(4).Color = c{2};
box off
set(gca,'XTick',x_tick_vals);
set(gca,'YTick',y_tick_vals);
%conduction velocity vs. BCL (GRI and GII data does not extend to the
%shortest BCL, because alternans developed at faster BCL in these models)
%note that these results correspond to calculating conduction velocity by 
%using the time it takes the wave to propagate the distance 2dx

figure
subplot(4,2,1)
plot(GRI_dt001(:,2), GRI_dt001(:,3), GRN_dt001(:,2), GRN_dt001(:,3), GII_dt001(:,2), GII_dt001(:,3),...
    GIN_dt001(:,2), GIN_dt001(:,3), no_drug_dt001(:,2), no_drug_dt001(:,3), 'k');
ylim(CV_lim);
title('CV 1cell vs. BCL'); %conduction velocity computed using the time it
%takes the wave to propagate the distance dx

subplot(4,2,2)
plot(GRI_dt001(:,2), GRI_dt001(:,4), GRN_dt001(:,2), GRN_dt001(:,4), GII_dt001(:,2), GII_dt001(:,4),...
    GIN_dt001(:,2), GIN_dt001(:,4), no_drug_dt001(:,2), no_drug_dt001(:,4), 'k');
ylim(CV_lim);
title('CV 2cell vs. BCL');%conduction velocity computed using the time it
%takes the wave to propagate the distance 2dx

subplot(4,2,3)
plot(GRI_dt001(:,2), GRI_dt001(:,5), GRN_dt001(:,2), GRN_dt001(:,5), GII_dt001(:,2), GII_dt001(:,5),...
    GIN_dt001(:,2), GIN_dt001(:,5), no_drug_dt001(:,2), no_drug_dt001(:,5), 'k');
ylim(CV_lim);
title('CV 10cell vs. BCL'); %conduction velocity computed using the time it
%takes the wave to propagate the distance 10dx


subplot(4,2,4)
plot(GRI_dt001(:,2), GRI_dt001(:,6), GRN_dt001(:,2), GRN_dt001(:,6), GII_dt001(:,2), GII_dt001(:,6),...
    GIN_dt001(:,2), GIN_dt001(:,6), no_drug_dt001(:,2), no_drug_dt001(:,6), 'k');
ylim(UV_lim);
title('UV vs. BCL');
%Upstroke velocity at the middle cell.

subplot(4,2,5)
plot(GRI_dt001(:,2), GRI_dt001(:,7), GRN_dt001(:,2), GRN_dt001(:,7), GII_dt001(:,2), GII_dt001(:,7),...
    GIN_dt001(:,2), GIN_dt001(:,7), no_drug_dt001(:,2), no_drug_dt001(:,7), 'k');
title('APD90 vs. BCL');
%APD90 at the middle cell

subplot(4,2,6)
plot(GRI_dt001(:,2), GRI_dt001(:,8), GRN_dt001(:,2), GRN_dt001(:,8), GII_dt001(:,2), GII_dt001(:,8),...
    GIN_dt001(:,2), GIN_dt001(:,8), no_drug_dt001(:,2), no_drug_dt001(:,8), 'k');
title('V_{thr} vs. BCL');
%Threshold potential at the middle cell (V of peak upstroke velocity).

subplot(4,2,7)
plot(GRI_dt001(:,2), GRI_dt001(:,9), GRN_dt001(:,2), GRN_dt001(:,9), GII_dt001(:,2), GII_dt001(:,9),...
    GIN_dt001(:,2), GIN_dt001(:,9), no_drug_dt001(:,2), no_drug_dt001(:,9), 'k');
title('Cycle Number vs. BCL');
%Number of action potentials that the middle cell has gone through
%(ensuring that ratio of stimuli to propagating waves is 1:1)


oldFolder = cd('1D_modified_Ten_Tusscher_with_drug_output/CV_results_12_22_17');

no_drug_dt001_480 = load('Simple_nodrug_CV_v_BCL_dt1eneg3_dx1eneg2_params_480s.txt');

GRI_dt001_480 = load('Simple_HH_guarded_receptor_inactive_CV_v_BCL_dt1eneg3_dx1eneg2_kon_1e1_kd10e-6_Drug20e-6_params_480s.txt');
GRN_dt001_480 = load('Simple_HH_guarded_receptor_noninactive_CV_v_BCL_dt1eneg3_dx1eneg2_kon_1e1_kd10e-6_Drug20e-6_params_480s.txt');
GII_dt001_480 = load('Simple_HH_gate_immobilization_inactive_CV_v_BCL_dt1eneg3_dx1eneg2_kon_1e1_kd10e-6_Drug20e-6_params_480s.txt');
GIN_dt001_480 = load('Simple_HH_gate_immobilization_noninactive_CV_v_BCL_dt1eneg3_dx1eneg2_kon_1e1_kd10e-6_Drug20e-6_params_480s.txt');


no_drug_dt001_500 = load('Simple_nodrug_CV_v_BCL_dt1eneg3_dx1eneg2_params_500s.txt');

GRI_dt001_500 = load('Simple_HH_guarded_receptor_inactive_CV_v_BCL_dt1eneg3_dx1eneg2_kon_1e1_kd10e-6_Drug20e-6_params_500s.txt');
GRN_dt001_500 = load('Simple_HH_guarded_receptor_noninactive_CV_v_BCL_dt1eneg3_dx1eneg2_kon_1e1_kd10e-6_Drug20e-6_params_500s.txt');
GII_dt001_500 = load('Simple_HH_gate_immobilization_inactive_CV_v_BCL_dt1eneg3_dx1eneg2_kon_1e1_kd10e-6_Drug20e-6_params_500s.txt');
GIN_dt001_500 = load('Simple_HH_gate_immobilization_noninactive_CV_v_BCL_dt1eneg3_dx1eneg2_kon_1e1_kd10e-6_Drug20e-6_params_500s.txt');

cd(oldFolder);


figure %for GRI, ensuring the system is at steady state by comparing values
%at each cell following the 480th and 500th stimuli
subplot(4,2,1)
plot(GRI_dt001_480(:,1), GRI_dt001_480(:,2), GRI_dt001_500(:,1), GRI_dt001_500(:,2));
ylim(CV_lim);
title('CV 1cell GRI');%conduction velocity computed using the time it
%takes the wave to propagate the distance dx

subplot(4,2,2)
plot(GRI_dt001_480(:,1), GRI_dt001_480(:,3), GRI_dt001_500(:,1), GRI_dt001_500(:,3));
ylim(CV_lim);
title('CV 2cell GRI');%conduction velocity computed using the time it
%takes the wave to propagate the distance 2dx

subplot(4,2,3)
plot(GRI_dt001_480(:,1), GRI_dt001_480(:,4), GRI_dt001_500(:,1), GRI_dt001_500(:,4));
ylim(CV_lim);
title('CV 10cell GRI');%conduction velocity computed using the time it
%takes the wave to propagate the distance 10dx

subplot(4,2,4)
plot(GRI_dt001_480(:,1), GRI_dt001_480(:,5), GRI_dt001_500(:,1), GRI_dt001_500(:,5));
ylim(UV_lim);
title('UV GRI');
%Upstroke velocity

subplot(4,2,5)
plot(GRI_dt001_480(:,1), GRI_dt001_480(:,6), GRI_dt001_500(:,1), GRI_dt001_500(:,6));
title('APD90 GRI');%APD90

subplot(4,2,6)
plot(GRI_dt001_480(:,1), GRI_dt001_480(:,7), GRI_dt001_500(:,1), GRI_dt001_500(:,7));
title('V_{thr} GRI');%Threshold potential (V of peak upstroke velocity).

subplot(4,2,7)
plot(GRI_dt001_480(:,1), GRI_dt001_480(:,8), GRI_dt001_500(:,1), GRI_dt001_500(:,8));
title('Cycle Number GRI');
%Number of action potentials that the middle cell has gone through
%(ensuring that ratio of stimuli to propagating waves is 1:1) (the cycle
%number drops by 1 around cell 400, because the last wave has not 
%propagated all the way through the cable yet)

%%
figure %for GRN, ensuring the system is at steady state by comparing values
%at each cell following the 480th and 500th stimuli
subplot(4,2,1)
plot(GRN_dt001_480(:,1), GRN_dt001_480(:,2), GRN_dt001_500(:,1), GRN_dt001_500(:,2));
ylim(CV_lim);
title('CV 1cell GRN');%conduction velocity computed using the time it
%takes the wave to propagate the distance dx

subplot(4,2,2)
plot(GRN_dt001_480(:,1), GRN_dt001_480(:,3), GRN_dt001_500(:,1), GRN_dt001_500(:,3));
ylim(CV_lim);
title('CV 2cell GRN');%conduction velocity computed using the time it
%takes the wave to propagate the distance 2dx

subplot(4,2,3)
plot(GRN_dt001_480(:,1), GRN_dt001_480(:,4), GRN_dt001_500(:,1), GRN_dt001_500(:,4));
ylim(CV_lim);
title('CV 10cell GRN');%conduction velocity computed using the time it
%takes the wave to propagate the distance 10dx


subplot(4,2,4)
plot(GRN_dt001_480(:,1), GRN_dt001_480(:,5), GRN_dt001_500(:,1), GRN_dt001_500(:,5));
ylim(UV_lim);
title('UV GRN');%Upstroke velocity

subplot(4,2,5)
plot(GRN_dt001_480(:,1), GRN_dt001_480(:,6), GRN_dt001_500(:,1), GRN_dt001_500(:,6));
title('APD90 GRN'); %APD90

subplot(4,2,6)
plot(GRN_dt001_480(:,1), GRN_dt001_480(:,7), GRN_dt001_500(:,1), GRN_dt001_500(:,7));
title('V_{thr} GRN');%Threshold potential (V of peak upstroke velocity).

subplot(4,2,7)
plot(GRN_dt001_480(:,1), GRN_dt001_480(:,8), GRN_dt001_500(:,1), GRN_dt001_500(:,8));
title('Cycle Number GRN');
%Number of action potentials that the middle cell has gone through
%(ensuring that ratio of stimuli to propagating waves is 1:1)

%% 
figure %for GII, ensuring the system is at steady state by comparing values
%at each cell following the 480th and 500th stimuli
subplot(4,2,1)
plot(GII_dt001_480(:,1), GII_dt001_480(:,2), GII_dt001_500(:,1), GII_dt001_500(:,2));
ylim(CV_lim);
title('CV 1cell GII');%conduction velocity computed using the time it
%takes the wave to propagate the distance dx

subplot(4,2,2)
plot(GII_dt001_480(:,1), GII_dt001_480(:,3), GII_dt001_500(:,1), GII_dt001_500(:,3));
ylim(CV_lim);
title('CV 2cell GII');%conduction velocity computed using the time it
%takes the wave to propagate the distance 2dx

subplot(4,2,3)
plot(GII_dt001_480(:,1), GII_dt001_480(:,4), GII_dt001_500(:,1), GII_dt001_500(:,4));
ylim(CV_lim);
title('CV 10cell GII');%conduction velocity computed using the time it
%takes the wave to propagate the distance 10dx

subplot(4,2,4)
plot(GII_dt001_480(:,1), GII_dt001_480(:,5), GII_dt001_500(:,1), GII_dt001_500(:,5));
ylim(UV_lim);
title('UV GII');%Upstroke velocity

subplot(4,2,5)
plot(GII_dt001_480(:,1), GII_dt001_480(:,6), GII_dt001_500(:,1), GII_dt001_500(:,6));

title('APD90 GII');%APD90

subplot(4,2,6)
plot(GII_dt001_480(:,1), GII_dt001_480(:,7), GII_dt001_500(:,1), GII_dt001_500(:,7));

title('V_{thr} GII');%Threshold potential (V of peak upstroke velocity).

subplot(4,2,7)
plot(GII_dt001_480(:,1), GII_dt001_480(:,8), GII_dt001_500(:,1), GII_dt001_500(:,8));
title('Cycle Number GII');%Number of action potentials that the middle cell has gone through
%(ensuring that ratio of stimuli to propagating waves is 1:1)

%%
figure%for GIN, ensuring the system is at steady state by comparing values
%at each cell following the 480th and 500th stimuli
subplot(4,2,1)
plot(GIN_dt001_480(:,1), GIN_dt001_480(:,2), GIN_dt001_500(:,1), GIN_dt001_500(:,2));
ylim(CV_lim);
title('CV 1cell GIN');%conduction velocity computed using the time it
%takes the wave to propagate the distance dx

subplot(4,2,2)
plot(GIN_dt001_480(:,1), GIN_dt001_480(:,3), GIN_dt001_500(:,1), GIN_dt001_500(:,3));
ylim(CV_lim);
title('CV 2cell GIN');%conduction velocity computed using the time it
%takes the wave to propagate the distance 2dx

subplot(4,2,3)
plot(GIN_dt001_480(:,1), GIN_dt001_480(:,4), GIN_dt001_500(:,1), GIN_dt001_500(:,4));
ylim(CV_lim);
title('CV 10cell GIN');%conduction velocity computed using the time it
%takes the wave to propagate the distance 10dx

subplot(4,2,4)
plot(GIN_dt001_480(:,1), GIN_dt001_480(:,5), GIN_dt001_500(:,1), GIN_dt001_500(:,5));
ylim(UV_lim);
title('UV GIN');%Upstroke velocity

subplot(4,2,5)
plot(GIN_dt001_480(:,1), GIN_dt001_480(:,6), GIN_dt001_500(:,1), GIN_dt001_500(:,6));
title('APD90 GIN');%APD90

subplot(4,2,6)
plot(GIN_dt001_480(:,1), GIN_dt001_480(:,7), GIN_dt001_500(:,1), GIN_dt001_500(:,7));
title('V_{thr} GIN');%Threshold potential (V of peak upstroke velocity).

subplot(4,2,7)
plot(GIN_dt001_480(:,1), GIN_dt001_480(:,8), GIN_dt001_500(:,1), GIN_dt001_500(:,8));
title('Cycle Number GIN');%Number of action potentials that the middle cell has gone through
%(ensuring that ratio of stimuli to propagating waves is 1:1)

%%
figure %for no drug, ensuring the system is at steady state by comparing values
%at each cell following the 480th and 500th stimuli
subplot(4,2,1)
plot(no_drug_dt001_480(:,1), no_drug_dt001_480(:,2), no_drug_dt001_500(:,1), no_drug_dt001_500(:,2));
ylim(CV_lim);
title('CV 1cell no drug');%conduction velocity computed using the time it
%takes the wave to propagate the distance dx

subplot(4,2,2)
plot(no_drug_dt001_480(:,1), no_drug_dt001_480(:,3), no_drug_dt001_500(:,1), no_drug_dt001_500(:,3));
ylim(CV_lim);
title('CV 2cell no drug');%conduction velocity computed using the time it
%takes the wave to propagate the distance 2dx

subplot(4,2,3)
plot(no_drug_dt001_480(:,1), no_drug_dt001_480(:,4), no_drug_dt001_500(:,1), no_drug_dt001_500(:,4));
ylim(CV_lim);
title('CV 10cell no drug');%conduction velocity computed using the time it
%takes the wave to propagate the distance 10dx

subplot(4,2,4)
plot(no_drug_dt001_480(:,1), no_drug_dt001_480(:,5), no_drug_dt001_500(:,1), no_drug_dt001_500(:,5));
ylim(UV_lim);
title('UV no drug');%Upstroke velocity

subplot(4,2,5)
plot(no_drug_dt001_480(:,1), no_drug_dt001_480(:,6), no_drug_dt001_500(:,1), no_drug_dt001_500(:,6));
title('APD90 no drug'); %APD90

subplot(4,2,6)
plot(no_drug_dt001_480(:,1), no_drug_dt001_480(:,7), no_drug_dt001_500(:,1), no_drug_dt001_500(:,7));
title('V_{thr} no drug'); %Threshold potential (V of peak upstroke velocity).

subplot(4,2,7)
plot(no_drug_dt001_480(:,1), no_drug_dt001_480(:,8), no_drug_dt001_500(:,1), no_drug_dt001_500(:,8));
title('Cycle Number no drug'); %Number of action potentials that the middle cell has gone through
%(ensuring that ratio of stimuli to propagating waves is 1:1)

%Calculating the maximum relative changes in cell parameters to ensure that
%steady state was reached
no_drug_max_rel_change = max(max(abs(no_drug_dt001_480(10:end, 2:end-1) - no_drug_dt001_500(10:end, 2:end-1))./...
    no_drug_dt001_500(10:end, 2:end-1)))

GRI_max_rel_change = max(max(abs(GRI_dt001_480(10:end, 2:end-1) - GRI_dt001_500(10:end, 2:end-1))./...
    GRI_dt001_500(10:end, 2:end-1)))

GRN_max_rel_change = max(max(abs(GRN_dt001_480(10:end, 2:end-1) - GRN_dt001_500(10:end, 2:end-1))./...
    GRN_dt001_500(10:end, 2:end-1)))

GII_max_rel_change = max(max(abs(GII_dt001_480(10:end, 2:end-1) - GII_dt001_500(10:end, 2:end-1))./...
    GII_dt001_500(10:end, 2:end-1)))

GIN_max_rel_change = max(max(abs(GIN_dt001_480(10:end, 2:end-1) - GIN_dt001_500(10:end, 2:end-1))./...
    GIN_dt001_500(10:end, 2:end-1))) %relative changes in values (other 
%than cycle number) between 480s and 500s 
%to check if at steady state.  Ignoring the first 10 cells, since these may
%have boundry effects.

no_drug_max_rel_change_60 = max(max(abs(no_drug_dt001_480(60, 2:end-1) - no_drug_dt001_500(60, 2:end-1))./...
    no_drug_dt001_500(60, 2:end-1)))

GRI_max_rel_change_60 = max(max(abs(GRI_dt001_480(60, 2:end-1) - GRI_dt001_500(60, 2:end-1))./...
    GRI_dt001_500(60, 2:end-1)))

GRN_max_rel_change_60 = max(max(abs(GRN_dt001_480(60, 2:end-1) - GRN_dt001_500(60, 2:end-1))./...
    GRN_dt001_500(60, 2:end-1)))

GII_max_rel_change_60 = max(max(abs(GII_dt001_480(60, 2:end-1) - GII_dt001_500(60, 2:end-1))./...
    GII_dt001_500(60, 2:end-1)))

GIN_max_rel_change_60 = max(max(abs(GIN_dt001_480(60, 2:end-1) - GIN_dt001_500(60, 2:end-1))./...
    GIN_dt001_500(60, 2:end-1)))


figure
p1 = plot(GRI_dt001(:,2), GRI_dt001(:,6), '-', GRN_dt001(:,2), GRN_dt001(:,6), '--',...
    GII_dt001(:,2), GII_dt001(:,6), '-', ...
    GIN_dt001(:,2), GIN_dt001(:,6), '--', no_drug_dt001(:,2), no_drug_dt001(:,6), 'k');
ylim(UV_lim);
title('UV vs. BCL 1D'); %Upstroke velocity vs. BCL

p1(1).Color = c{1};
p1(2).Color = c{1};
p1(3).Color = c{2};
p1(4).Color = c{2};

%% For Graphical Abstract
load Upstroke_v_Period_3_tau_b.mat

figure
yyaxis left
p1 = plot(GRFI_1e1(7:end,1), GRFI_1e1(7:end,3), '-',...
    GIFI_1e1(7:end,1), GIFI_1e1(7:end,3), '-', ...
    nodrug(7:end,1), nodrug(7:end,3), '-k');
ylim([100, 310]);
ylabel('Peak Upstroke Velocity (mV/ms)')
xlabel('Basic Cycle Length (ms)')
p1(1).Color = c{1};
p1(2).Color = c{2};
set(gca,'XTick',400:200:1000);
set(gca,'YTick', 100:100:300);

yyaxis right
p2 = plot(GRI_dt001(1:end-9,2), GRI_dt001(1:end-9,4), '-.',...
    GII_dt001(1:end-3,2), GII_dt001(1:end-3,4), '-.',...
    no_drug_dt001(:,2), no_drug_dt001(:,4), '-.k');
ylim([0,0.062]);
ylabel('Conduction Velocity (cm/ms)')
p2(1).Color = c{1};
p2(2).Color = c{2};

box off
set(gca,'YTick',y_tick_vals, 'FontName', 'Arial', ...
        'FontSize', 5, ...
        'LineWidth', 0.5);

