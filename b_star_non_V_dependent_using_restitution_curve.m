%Steffen Docken (Lewis Lab). All rights reserved.
%  Published in the Journal of Theoretical Biology under the title "Rate-dependent effects of state-specific sodium channel blockers in cardiac tissue: Insights from idealized models"
%9-27-17
%This code calculates the b* values for drug models that are non-V-dependent
%(i.e., neutral drugs) given the square wave approximation of
%a train of APs and using the APD and DI values from the restitution curves
%for each of the drug binding models

clear 
close all


Drug = 20e-6; %Drug concentration in M

oldFolder = cd('Ten_Tusscher_model_with_drug_output/Upstroke_vs_Period_nonVdepend_2_28_18/');

nodrug = load('Simple_nodrug_drug_Upstroke_v_Period_stimdur06.txt');

GIN_nonVdepend_1e1 = load('Simple_HH_gate_immobilization_noninactive_drug_Upstroke_v_Period_Diffval1e1_kd10e-6_Drug20e-6_stimdur06.txt');

GII_nonVdepend_1e1 = load('Simple_HH_gate_immobilization_inactive_drug_Upstroke_v_Period_Diffval1e1_kd10e-6_Drug20e-6_stimdur06.txt');

GRN_nonVdepend_1e1 = load('Simple_HH_guarded_receptor_noninactive_drug_Upstroke_v_Period_Diffval1e1_kd10e-6_Drug20e-6_stimdur06.txt');

GRI_nonVdepend_1e1 = load('Simple_HH_guarded_receptor_inactive_drug_Upstroke_v_Period_Diffval1e1_kd10e-6_Drug20e-6_stimdur06.txt');

cd(oldFolder);%loading results for drug binding models where drug binding 
%is non-V-dependent


V_DI = -80;
V_A = 20; %in mV

R = 8314.472;	% J/mol*K
T = 310;% K
F = 96485.3415; 

Diffusion = 1e1; %in M^-1ms^-1 

k_D0 = 10e-6; %in M

inact_noninact_vec = [0, 1]; %inactive state binding

transparency_param = .5;

N = length(nodrug(:,1)); %number of BCLs used

GR_restitution = zeros(N,4,2);
GI_restitution = zeros(N,4,2);

GR_restitution(:,:,1) = GRI_nonVdepend_1e1; %inact binding and Diff =1e1
GR_restitution(:,:,2) = GRN_nonVdepend_1e1; %non-inact binding and Diff =1e1

GI_restitution(:,:,1) = GII_nonVdepend_1e1; %inact binding and Diff =1e1
GI_restitution(:,:,2) = GIN_nonVdepend_1e1; %non-inact binding and Diff =1e1


b_star_gr_mat = zeros(N, 2);
b_star_gi_mat = zeros(N, 2); %1st column will correspond to inact, 2nd 
%to non-inact binding

k_on = Diffusion;
k_off_0 = Diffusion*k_D0; %defining params from options

for ll = 1:2 %looping through state binding
    inact_noninact = inact_noninact_vec(ll);

    %% guarded receptor binding
    gi_gr_nons = 1; 
    %HH formulation

    for kk = 1:N
        [~,b_star_gr_mat(kk,ll)] = HH_b_star_neutral_drug(GR_restitution(kk,4,ll),...
            GR_restitution(kk,1,ll)-GR_restitution(kk,4,ll), V_A,...
            V_DI, k_on, k_off_0, Drug, inact_noninact, ...
            gi_gr_nons);  %calculating b_star for each BCL
    end

    %% gate immobilization binding
    gi_gr_nons = 0; 
    %HH formulation

    for kk = 1:N
        [~,b_star_gi_mat(kk,ll)] = HH_b_star_neutral_drug(GI_restitution(kk,4,ll),...
            GI_restitution(kk,1,ll)-GI_restitution(kk,4,ll), V_A,...
            V_DI, k_on, k_off_0, Drug, inact_noninact, ...
            gi_gr_nons); %calculating b_star for each BCL
    end

end


xlim_vec = [300, 1000];
x_tick_vals = [400, 600, 800, 1000];

ylim_vec = [0,1];
ylim_vec_nonnorm = [100, 350];
y_tick_vals = [0, 0.25, 0.5, 0.75, 1];
y_tick_vals_nonnorm = [100, 200, 300];

h = plot(1:10,1:10,1:10,2:11,1:10,3:12,1:10,4:13,1:10,5:14,1:10,6:15);
c = get(h,'Color');


figure
p1 = plot(nodrug(7:end,1), 1-b_star_gr_mat(7:end,1), '-', nodrug(7:end,1), 1-b_star_gr_mat(7:end,2)+.01, '--',...
    nodrug(7:end,1), 1-b_star_gi_mat(7:end,1), '-', nodrug(7:end,1), 1-b_star_gi_mat(7:end,2), '--', ...
    nodrug(7:end,1), GR_restitution(7:end,3,1)./nodrug(7:end,3), '-',...
    nodrug(7:end,1), GR_restitution(7:end,3,2)./nodrug(7:end,3)+.01, '--',...
    nodrug(7:end,1), GI_restitution(7:end,3,1)./nodrug(7:end,3), '-',...
    nodrug(7:end,1), GI_restitution(7:end,3,2)./nodrug(7:end,3), '--',...
    'MarkerSize', 3);
%title('k_{on} = 1e1')
xlim(xlim_vec);
ylim(ylim_vec);
p1(1).Color = [c{1}, transparency_param];
p1(2).Color = [c{1}, transparency_param];
p1(3).Color = [c{2}, transparency_param];
p1(4).Color = [c{2}, transparency_param];
p1(5).Color = c{1};
p1(6).Color = c{1};
p1(7).Color = c{2};
p1(8).Color = c{2};
box off;
set(gca,'XTick',x_tick_vals);
set(gca,'YTick',y_tick_vals, 'YTickLabel', []);
yyaxis right
ax = gca;
ax.YColor = transparency_param*[1, 1, 1];
ax.YTick = y_tick_vals;
ax.YTickLabel = []; %comparing the steady state fraction of
%channels not bound to drug during the upstroke using the 1D map to the 
%peak upstroke velocity of each
%drug model normalized to the peak upstroke velocity in the no drug model.
