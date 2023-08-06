%Steffen Docken (Lewis Lab). All rights reserved.
%  Published in the Journal of Theoretical Biology under the title "Rate-dependent effects of state-specific sodium channel blockers in cardiac tissue: Insights from idealized models"
%8-24-16
%Code to compute the dynamics of h and b throughout a square wave V-clamp
%protocol

V_DI = -80;
V_A = 20; %in mV

Temp = 310;

a_h = @(V)a_h_6_14_2016(V, Temp);
b_h = @(V)b_h_6_14_2016(V, Temp);


Drug = 20e-6; %Drug concentration in M

Diffusion = 1e2; %in M^-1ms^-1 first value is approximately that of 
%Flecainide and second is Lidocaine from the Moreno et al paper.

k_D0 = 10e-6; %in ms^-1 first value is approximately that of 
%Flecainide and second is Lidocaine from the Moreno et al paper.

%% getting APD and DI
load Upstroke_v_Period_3_tau_b.mat
BCL_750_ind = find(nodrug(:,1) == 750, 1);
APD = round(nodrug(BCL_750_ind,4));
DI = 750-APD;

k_on = Diffusion;
k_off_0 = Diffusion*k_D0; %defining drug binding rates


inact_noninact = [0, 1]; %inact (0) or noninact (1) state binding

dt = 0.1;
t_vec = 0:dt:2250; %time vec in ms

h = plot(1:10,1:10,1:10,2:11,1:10,3:12,1:10,4:13,1:10,5:14,1:10,6:15);
c = get(h,'Color');

transparency_param = .5;

b_nons_vec = zeros(length(t_vec),2); %vector for b throughout 3 Duty cycles
h_nons_vec = zeros(length(t_vec),2); %vector for h throughout 3 Duty cycles

b_gr_vec = zeros(length(t_vec),2); %vector for b throughout 3 Duty cycles
h_gr_vec = zeros(length(t_vec),2); %vector for h throughout 3 Duty cycles

b_gi_vec = zeros(length(t_vec),2); %vector for b throughout 3 Duty cycles
h_gi_vec = zeros(length(t_vec),2); %vector for h throughout 3 Duty cycles

hb_gi_vec = zeros(2, length(t_vec),2); %array for b and h throughout 3 Duty cycles

V_vec = zeros(1, length(t_vec)); %Vector to hold V values throughout time
V_vec(1) = V_A;

for ll = 1:2
    %non-state-specific plots
    GIHH_GRHH_nons = 2; %non-state-specific binding
    
    [h_ss_nons_init, b_ss_nons_init] = HH_b_star(APD, DI, V_A, V_DI, k_on, k_off_0, Drug, ...
    inact_noninact(ll), GIHH_GRHH_nons);%calculating the fraction 
%of channels bound to drug at the end of an AP

    b_nons_vec(1, ll) = b_ss_nons_init;
    h_nons_vec(1, ll) = h_ss_nons_init;

    for ii = 0:2
        for jj = 1:APD/dt
            [h_nons_vec(ii*(APD+DI)/dt + 1 + jj,ll), b_nons_vec(ii*(APD+DI)/dt + 1 + jj,ll)] = ...
                HH_analytic(V_A, dt, a_h(V_A)/(a_h(V_A) + b_h(V_A)),...
                b_nons_vec(ii*(APD+DI)/dt + jj,ll), k_on, k_off_0, Drug,...
                inact_noninact(ll), GIHH_GRHH_nons);
        end %updating through the AP

        for jj = 1:DI/dt
            [h_nons_vec((ii*(APD+DI) + APD)/dt + 1 + jj,ll), b_nons_vec((ii*(APD+DI) ...
                + APD)/dt + 1 + jj,ll)] = ...
                HH_analytic(V_DI, dt, a_h(V_DI)/(a_h(V_DI) + b_h(V_DI)),...
                b_nons_vec((ii*(APD+DI) + APD)/dt + jj,ll), k_on, k_off_0, Drug,...
                inact_noninact(ll), GIHH_GRHH_nons);
        end %updating through the DI
    end

    %guarded receptor plots

    GIHH_GRHH_nons = 1; %guarded Receptor model with the HH formulation
    %GRHH_nons = 1; %non-state-specific binding

    [h_ss_gr_init, b_ss_gr_init] = HH_b_star(APD, DI, V_A, V_DI, k_on, k_off_0, Drug, ...
        inact_noninact(ll), GIHH_GRHH_nons);%calculating the fraction 
%of channels bound to drug at the end of an AP

    b_gr_vec(1,ll) = b_ss_gr_init;
    h_gr_vec(1,ll) = h_ss_gr_init;

    for ii = 0:2
        for jj = 1:APD/dt
            [h_gr_vec(ii*(APD+DI)/dt + 1 + jj,ll), b_gr_vec(ii*(APD+DI)/dt + 1 + jj,ll)] = ...
                HH_analytic(V_A, dt, a_h(V_A)/(a_h(V_A) + b_h(V_A)),...
                b_gr_vec(ii*(APD+DI)/dt + jj,ll), k_on, k_off_0, Drug,...
                inact_noninact(ll), GIHH_GRHH_nons);
        end %updating through AP

        for jj = 1:DI/dt
            [h_gr_vec((ii*(APD+DI) + APD)/dt + 1 + jj,ll), b_gr_vec((ii*(APD+DI) ...
                + APD)/dt + 1 + jj,ll)] = ...
                HH_analytic(V_DI, dt, a_h(V_DI)/(a_h(V_DI) + b_h(V_DI)),...
                b_gr_vec((ii*(APD+DI) ...
                + APD)/dt + jj,ll), k_on, k_off_0, Drug,...
                inact_noninact(ll), GIHH_GRHH_nons);
        end %updating through DI
    end
    
    %Gate Immobilization plots

    GIHH_GRHH_nons = 0; %guarded Receptor model with the HH formulation
    %GRHH_nons = 1; %non-state-specific binding

    [h_ss_gi_init, b_ss_gi_init] = HH_b_star(APD, DI, V_A, V_DI, k_on, k_off_0, Drug, ...
        inact_noninact(ll), GIHH_GRHH_nons);%calculating the fraction 
%of channels bound to drug at the end of an AP

    b_gi_vec(1,ll) = b_ss_gi_init;
    h_gi_vec(1,ll) = h_ss_gi_init;

    for ii = 0:2
        for jj = 1:APD/dt
            [h_gi_vec(ii*(APD+DI)/dt + 1 + jj,ll), b_gi_vec(ii*(APD+DI)/dt + 1 + jj,ll)] = ...
                HH_analytic(V_A, dt, a_h(V_A)/(a_h(V_A) + b_h(V_A)),...
                b_gi_vec(ii*(APD+DI)/dt + jj,ll), k_on, k_off_0, Drug,...
                inact_noninact(ll), GIHH_GRHH_nons);
            
            V_vec(ii*(APD+DI)/dt + 1 + jj) = V_A; %storing V
        end %updating through AP

        for jj = 1:DI/dt
            [h_gi_vec((ii*(APD+DI) + APD)/dt + 1 + jj,ll), b_gi_vec((ii*(APD+DI) ...
                + APD)/dt + 1 + jj,ll)] = ...
                HH_analytic(V_DI, dt, a_h(V_DI)/(a_h(V_DI) + b_h(V_DI)),...
                b_gi_vec((ii*(APD+DI) ...
                + APD)/dt + jj,ll), k_on, k_off_0, Drug,...
                inact_noninact(ll), GIHH_GRHH_nons);
            
            V_vec((ii*(APD+DI) + APD)/dt + 1 + jj) = V_DI; %storing V
        end %updating through DI
    end

end
%%
xlim_vec = [600,2250];

%BCL 750ms %%files may need to be updated
Simple_nodrug_BCL750 = load('Ten_Tusscher_model_with_drug_output/APDs_all_HH_forms_2_26_18/Simple_nodrug_drug_APs_BCL750_Stimdur06.txt');
Simple_GRI_kon_1e2_BCL750 = load('Ten_Tusscher_model_with_drug_output/APDs_all_HH_forms_2_26_18/Simple_HH_guarded_receptor_inactive_drug_APs_Diffval_1e2_kd10e-6_Drug20e-6_BCL750_stimdur06.txt');
Simple_GRN_kon_1e2_BCL750 = load('Ten_Tusscher_model_with_drug_output/APDs_all_HH_forms_2_26_18/Simple_HH_guarded_receptor_noninactive_drug_APs_Diffval_1e2_kd10e-6_Drug20e-6_BCL750_stimdur06.txt');
Simple_GII_kon_1e2_BCL750 = load('Ten_Tusscher_model_with_drug_output/APDs_all_HH_forms_2_26_18/Simple_HH_gate_immobilization_inactive_drug_APs_Diffval_1e2_kd10e-6_Drug20e-6_BCL750_stimdur06.txt');
Simple_GIN_kon_1e2_BCL750 = load('Ten_Tusscher_model_with_drug_output/APDs_all_HH_forms_2_26_18/Simple_HH_gate_immobilization_noninactive_drug_APs_Diffval_1e2_kd10e-6_Drug20e-6_BCL750_stimdur06.txt');

Vlim = [-100, 40];
gatelim = [0,1];

V_ticks = [-100, -50, 0];
gate_ticks = [0, 0.5, 1];

xrange= [381150, 382750];
xrange2 = [257400, 259000];



%Analytic soln
%%
grey = [0.5 0.5 0.5];
lightred = [1 .4 .4];

figure(6) %comparison of V for full modified ten Tusscher et al model and square wave
p1=plot(t_vec, V_vec, Simple_nodrug_BCL750(:,1) - (381100 - 600), Simple_nodrug_BCL750(:,2));
xlim(xlim_vec); %the - (381100 - 600) is to get the time range for the full
%cell simulations the same as for the square wave approximation
p1(1).Color = [c{5}, transparency_param];
p1(2).Color = c{5};
ylim(Vlim);
box off
set(gca,'XTick',[], 'XTicklabels', []);
set(gca,'YTick',V_ticks, 'YTicklabels', []);

figure(7)
p1 = plot(t_vec, h_gr_vec(:,1), '-', Simple_nodrug_BCL750(:,1) - (381100 - 600),...
    Simple_nodrug_BCL750(:,12),'-');
p1(1).Color = [c{4}, transparency_param];
p1(2).Color = c{4};
xlim(xlim_vec);
ylim(gatelim);
box off
set(gca,'XTick',[], 'XTicklabels', []);
set(gca,'YTick',gate_ticks, 'YTicklabels', []);
%%
%GRI
figure(8)%comparison of h for full modified ten Tusscher et al model and h_infty
%for a square wave
p2 = plot(t_vec, b_gr_vec(:,1), '-', t_vec((APD + DI)/dt +1), b_gr_vec((APD + DI)/dt +1,1), '*',...
    t_vec(2*(APD + DI)/dt +1), b_gr_vec(2*(APD + DI)/dt +1,1), '*', Simple_GRI_kon_1e2_BCL750(:,1) - (381100 - 600), ...
    Simple_GRI_kon_1e2_BCL750(:,14),'-', 'Markersize', 6);
p2(1).Color = [c{1}, transparency_param];
p2(2).Color = [c{1}, transparency_param];
p2(3).Color = [c{1}, transparency_param];
p2(4).Color = c{1};
xlim(xlim_vec);
ylim(gatelim);
box off
set(gca,'XTick',[], 'XTicklabels', []);
set(gca,'YTick',gate_ticks, 'YTicklabels', []);

%GRN
figure(9)
p2 = plot(t_vec, b_gr_vec(:,2), '--', t_vec((APD + DI)/dt +1), b_gr_vec((APD + DI)/dt +1,2), '*',...
    t_vec(2*(APD + DI)/dt +1), b_gr_vec(2*(APD + DI)/dt +1,2), '*', Simple_GRN_kon_1e2_BCL750(:,1) - (381100 - 600),...
    Simple_GRN_kon_1e2_BCL750(:,14), '--', 'Markersize', 6);
p2(1).Color = [c{1}, transparency_param];
p2(2).Color = [c{1}, transparency_param];
p2(3).Color = [c{1}, transparency_param];
p2(4).Color = c{1};
xlim(xlim_vec);
ylim(gatelim);
box off
set(gca,'XTick',[], 'XTicklabels', []);
set(gca,'YTick',gate_ticks, 'YTicklabels', []);

%GII
figure(10)
p2 = plot(t_vec, b_gi_vec(:,1), '-', t_vec((APD + DI)/dt +1), b_gi_vec((APD + DI)/dt +1,1), '*',...
    t_vec(2*(APD + DI)/dt +1), b_gi_vec(2*(APD + DI)/dt +1,1), '*', Simple_GII_kon_1e2_BCL750(:,1) - (381100 - 600),...
    Simple_GII_kon_1e2_BCL750(:,14), '-', 'Markersize', 6);
p2(1).Color = [c{2}, transparency_param];
p2(2).Color = [c{2}, transparency_param];
p2(3).Color = [c{2}, transparency_param];
p2(4).Color = c{2};
xlim(xlim_vec);
ylim(gatelim);
box off
set(gca,'XTick',[], 'XTicklabels', []);
set(gca,'YTick',gate_ticks, 'YTicklabels', []);

%GIN
figure(11)
p2 = plot(t_vec, b_gi_vec(:,2), '--', t_vec((APD + DI)/dt +1), b_gi_vec((APD + DI)/dt +1,2), '*',...
    t_vec(2*(APD + DI)/dt +1), b_gi_vec(2*(APD + DI)/dt +1,2), '*', Simple_GIN_kon_1e2_BCL750(:,1) - (381100 - 600),...
    Simple_GIN_kon_1e2_BCL750(:,14), '--', 'Markersize', 6);
p2(1).Color = [c{2}, transparency_param];
p2(2).Color = [c{2}, transparency_param];
p2(3).Color = [c{2}, transparency_param];
p2(4).Color = c{2};
xlim(xlim_vec);
ylim(gatelim);
box off
set(gca,'XTick',[], 'XTicklabels', []);
set(gca,'YTick',gate_ticks, 'YTicklabels', []);

%% V dynamics for graphical abstract
figure(12) %comparison of V for full modified ten Tusscher et al model and square wave
p1=plot(Simple_nodrug_BCL750(:,1) - (381100 - 600), Simple_nodrug_BCL750(:,2), 'k');
xlim(xlim_vec); %the - (381100 - 600) is to get the time range for the full
%cell simulations the same as for the square wave approximation
%p1(1).Color = c{4};
ylim(Vlim);
xlabel('time')
ylabel('V (mV)')
box off
set(gca,'XTick',[], 'XTicklabels', []);
set(gca,'YTick',V_ticks, 'FontName', 'Arial', ...
        'FontSize', 5, ...
        'LineWidth', 0.5);
