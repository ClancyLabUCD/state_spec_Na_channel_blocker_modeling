%Steffen Docken (Lewis Lab). All rights reserved.
%  Published in the Journal of Theoretical Biology under the title "Rate-dependent effects of state-specific sodium channel blockers in cardiac tissue: Insights from idealized models"
%1-30-17
%Code to calculate and graph tau_I_O_approx and I_O_infty_approx for
%different drug binding models and binding rate constants.  

clear
close all

transparency_param = .5;

V_stim = -20; %mV
V_hold = -80; %mV

T = 310; %K

xlim_vec = [5e-4, 5e1];
ylim_vec_I_O = [-0.1,1];
ylim_vec_tau = [2e-3,4e-1];

x_tick_vals = [0.01, 1];
y_tick_vals_tau = [0.01, 0.1];
y_tick_vals_I_O = [0, 0.50, 1];

Drug = 20e-6; %concentration of drug in M

k_on_init = 1e5; %M^-1ms^-1
k_D0 = 10e-6; %ms^-1

f_vec = logspace(-2, 2, 21); %vector of factors k_on will be multiplied by. from
%10^-2 to 10^2

N_f = length(f_vec);

tau_b_vec = zeros(N_f,1); %will hold tau_b for each f value

tau_I_O_approx_GIHHI_vec = zeros(N_f,1); %vec for Gate Immobilization 
%inactivated state binding, HH formulation

tau_I_O_approx_GIFI_vec = zeros(N_f,1); %vec for Gate Immobilization 
%inactivated state binding, Full formulation

tau_I_O_approx_GIHHN_vec = zeros(N_f,1);%vec for Gate Immobilization 
%noninactivated state binding, HH formulation

tau_I_O_approx_GIFN_vec = zeros(N_f,1);%vec for Gate Immobilization 
%noninactivated state binding, Full formulation

tau_I_O_approx_GRHHI_vec = zeros(N_f,1); %vec for Guarded Receptor 
%inactivated state binding, HH formulation

tau_I_O_approx_GRFI_vec = zeros(N_f,1);%vec for Guarded Receptor
%inactivated state binding, Full formulation

tau_I_O_approx_GRHHN_vec = zeros(N_f,1);%vec for Guarded Receptor 
%noninactivated state binding, HH formulation

tau_I_O_approx_GRFN_vec = zeros(N_f,1);%vec for Guarded Receptor 
%noninactivated state binding, Full formulation


a_h = @(V)a_h_6_14_2016(V, T);
b_h = @(V)b_h_6_14_2016(V, T);

h_infty_nons_hold = a_h(V_hold)/(a_h(V_hold) + b_h(V_hold)); %steady state 
%value of h at the holding voltage if drug binding is not state specific

tau_h = 1/(a_h(V_stim) + b_h(V_stim));%time constant of the I_O
%state at the stimulus potential when no drug is present
I_O_infty_nodrug = a_h(V_stim)/(a_h(V_stim) + b_h(V_stim));%steady state
%fraction of channels in the I_O state at the stimulus potential when no 
%drug is present

for ii = 1:N_f
    k_on = k_on_init*f_vec(ii);
    
    k_off_0 = k_on*k_D0; %setting binding rates for this loop iteration
    
    [~, ~, ~, tau_b_vec(ii)] = HH_infty_tau(V_stim, I_O_infty_nodrug,...
        k_on, k_off_0, Drug, 1, 2);
    
    %% Gate Immobilization Inactivated state
    
    %HH Model
    [h_infty_hold, ~, ~, ~] = HH_infty_tau(V_hold, ...
        h_infty_nons_hold, k_on, k_off_0, Drug, 0, 0); %calculating
    %h_infty and b_infty for end of period where cell is held at holding
    %potential
    
    tau_I_O_approx_GIHHI_vec(ii) =...
        I_O_approx_dynamics(V_stim, h_infty_hold,...
        Drug, k_on, k_off_0, 0, 0, 0);
    
    %Full formulation
    [~, ~, ~, X_ss] = Model_Matrices(V_hold, k_on, k_off_0, ...
        Drug, 0, 0);
    h_infty_hold = X_ss(1);%calculating h_infty for end of 
    %period where cell is held at holding potential
    
    tau_I_O_approx_GIFI_vec(ii) =...
        I_O_approx_dynamics(V_stim, h_infty_hold,...
        Drug, k_on, k_off_0, 0, 0, 1);
    
    %% Gate Immobilization Non-inactivated state
    
    %HH Model
    [h_infty_hold, ~, ~, ~] = HH_infty_tau(V_hold, ...
        h_infty_nons_hold, k_on, k_off_0, Drug, 1, 0); %calculating
    %h_infty and b_infty for end of period where cell is held at holding
    %potential
    
    tau_I_O_approx_GIHHN_vec(ii) =...
        I_O_approx_dynamics(V_stim, h_infty_hold, ...
        Drug, k_on, k_off_0, 1, 0, 0);
    
    %Full formulation
    [~, ~, ~, X_ss] = Model_Matrices(V_hold, k_on, k_off_0, ...
        Drug, 1, 0);
    h_infty_hold = X_ss(1);%calculating h_infty for end of 
    %period where cell is held at holding potential
    
    tau_I_O_approx_GIFN_vec(ii) =...
        I_O_approx_dynamics(V_stim, h_infty_hold, ...
        Drug, k_on, k_off_0, 1, 0, 1);
    
    %% Guarded Receptor Inactivated state
    
    %HH Model
    [h_infty_hold, ~, ~, ~] = HH_infty_tau(V_hold, ...
        h_infty_nons_hold, k_on, k_off_0, Drug, 0, 1); %calculating
    %h_infty and b_infty for end of period where cell is held at holding
    %potential
    
    tau_I_O_approx_GRHHI_vec(ii) =...
        I_O_approx_dynamics(V_stim, h_infty_hold, ...
        Drug, k_on, k_off_0, 0, 1, 0);
    
    %Full formulation
    [~, ~, ~, X_ss] = Model_Matrices(V_hold, k_on, k_off_0, ...
        Drug, 0, 1);
    h_infty_hold = X_ss(1);%calculating h_infty for end of 
    %period where cell is held at holding potential
    
    tau_I_O_approx_GRFI_vec(ii) =...
        I_O_approx_dynamics(V_stim, h_infty_hold, ...
        Drug, k_on, k_off_0, 0, 1, 1);
    
    %% Guarded Receptor Non-inactivated state
    
    %HH Model
    [h_infty_hold, ~, ~, ~] = HH_infty_tau(V_hold, ...
        h_infty_nons_hold, k_on, k_off_0, Drug, 1, 1); %calculating
    %h_infty and b_infty for end of period where cell is held at holding
    %potential
    
    tau_I_O_approx_GRHHN_vec(ii) =...
        I_O_approx_dynamics(V_stim, h_infty_hold, ...
        Drug, k_on, k_off_0, 1, 1, 0);
    
    %Full formulation
    [~, ~, ~, X_ss] = Model_Matrices(V_hold, k_on, k_off_0, ...
        Drug, 1, 1);
    h_infty_hold = X_ss(1);%calculating h_infty for end of 
    %period where cell is held at holding potential
    
    tau_I_O_approx_GRFN_vec(ii) =...
        I_O_approx_dynamics(V_stim, h_infty_hold, ...
        Drug, k_on, k_off_0, 1, 1, 1);
    
end

h = plot(1:10,1:10,1:10,2:11,1:10,3:12,1:10,4:13,1:10,5:14,1:10,6:15);
c = get(h,'Color'); %figure so I can reference colors in the following figures.



figure(11)
p1 = loglog(tau_b_vec, tau_I_O_approx_GRHHI_vec, '-', tau_b_vec, ...
    tau_I_O_approx_GRFI_vec, '-', tau_b_vec,...
    tau_I_O_approx_GRHHN_vec, '--', tau_b_vec, tau_I_O_approx_GRFN_vec, '--',...
    tau_h*[1,1],[2e-3, 2.5e-3], 'r', 'markersize', 3);
xlim(xlim_vec);
ylim(ylim_vec_tau);
set(gca,'XTick',x_tick_vals);
set(gca,'YTick',y_tick_vals_tau);
%title('GR tau')
%legend('GRHHI', 'GRFI', 'GRHHN', 'GRFN')

p1(1).Color = c{1};
p1(2).Color = [c{1}, transparency_param];
p1(3).Color = c{1};
p1(4).Color = [c{1}, transparency_param];
box off

figure(13)
p1 = loglog(tau_b_vec, tau_I_O_approx_GIHHI_vec, '-', tau_b_vec, ...
    tau_I_O_approx_GIFI_vec, '-', tau_b_vec,...
    tau_I_O_approx_GIHHN_vec, '--', tau_b_vec, tau_I_O_approx_GIFN_vec, '--',...
    tau_h*[1,1],[2e-3, 2.5e-3], 'r', 'markersize', 3);
xlim(xlim_vec);
ylim(ylim_vec_tau);
set(gca,'XTick',x_tick_vals);
set(gca,'YTick',y_tick_vals_tau);
%title('GI tau');
%legend('GIHHI', 'GIFI', 'GIHHN', 'GIFN')

p1(1).Color = c{2};
p1(2).Color = [c{2}, transparency_param];
p1(3).Color = c{2};
p1(4).Color = [c{2}, transparency_param];
box off
