%Steffen Docken (Lewis Lab). All rights reserved.
%  Published in the Journal of Theoretical Biology under the title "Rate-dependent effects of state-specific sodium channel blockers in cardiac tissue: Insights from idealized models"
%7-4-16
%Code to compute steady-state b as a funciton of V for each model

clear

V_vec = -120:2:40;

N = length(V_vec);

R = 8314.472;	% mJ/mol*K
T = 310;% K
F = 96485.3415; %Col/mol

h = plot(1:10,1:10,1:10,2:11,1:10,3:12,1:10,4:13,1:10,5:14,1:10,6:15);
c = get(h,'Color'); %This graph is just so I can reference the colors in 
%future plots

a_h = @(V)a_h_6_14_2016(V, T);
b_h = @(V)b_h_6_14_2016(V, T);


Drug = 20e-6; %Drug concentration in M

Diffusion = 1e2; %in M^-1ms^-1

k_D0 = 10e-6;%in M 

k_on = Diffusion;
k_D0 = Diffusion*k_D0; %defining params from options

k_off = @(V)(k_D0*exp(-0.7*V*F/(R*T))); %V-dependent %(k_D0); %non-V-dependent 

%For Guarded Receptor model
h_infty = @(V)(a_h(V)./(a_h(V)+b_h(V))); %steady state value of h

tau_h = @(V)(1./(a_h(V)+b_h(V))); %time constant of h

b_infty_equal = @(V,D)(D*k_on./(D*k_on + k_off(V))); %steady state value of b

tau_b_equal = @(V,D)(1./(D*k_on + k_off(V))); %time constant of b

%for Gate immobilization model with binding to the inactive state
b_infty_gi_inact = @(V,D)(b_infty_equal(V,D).*(1-h_infty(V))./(1-b_infty_equal(V,D).*...
    h_infty(V))); %steady state value of b

tau_b_gi_inact = @(V,D)(tau_b_equal(V,D)./(1-h_infty(V).*b_infty_equal(V,D)));
%time constant of drug binding (assuming h = h_infty)

h_infty_bar = @(V)(b_h(V)./(a_h(V) + b_h(V))); %steady state fraction of 
%channels inactivated in the normal H-H model

%for Gate immobilization model with binding to non-inactivated channels
b_infty_gi_noninact = @(V,D)(b_infty_equal(V,D).*(1-h_infty_bar(V))./(1-b_infty_equal(V,D).*...
    h_infty_bar(V))); %steady state fraction of channels bound to drug

tau_b_gi_noninact = @(V,D)(tau_b_equal(V,D)./(1-h_infty_bar(V).*b_infty_equal(V,D)));
%time constant of drug binding (assuming h = h_infty)

ylim_tau_b_vec = [10, 1e7];
y_tick_vals = [0, 0.25, 0.5, 0.75, 1];
y_tau_b_tick_vals = [10, 1e3, 1e5, 1e7];

V_DI_plot = [-80, -80];
V_AP_plot = [20, 20];

grey = [0.5 0.5 0.5];

xrange = [-100,40];
yrange = [0,1];
x_tick_vals = [-100, -80, -60, -40, -20, 0, 20, 40];

figure(6)
p1 = plot(V_vec, b_infty_equal(V_vec, Drug), 'k',...
    V_vec, b_infty_equal(V_vec, Drug) + 0.01, '-',...
    V_vec, b_infty_equal(V_vec, Drug)+.02, '--', ...
    V_vec, b_infty_gi_inact(V_vec, Drug), '-',...
    V_vec, b_infty_gi_noninact(V_vec, Drug), '--',...
    V_DI_plot, [0,1], '--', V_AP_plot, [0,1], '--');
%xlabel('V (mV)');
%ylabel('b_\infty');
%legend('GR models', 'GR models', 'GII', 'GIN');
p1(2).Color = c{1};
p1(3).Color = c{1};
p1(4).Color = c{2};
p1(5).Color = c{2};
p1(6).Color = grey;
p1(7).Color = grey;
xlim(xrange);
ylim(yrange);
box off
set(gca,'XTick', x_tick_vals, 'XTickLabel', []);
set(gca,'YTick',[0, 0.25, 0.5, 0.75, 1]);
%comparison of b_infty for all models 

tau_gi_inact = zeros(N,1);
tau_gi_noninact = zeros(N,1); %arrays to hold the time constant of 
%drug binding in the GII and GIN models

tau_grS_inact = zeros(N,1);
tau_grS_noninact = zeros(N,1);%arrays to hold the time constant of 
%drug binding in the GRI and GRN models

tau_equal = zeros(N,1); %time constant of drug binding for a non-state 
%specific blocker

gi_gr = 0; %gate immobilization model

for ii = 1:length(V_vec)
    
    tau_gi_inact(ii) = tau_b_gi_inact(V_vec(ii),Drug);
    tau_gi_noninact(ii) = tau_b_gi_noninact(V_vec(ii),Drug);
    %arrays to hold the time constant of drug binding in the GII and GIN 
    %models
    
    tau_grS_inact(ii) = tau_b_equal(V_vec(ii),Drug)/(1-h_infty(V_vec(ii)));
    tau_grS_noninact(ii) = tau_b_equal(V_vec(ii),Drug)/(h_infty(V_vec(ii)));
    %arrays to hold the time constant of drug binding in the GRI and GRN 
    %models
    
    tau_equal(ii) = tau_b_equal(V_vec(ii), Drug);%time constant of drug 
    %binding for a non-state specific blocker
end


figure(12)
p2 = semilogy(V_vec, tau_b_equal(V_vec, Drug), 'k', V_vec, tau_grS_inact, '-', V_vec, tau_grS_noninact, '--',...
    V_vec, tau_gi_inact, '-', V_vec, tau_gi_noninact, '--',...
    V_DI_plot, ylim_tau_b_vec, '--', V_AP_plot, ylim_tau_b_vec, '--')%,...
%xlabel('V(mV)');
%ylabel('\tau_b (ms)')
xlim(xrange);
ylim(ylim_tau_b_vec);
%legend('no Drug', 'GRI', 'GRN', 'GII', 'GIN');
box off
set(gca,'XTick',x_tick_vals, 'XTicklabels', []);%x_tick_labels);
set(gca,'YTick',y_tau_b_tick_vals);

p2(2).Color = c{1};
p2(3).Color = c{1};
p2(4).Color = c{2};
p2(5).Color = c{2};
p2(6).Color = grey;
p2(7).Color = grey;
%comparing the time constant of drug binding in the various drug binding
%models

