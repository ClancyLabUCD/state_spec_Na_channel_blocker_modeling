%Steffen Docken (Lewis Lab). All rights reserved.
%  Published in the Journal of Theoretical Biology under the title "Rate-dependent effects of state-specific sodium channel blockers in cardiac tissue: Insights from idealized models"
%1-31-17
%Code to visualize Upstroke Velocity as a function of tau_b
clear
close all

xlim_vec = [1e-3, 1e4];
x_tick_vals = [0.01, 1, 100];

ylim_vec = [80,300];
y_tick_vals = [100, 200, 300]; %figure axis specifications

Drug = 20e-6; %concentration of drug in M

k_D0 = 10e-6; %ms^-1

transparency_param = .5;

%load Upstroke_kon1e5_kd10eneg6_Drug20eneg6_dt01_Stimdur06
load Upstroke_vs_taub_BCL1000_data_12_7_17.mat

V_stim = -20; %Transmembrane potential at which tau_b is calculated.


tau_b_vec = zeros(length(GRHHI(:,1)),1); %vector to hold tau_b(-20) values

for ii = 1:length(GRHHI(:,1))
    k_on = GRHHI(ii,1);
    
    k_off_0 = k_on*k_D0; %setting binding rates for this loop iteration
    
    [~, ~, tau_h_val, tau_b_vec(ii)] = HH_infty_tau(V_stim, 1,...
        k_on, k_off_0, Drug, 1, 2); %This calculation for tau_b 
    %assumes there has been a rapid increase from a hyperpolarized
    %potential to V = -20mV, and therefore h is approximately 1.
end
h = plot(1:10,1:10,1:10,2:11,1:10,3:12,1:10,4:13,1:10,5:14,1:10,6:15);
c = get(h,'Color');

figure(1)
p1 = semilogx(tau_b_vec, nodrug(1,4)*ones(length(tau_b_vec),1),'k',...
    tau_b_vec, GRFI(:,4), '-', tau_b_vec, GRFN(:,4), '--',...
    tau_b_vec, GRHHI(:,4), '-', tau_b_vec, GRHHN(:,4), '--', ...
    tau_h_val*[1,1], [80,83], 'r', 'markersize', 3);
title('GR'); %Peak Upstroke Velocity vs. tau_b for Guarded Receptor models
p1(2).Color = [c{1}, transparency_param];
p1(3).Color = [c{1}, transparency_param];
p1(4).Color = c{1};
p1(5).Color = c{1};
legend('nodrug', 'GRFI', 'GRFN', 'GRHHI', 'GRHHN')
xlim(xlim_vec);
ylim(ylim_vec); 
set(gca,'XTick',x_tick_vals);%, 'XTickLabel', []);
set(gca,'YTick',y_tick_vals);
box off

figure(2)
p1 = semilogx(tau_b_vec, nodrug(1,4)*ones(length(tau_b_vec),1),'k',...
    tau_b_vec, GIFI(:,4), '-', tau_b_vec, GIFN(:,4), '--',...
    tau_b_vec, GIHHI(:,4), '-', tau_b_vec, GIHHN(:,4), '--', ...
    tau_h_val*[1,1], [80,83], 'r', 'markersize', 3);
title('GI'); %Peak Upstroke Velocity vs. tau_b for Gate Immobilization models
p1(2).Color = [c{2}, transparency_param];
p1(3).Color = [c{2}, transparency_param];
p1(4).Color = c{2};
p1(5).Color = c{2};
xlim(xlim_vec);
ylim(ylim_vec); 
set(gca,'XTick',x_tick_vals);%, 'XTickLabel', []);
set(gca,'YTick',y_tick_vals);
box off


