%Steffen Docken (Lewis Lab). All rights reserved.
%  Published in the Journal of Theoretical Biology under the title "Rate-dependent effects of state-specific sodium channel blockers in cardiac tissue: Insights from idealized models"
%3-23-17
%This code calculates the b* values given the square wave approximation of
%a train of APs and using the APD and DI values from the restitution curves
%for each of the drug binding models

clear
close all

load Upstroke_v_Period_3_tau_b.mat

V_DI = -80;
V_A = 20; %in mV

R = 8314.472;	% J/mol*K
T = 310;% K
F = 96485.3415; 

Drug = 20e-6; %Drug concentration in M

Diffusion_vec = [1e1, 1e3]; %in M^-1ms^-1 

k_D0 = 10e-6; %in M

inact_noninact_vec = [0, 1]; %0 = inactive state binding, 1 = noninactive
%state binding

transparency_param = .5;

N = length(nodrug(:,1)); %number of BCLs used

GR_restitution = zeros(N,4,2,2);
GI_restitution = zeros(N,4,2,2); %These arrays will hold the data from 
%Upstroke_v_Period_3_tau_b.mat

GR_restitution(:,:,1,1) = GRHHI_1e1; %inact binding and Diff =1e1
GR_restitution(:,:,1,2) = GRHHI_1e3; %inact binding and Diff =1e3
GR_restitution(:,:,2,1) = GRHHN_1e1; %non-inact binding and Diff =1e1
GR_restitution(:,:,2,2) = GRHHN_1e3; %non-inact binding and Diff =1e3

GI_restitution(:,:,1,1) = GIHHI_1e1; %inact binding and Diff =1e1
GI_restitution(:,:,1,2) = GIHHI_1e3; %inact binding and Diff =1e3
GI_restitution(:,:,2,1) = GIHHN_1e1; %non-inact binding and Diff =1e1
GI_restitution(:,:,2,2) = GIHHN_1e3; %non-inact binding and Diff =1e3


b_star_gr_mat = zeros(N, 2, 2);
b_star_gi_mat = zeros(N, 2, 2); %These arrays will hold b_star for the 
%various drug binding models with different BCL (each row)
%1st column will correspond to inact, 2nd 
%to non-inact binding, 1st 3rd-dim will correspond to Diff = 1e1, 2nd to
%Diff = 1e3

for ii = 1:2 %looping through Diffusion values
    k_on = Diffusion_vec(ii);
    k_off_0 = Diffusion_vec(ii)*k_D0; %defining params from options
    
    for ll = 1:2 %looping through state binding
        inact_noninact = inact_noninact_vec(ll);
        
        %% guarded receptor binding
        gi_gr_nons = 1; 
        %HH formulation

        for kk = 1:N
            [~,b_star_gr_mat(kk,ll,ii)] = HH_b_star(GR_restitution(kk,4,ll,ii),...
                GR_restitution(kk,1,ll, ii)-GR_restitution(kk,4,ll,ii), V_A,...
                V_DI, k_on, k_off_0, Drug, inact_noninact, ...
                gi_gr_nons); %calculating b_star for each BCL
        end
        
        %% gate immobilization binding
        gi_gr_nons = 0; 
        %HH formulation

        for kk = 1:N
            [~,b_star_gi_mat(kk,ll,ii)] = HH_b_star(GI_restitution(kk,4,ll,ii),...
                GI_restitution(kk,1,ll, ii)-GI_restitution(kk,4,ll,ii), V_A,...
                V_DI, k_on, k_off_0, Drug, inact_noninact, ...
                gi_gr_nons); %calculating b_star for each BCL
        end
        
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


figure(5);
p1 = plot(nodrug(7:end,1), GR_restitution(7:end,3,1,1)./nodrug(7:end,3), '-',...
    nodrug(7:end,1), GR_restitution(7:end,3,2,1)./nodrug(7:end,3), '--',...
    nodrug(7:end,1), GI_restitution(7:end,3,1,1)./nodrug(7:end,3), '-',...
    nodrug(7:end,1), GI_restitution(7:end,3,2,1)./nodrug(7:end,3), '--',...
    nodrug(7:end,1), 1-b_star_gr_mat(7:end,1,1), '-', nodrug(7:end,1), 1-b_star_gr_mat(7:end,2,1), '--',...
    nodrug(7:end,1), 1-b_star_gi_mat(7:end,1,1), '-', nodrug(7:end,1), 1-b_star_gi_mat(7:end,2,1), '--');
%title('k_{on} = 1e1')
ylim(ylim_vec);
p1(1).Color = c{1};
p1(2).Color = c{1};
p1(3).Color = c{2};
p1(4).Color = c{2};
p1(5).Color = [c{1}, transparency_param];
p1(6).Color = [c{1}, transparency_param];
p1(7).Color = [c{2}, transparency_param];
p1(8).Color = [c{2}, transparency_param];
box off;
set(gca,'XTick',x_tick_vals);
set(gca,'YTick',y_tick_vals, 'YTickLabels',[]);
yyaxis right
ax = gca;
ax.YColor = transparency_param*[1, 1, 1];
ax.YTick = y_tick_vals;
ax.YTickLabel = [];
%b_star for the various drug binding models (solid lines) compared to peak
%upstroke velocity normalized by peak upstroke velocity in the no drug case
%(-0) for k_on = 10 M^-1ms^-1

figure(6)
p1 = plot(nodrug(7:end,1), GR_restitution(7:end,3,1,2)./nodrug(7:end,3), '-',...
    nodrug(7:end,1), GR_restitution(7:end,3,2,2)./nodrug(7:end,3), '--',...
    nodrug(7:end,1), GI_restitution(7:end,3,1,2)./nodrug(7:end,3), '-',...
    nodrug(7:end,1), GI_restitution(7:end,3,2,2)./nodrug(7:end,3), '--', ...
    nodrug(7:end,1), 1-b_star_gr_mat(7:end,1,2), '-', nodrug(7:end,1), 1-b_star_gr_mat(7:end,2,2), '--',...
    nodrug(7:end,1), 1-b_star_gi_mat(7:end,1,2), '-', nodrug(7:end,1), 1-b_star_gi_mat(7:end,2,2), '--');
%title('k_{on} = 1e3')
ylim(ylim_vec);
p1(1).Color = c{1};
p1(2).Color = c{1};
p1(3).Color = c{2};
p1(4).Color = c{2};
p1(5).Color = [c{1}, transparency_param];
p1(6).Color = [c{1}, transparency_param];
p1(7).Color = [c{2}, transparency_param];
p1(8).Color = [c{2}, transparency_param];
box off;
set(gca,'XTick',x_tick_vals);
set(gca,'YTick',y_tick_vals, 'YTickLabels',[]);
yyaxis right
ax = gca;
ax.YColor = transparency_param*[1, 1, 1];
ax.YTick = y_tick_vals;
ax.YTickLabel = [];
%b_star for the various drug binding models (solid lines) compared to peak
%upstroke velocity normalized by peak upstroke velocity in the no drug case
%(-0) for k_on = 10^3 M^-1ms^-1

