%Steffen Docken (Lewis Lab). All rights reserved.
%  Published in the Journal of Theoretical Biology under the title "Rate-dependent effects of state-specific sodium channel blockers in cardiac tissue: Insights from idealized models"
%1-25-17
%Code to illustrate Upstroke vs. Period data for different models

load Upstroke_v_Period_3_tau_b.mat
load Upstroke_rest_30s_3_tau_b.mat
%load Upstroke_v_Period_kon1e5_kd10eneg6_Drug20eneg6_dt01_Stimdur06.mat

fake_rest = 1.05e3; %a value used for the BCL of 30s to make all data fit 
%in one graph

transparency_param = .5;

nodrug_BL(end-1:end,1) = fake_rest;

GRHHI_1e1_BL(:,1) = fake_rest;
GRHHI_1e3_BL(:,1) = fake_rest;
GRHHI_1e5_BL(:,1) = fake_rest;

GRHHN_1e1_BL(:,1) = fake_rest;
GRHHN_1e3_BL(:,1) = fake_rest;
GRHHN_1e5_BL(:,1) = fake_rest;

GIHHI_1e1_BL(:,1) = fake_rest;
GIHHI_1e3_BL(:,1) = fake_rest;
GIHHI_1e5_BL(:,1) = fake_rest;

GIHHN_1e1_BL(:,1) = fake_rest;
GIHHN_1e3_BL(:,1) = fake_rest;
GIHHN_1e5_BL(:,1) = fake_rest;

GRFI_1e1_BL(:,1) = fake_rest;
GRFI_1e3_BL(:,1) = fake_rest;
GRFI_1e5_BL(:,1) = fake_rest;

GRFN_1e1_BL(:,1) = fake_rest;
GRFN_1e3_BL(:,1) = fake_rest;
GRFN_1e5_BL(:,1) = fake_rest;

GIFI_1e1_BL(:,1) = fake_rest;
GIFI_1e3_BL(:,1) = fake_rest;
GIFI_1e5_BL(:,1) = fake_rest;

GIFN_1e1_BL(:,1) = fake_rest;
GIFN_1e3_BL(:,1) = fake_rest;
GIFN_1e5_BL(:,1) = fake_rest;

h = plot(1:10,1:10,1:10,2:11,1:10,3:12,1:10,4:13,1:10,5:14,1:10,6:15);
c = get(h,'Color'); %This graph is just so I can reference the colors in 
%future plots

xlim_vec = [250, 1050];
ylim_vec = [100, 350];
x_tick_vals = [400, 600, 800, 1000];
y_tick_vals = [100, 200, 300];

figure(1)
p1 = plot(GRFI_1e1(7:end,1), GRFI_1e1(7:end,3), '-', GRFN_1e1(7:end,1), GRFN_1e1(7:end,3), '--', ...
    GIFI_1e1(7:end,1), GIFI_1e1(7:end,3), '-', GIFN_1e1(7:end,1), GIFN_1e1(7:end,3), '--',...
    nodrug(7:end,1), nodrug(7:end,3), '-k',...
    GRHHI_1e1(7:end,1), GRHHI_1e1(7:end,3)+2, '-', GRHHN_1e1(7:end,1), GRHHN_1e1(7:end,3)+2, '--', ...
    GIHHI_1e1(7:end,1), GIHHI_1e1(7:end,3)+2, '-', GIHHN_1e1(7:end,1), GIHHN_1e1(7:end,3)+2, '--',...
    GRFI_1e1_BL(1,1), GRFI_1e1_BL(1,3), 'O', GRFN_1e1_BL(1,1), GRFN_1e1_BL(1,3), '.', ...
    GIFI_1e1_BL(1,1), GIFI_1e1_BL(1,3), 'O', GIFN_1e1_BL(1,1), GIFN_1e1_BL(1,3), '.',...
    GRHHI_1e1_BL(1,1), GRHHI_1e1_BL(1,3), 'O', GRHHN_1e1_BL(1,1), GRHHN_1e1_BL(1,3), '.', ...
    GIHHI_1e1_BL(1,1), GIHHI_1e1_BL(1,3), 'O', GIHHN_1e1_BL(1,1), GIHHN_1e1_BL(1,3), '.',...
    nodrug_BL(1,1), nodrug_BL(1,3), 'kO', 'MarkerSize', 3);
xlim(xlim_vec);
ylim(ylim_vec);
set(gca,'XTick',x_tick_vals);
set(gca,'YTick',y_tick_vals, 'YTicklabels', []);
box off

p1(1).Color = [c{1}, transparency_param];
p1(2).Color = [c{1}, transparency_param];
p1(3).Color = [c{2}, transparency_param];
p1(4).Color = [c{2}, transparency_param];
p1(6).Color = c{1};
p1(7).Color = c{1};
p1(8).Color = c{2};
p1(9).Color = c{2};

p1(10).Color = [c{1}, transparency_param];
p1(11).Color = [c{1}, transparency_param];
p1(12).Color = [c{2}, transparency_param];
p1(13).Color = [c{2}, transparency_param];
p1(14).Color = c{1};
p1(15).Color = c{1};
p1(16).Color = c{2};
p1(17).Color = c{2};


%title('k_on = 1e1')
%legend('GRI', 'GRN', 'GII', 'GIN', 'no drug');

%%

figure(2)
p1 = plot(GRFI_1e3(7:end,1), GRFI_1e3(7:end,3), '-', GRFN_1e3(7:end,1), GRFN_1e3(7:end,3), '--', ...
    GIFI_1e3(7:end,1), GIFI_1e3(7:end,3), '-', GIFN_1e3(7:end,1), GIFN_1e3(7:end,3), '--',...
    nodrug(7:end,1), nodrug(7:end,3), 'k-',...
    GRHHI_1e3(7:end,1), GRHHI_1e3(7:end,3), '-', GRHHN_1e3(7:end,1), GRHHN_1e3(7:end,3), '--', ...
    GIHHI_1e3(7:end,1), GIHHI_1e3(7:end,3), '-', GIHHN_1e3(7:end,1), GIHHN_1e3(7:end,3), '--',...
    GRFI_1e3_BL(1,1), GRFI_1e3_BL(1,3), 'O', GRFN_1e3_BL(1,1), GRFN_1e3_BL(1,3), '.', ...
    GIFI_1e3_BL(1,1), GIFI_1e3_BL(1,3), 'O', GIFN_1e3_BL(1,1), GIFN_1e3_BL(1,3), '.',...
    GRHHI_1e3_BL(1,1), GRHHI_1e3_BL(1,3), 'O', GRHHN_1e3_BL(1,1), GRHHN_1e3_BL(1,3), '.', ...
    GIHHI_1e3_BL(1,1), GIHHI_1e3_BL(1,3), 'O', GIHHN_1e3_BL(1,1), GIHHN_1e3_BL(1,3), '.',...
    nodrug_BL(1,1), nodrug_BL(1,3), 'kO', 'MarkerSize', 3);
xlim(xlim_vec);
ylim(ylim_vec);
set(gca,'XTick',x_tick_vals);
set(gca,'YTick',y_tick_vals, 'YTicklabels', []);
box off

p1(1).Color = [c{1}, transparency_param];
p1(2).Color = [c{1}, transparency_param];
p1(3).Color = [c{2}, transparency_param];
p1(4).Color = [c{2}, transparency_param];
p1(6).Color = c{1};
p1(7).Color = c{1};
p1(8).Color = c{2};
p1(9).Color = c{2};

p1(10).Color = [c{1}, transparency_param];
p1(11).Color = [c{1}, transparency_param];
p1(12).Color = [c{2}, transparency_param];
p1(13).Color = [c{2}, transparency_param];
p1(14).Color = c{1};
p1(15).Color = c{1};
p1(16).Color = c{2};
p1(17).Color = c{2};

%title('k_on = 1e3')
%legend('GRI', 'GRN', 'GII', 'GIN', 'no drug');

%%

figure(3)
p1 = plot(GRFI_1e5(7:end,1), GRFI_1e5(7:end,3), '-', GRFN_1e5(7:end,1), GRFN_1e5(7:end,3), '--', ...
    GIFI_1e5(7:end,1), GIFI_1e5(7:end,3), '-', GIFN_1e5(7:end,1), GIFN_1e5(7:end,3), '--',...
    nodrug(7:end,1), nodrug(7:end,3), 'k-',...
    GRHHI_1e5(7:end,1), GRHHI_1e5(7:end,3), '-', GRHHN_1e5(7:end,1), GRHHN_1e5(7:end,3), '--', ...
    GIHHI_1e5(7:end,1), GIHHI_1e5(7:end,3), '-', GIHHN_1e5(7:end,1), GIHHN_1e5(7:end,3), '--',...
    GRFI_1e5_BL(1,1), GRFI_1e5_BL(1,3), 'O', GRFN_1e5_BL(1,1), GRFN_1e5_BL(1,3), '.', ...
    GIFI_1e5_BL(1,1), GIFI_1e5_BL(1,3), 'O', GIFN_1e5_BL(1,1), GIFN_1e5_BL(1,3), '.',...
    GRHHI_1e5_BL(1,1), GRHHI_1e5_BL(1,3), 'O', GRHHN_1e5_BL(1,1), GRHHN_1e5_BL(1,3), '.', ...
    GIHHI_1e5_BL(1,1), GIHHI_1e5_BL(1,3), 'O', GIHHN_1e5_BL(1,1), GIHHN_1e5_BL(1,3), '.',...
    nodrug_BL(1,1), nodrug_BL(1,3), 'kO', 'MarkerSize', 3);
xlim(xlim_vec);
ylim(ylim_vec);
set(gca,'XTick',x_tick_vals);
set(gca,'YTick',y_tick_vals, 'YTicklabels', []);
box off

p1(1).Color = [c{1}, transparency_param];
p1(2).Color = [c{1}, transparency_param];
p1(3).Color = [c{2}, transparency_param];
p1(4).Color = [c{2}, transparency_param];
p1(6).Color = c{1};
p1(7).Color = c{1};
p1(8).Color = c{2};
p1(9).Color = c{2};

p1(10).Color = [c{1}, transparency_param];
p1(11).Color = [c{1}, transparency_param];
p1(12).Color = [c{2}, transparency_param];
p1(13).Color = [c{2}, transparency_param];
p1(14).Color = c{1};
p1(15).Color = c{1};
p1(16).Color = c{2};
p1(17).Color = c{2};

%title('k_on = 1e5')
%legend('GRI', 'GRN', 'GII', 'GIN', 'no drug');

%%


figure(4) % a simpler version of figure(1)
p1 = plot(GRFI_1e1(7:end,1), GRFI_1e1(7:end,3), '-', GRFN_1e1(7:end,1), GRFN_1e1(7:end,3), '--', ...
    GIFI_1e1(7:end,1), GIFI_1e1(7:end,3), '-', GIFN_1e1(7:end,1), GIFN_1e1(7:end,3), '--',...
    nodrug(7:end,1), nodrug(7:end,3), '-k',...
    GRFI_1e1_BL(1,1), GRFI_1e1_BL(1,3), 'O', GRFN_1e1_BL(1,1), GRFN_1e1_BL(1,3), '.', ...
    GIFI_1e1_BL(1,1), GIFI_1e1_BL(1,3), 'O', GIFN_1e1_BL(1,1), GIFN_1e1_BL(1,3), '.',...
    nodrug_BL(1,1), nodrug_BL(1,3), 'kO', 'MarkerSize', 3);
xlim(xlim_vec);
ylim(ylim_vec);
set(gca,'XTick',x_tick_vals);
set(gca,'YTick',y_tick_vals);
xlabel('Basic Cycle Length (ms)')
ylabel('Peak Upstroke Velocity (mV/ms)')
box off

p1(1).Color = [c{1}, transparency_param];
p1(2).Color = [c{1}, transparency_param];
p1(3).Color = [c{2}, transparency_param];
p1(4).Color = [c{2}, transparency_param];
p1(6).Color = [c{1}, transparency_param];
p1(7).Color = [c{1}, transparency_param];
p1(8).Color = [c{2}, transparency_param];
p1(9).Color = [c{2}, transparency_param];


%% Restitution Curve comparison for kon = 1e1 (APD vs. BCL)
figure(5)
p1 = plot(GRFI_1e1(7:end,1), GRFI_1e1(7:end,4), '-', GRFN_1e1(7:end,1), GRFN_1e1(7:end,4), '--', ...
    GIFI_1e1(7:end,1), GIFI_1e1(7:end,4), '-', GIFN_1e1(7:end,1), GIFN_1e1(7:end,4), '--',...
    nodrug(7:end,1), nodrug(7:end,4), '-k',...
    GRHHI_1e1(7:end,1), GRHHI_1e1(7:end,4), '-', GRHHN_1e1(7:end,1), GRHHN_1e1(7:end,4), '--', ...
    GIHHI_1e1(7:end,1), GIHHI_1e1(7:end,4), '-', GIHHN_1e1(7:end,1), GIHHN_1e1(7:end,4), '--',...
    GRFI_1e1_BL(1,1), GRFI_1e1_BL(1,4), 'O', GRFN_1e1_BL(1,1), GRFN_1e1_BL(1,4), '.', ...
    GIFI_1e1_BL(1,1), GIFI_1e1_BL(1,4), 'O', GIFN_1e1_BL(1,1), GIFN_1e1_BL(1,4), '.',...
    GRHHI_1e1_BL(1,1), GRHHI_1e1_BL(1,4), 'O', GRHHN_1e1_BL(1,1), GRHHN_1e1_BL(1,4), '.', ...
    GIHHI_1e1_BL(1,1), GIHHI_1e1_BL(1,4), 'O', GIHHN_1e1_BL(1,1), GIHHN_1e1_BL(1,4), '.',...
    nodrug_BL(1,1), nodrug_BL(1,4), 'kO', 'MarkerSize', 3);
xlim(xlim_vec);
ylim([200, 320]);
set(gca,'XTick',x_tick_vals);
set(gca,'YTick',200:25:300, 'YTicklabels', []);
box off

p1(1).Color = [c{1}, transparency_param];
p1(2).Color = [c{1}, transparency_param];
p1(3).Color = [c{2}, transparency_param];
p1(4).Color = [c{2}, transparency_param];
p1(6).Color = c{1};
p1(7).Color = c{1};
p1(8).Color = c{2};
p1(9).Color = c{2};

p1(10).Color = [c{1}, transparency_param];
p1(11).Color = [c{1}, transparency_param];
p1(12).Color = [c{2}, transparency_param];
p1(13).Color = [c{2}, transparency_param];
p1(14).Color = c{1};
p1(15).Color = c{1};
p1(16).Color = c{2};
p1(17).Color = c{2};

max_APD_diff_1e1 = max([max(abs(GRFI_1e1(7:end,4) - nodrug(7:end,4)));...
    max(abs(GRFN_1e1(7:end,4) - nodrug(7:end,4)));...
    max(abs(GIFI_1e1(7:end,4) - nodrug(7:end,4)));...
    max(abs(GIFN_1e1(7:end,4) - nodrug(7:end,4)));...
    max(abs(GRHHI_1e1(7:end,4) - nodrug(7:end,4)));...
    max(abs(GRHHN_1e1(7:end,4) - nodrug(7:end,4)));...
    max(abs(GIHHI_1e1(7:end,4) - nodrug(7:end,4)));...
    max(abs(GIHHN_1e1(7:end,4) - nodrug(7:end,4)))])

%title('k_on = 1e1')
%legend('GRI', 'GRN', 'GII', 'GIN', 'no drug');

%%
figure(6)
p1 = plot(GRFI_1e3(7:end,1), GRFI_1e3(7:end,4), '-', GRFN_1e3(7:end,1), GRFN_1e3(7:end,4), '--', ...
    GIFI_1e3(7:end,1), GIFI_1e3(7:end,4), '-', GIFN_1e3(7:end,1), GIFN_1e3(7:end,4), '--',...
    nodrug(7:end,1), nodrug(7:end,4), 'k-',...
    GRHHI_1e3(7:end,1), GRHHI_1e3(7:end,4), '-', GRHHN_1e3(7:end,1), GRHHN_1e3(7:end,4), '--', ...
    GIHHI_1e3(7:end,1), GIHHI_1e3(7:end,4), '-', GIHHN_1e3(7:end,1), GIHHN_1e3(7:end,4), '--',...
    GRFI_1e3_BL(1,1), GRFI_1e3_BL(1,4), 'O', GRFN_1e3_BL(1,1), GRFN_1e3_BL(1,4), '.', ...
    GIFI_1e3_BL(1,1), GIFI_1e3_BL(1,4), 'O', GIFN_1e3_BL(1,1), GIFN_1e3_BL(1,4), '.',...
    GRHHI_1e3_BL(1,1), GRHHI_1e3_BL(1,4), 'O', GRHHN_1e3_BL(1,1), GRHHN_1e3_BL(1,4), '.', ...
    GIHHI_1e3_BL(1,1), GIHHI_1e3_BL(1,4), 'O', GIHHN_1e3_BL(1,1), GIHHN_1e3_BL(1,4), '.',...
    nodrug_BL(1,1), nodrug_BL(1,4), 'kO', 'MarkerSize', 3);
xlim(xlim_vec);
ylim([200, 320]);
set(gca,'XTick',x_tick_vals);
set(gca,'YTick',200:25:300, 'YTicklabels', []);
box off

p1(1).Color = [c{1}, transparency_param];
p1(2).Color = [c{1}, transparency_param];
p1(3).Color = [c{2}, transparency_param];
p1(4).Color = [c{2}, transparency_param];
p1(6).Color = c{1};
p1(7).Color = c{1};
p1(8).Color = c{2};
p1(9).Color = c{2};

p1(10).Color = [c{1}, transparency_param];
p1(11).Color = [c{1}, transparency_param];
p1(12).Color = [c{2}, transparency_param];
p1(13).Color = [c{2}, transparency_param];
p1(14).Color = c{1};
p1(15).Color = c{1};
p1(16).Color = c{2};
p1(17).Color = c{2};

max_APD_diff_1e3 = max([max(abs(GRFI_1e3(7:end,4) - nodrug(7:end,4)));...
    max(abs(GRFN_1e3(7:end,4) - nodrug(7:end,4)));...
    max(abs(GIFI_1e3(7:end,4) - nodrug(7:end,4)));...
    max(abs(GIFN_1e3(7:end,4) - nodrug(7:end,4)));...
    max(abs(GRHHI_1e3(7:end,4) - nodrug(7:end,4)));...
    max(abs(GRHHN_1e3(7:end,4) - nodrug(7:end,4)));...
    max(abs(GIHHI_1e3(7:end,4) - nodrug(7:end,4)));...
    max(abs(GIHHN_1e3(7:end,4) - nodrug(7:end,4)))])

%title('k_on = 1e3')
%legend('GRI', 'GRN', 'GII', 'GIN', 'no drug');

%%

figure(7)
p1 = plot(GRFI_1e5(7:end,1), GRFI_1e5(7:end,4), '-', GRFN_1e5(7:end,1), GRFN_1e5(7:end,4), '--', ...
    GIFI_1e5(7:end,1), GIFI_1e5(7:end,4), '-', GIFN_1e5(7:end,1), GIFN_1e5(7:end,4), '--',...
    nodrug(7:end,1), nodrug(7:end,4), 'k-',...
    GRHHI_1e5(7:end,1), GRHHI_1e5(7:end,4), '-', GRHHN_1e5(7:end,1), GRHHN_1e5(7:end,4), '--', ...
    GIHHI_1e5(7:end,1), GIHHI_1e5(7:end,4), '-', GIHHN_1e5(7:end,1), GIHHN_1e5(7:end,4), '--',...
    GRFI_1e5_BL(1,1), GRFI_1e5_BL(1,4), 'O', GRFN_1e5_BL(1,1), GRFN_1e5_BL(1,4), '.', ...
    GIFI_1e5_BL(1,1), GIFI_1e5_BL(1,4), 'O', GIFN_1e5_BL(1,1), GIFN_1e5_BL(1,4), '.',...
    GRHHI_1e5_BL(1,1), GRHHI_1e5_BL(1,4), 'O', GRHHN_1e5_BL(1,1), GRHHN_1e5_BL(1,4), '.', ...
    GIHHI_1e5_BL(1,1), GIHHI_1e5_BL(1,4), 'O', GIHHN_1e5_BL(1,1), GIHHN_1e5_BL(1,4), '.',...
    nodrug_BL(1,1), nodrug_BL(1,4), 'kO', 'MarkerSize', 3);
xlim(xlim_vec);
ylim([200, 320]);
set(gca,'XTick',x_tick_vals);
set(gca,'YTick',200:25:300, 'YTicklabels', []);
box off

p1(1).Color = [c{1}, transparency_param];
p1(2).Color = [c{1}, transparency_param];
p1(3).Color = [c{2}, transparency_param];
p1(4).Color = [c{2}, transparency_param];
p1(6).Color = c{1};
p1(7).Color = c{1};
p1(8).Color = c{2};
p1(9).Color = c{2};

p1(10).Color = [c{1}, transparency_param];
p1(11).Color = [c{1}, transparency_param];
p1(12).Color = [c{2}, transparency_param];
p1(13).Color = [c{2}, transparency_param];
p1(14).Color = c{1};
p1(15).Color = c{1};
p1(16).Color = c{2};
p1(17).Color = c{2};

%title('k_on = 1e5')
%legend('GRI', 'GRN', 'GII', 'GIN', 'no drug');

max_APD_diff_1e5 = max([max(abs(GRFI_1e5(7:end,4) - nodrug(7:end,4)));...
    max(abs(GRFN_1e5(7:end,4) - nodrug(7:end,4)));...
    max(abs(GIFI_1e5(7:end,4) - nodrug(7:end,4)));...
    max(abs(GIFN_1e5(7:end,4) - nodrug(7:end,4)));...
    max(abs(GRHHI_1e5(7:end,4) - nodrug(7:end,4)));...
    max(abs(GRHHN_1e5(7:end,4) - nodrug(7:end,4)));...
    max(abs(GIHHI_1e5(7:end,4) - nodrug(7:end,4)));...
    max(abs(GIHHN_1e5(7:end,4) - nodrug(7:end,4)))])
