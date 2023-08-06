%Steffen Docken (Lewis Lab). All rights reserved.
%  Published in the Journal of Theoretical Biology under the title "Rate-dependent effects of state-specific sodium channel blockers in cardiac tissue: Insights from idealized models"
%12-7-17
%Code to load and save data files as a large .mat file
%specifically Upstroke vs. Period data for three different tau_b values.
%k_D0 = 10e-6M, [D] = 20e-6M, and the stimulus lasts 0.6ms

oldFolder = cd('Ten_Tusscher_model_with_drug_output/Upstroke_vs_Period_3_tau_b_12_5_17/');

nodrug_BL = load('Simple_nodrug_drug_Upstroke_rest_30s_stimdur06.txt');

GRHHI_1e5_BL = load('Simple_HH_guarded_receptor_inactive_drug_Upstroke_rest_30s_Diffval1e5_kd10e-6_Drug20e-6_stimdur06.txt');
GRHHI_1e3_BL = load('Simple_HH_guarded_receptor_inactive_drug_Upstroke_rest_30s_Diffval1e3_kd10e-6_Drug20e-6_stimdur06.txt');
GRHHI_1e1_BL = load('Simple_HH_guarded_receptor_inactive_drug_Upstroke_rest_30s_Diffval1e1_kd10e-6_Drug20e-6_stimdur06.txt');

GRHHN_1e5_BL = load('Simple_HH_guarded_receptor_noninactive_drug_Upstroke_rest_30s_Diffval1e5_kd10e-6_Drug20e-6_stimdur06.txt');
GRHHN_1e3_BL = load('Simple_HH_guarded_receptor_noninactive_drug_Upstroke_rest_30s_Diffval1e3_kd10e-6_Drug20e-6_stimdur06.txt');
GRHHN_1e1_BL = load('Simple_HH_guarded_receptor_noninactive_drug_Upstroke_rest_30s_Diffval1e1_kd10e-6_Drug20e-6_stimdur06.txt');

GIHHI_1e5_BL = load('Simple_HH_gate_immobilization_inactive_drug_Upstroke_rest_30s_Diffval1e5_kd10e-6_Drug20e-6_stimdur06.txt');
GIHHI_1e3_BL = load('Simple_HH_gate_immobilization_inactive_drug_Upstroke_rest_30s_Diffval1e3_kd10e-6_Drug20e-6_stimdur06.txt');
GIHHI_1e1_BL = load('Simple_HH_gate_immobilization_inactive_drug_Upstroke_rest_30s_Diffval1e1_kd10e-6_Drug20e-6_stimdur06.txt');

GIHHN_1e5_BL = load('Simple_HH_gate_immobilization_noninactive_drug_Upstroke_rest_30s_Diffval1e5_kd10e-6_Drug20e-6_stimdur06.txt');
GIHHN_1e3_BL = load('Simple_HH_gate_immobilization_noninactive_drug_Upstroke_rest_30s_Diffval1e3_kd10e-6_Drug20e-6_stimdur06.txt');
GIHHN_1e1_BL = load('Simple_HH_gate_immobilization_noninactive_drug_Upstroke_rest_30s_Diffval1e1_kd10e-6_Drug20e-6_stimdur06.txt');

GRFI_1e5_BL = load('Simple_Full_guarded_receptor_inactive_drug_Upstroke_rest_30s_Diffval1e5_kd10e-6_Drug20e-6_stimdur06.txt');
GRFI_1e3_BL = load('Simple_Full_guarded_receptor_inactive_drug_Upstroke_rest_30s_Diffval1e3_kd10e-6_Drug20e-6_stimdur06.txt');
GRFI_1e1_BL = load('Simple_Full_guarded_receptor_inactive_drug_Upstroke_rest_30s_Diffval1e1_kd10e-6_Drug20e-6_stimdur06.txt');

GRFN_1e5_BL = load('Simple_Full_guarded_receptor_noninactive_drug_Upstroke_rest_30s_Diffval1e5_kd10e-6_Drug20e-6_stimdur06.txt');
GRFN_1e3_BL = load('Simple_Full_guarded_receptor_noninactive_drug_Upstroke_rest_30s_Diffval1e3_kd10e-6_Drug20e-6_stimdur06.txt');
GRFN_1e1_BL = load('Simple_Full_guarded_receptor_noninactive_drug_Upstroke_rest_30s_Diffval1e1_kd10e-6_Drug20e-6_stimdur06.txt');

GIFI_1e5_BL = load('Simple_Full_gate_immobilization_inactive_drug_Upstroke_rest_30s_Diffval1e5_kd10e-6_Drug20e-6_stimdur06.txt');
GIFI_1e3_BL = load('Simple_Full_gate_immobilization_inactive_drug_Upstroke_rest_30s_Diffval1e3_kd10e-6_Drug20e-6_stimdur06.txt');
GIFI_1e1_BL = load('Simple_Full_gate_immobilization_inactive_drug_Upstroke_rest_30s_Diffval1e1_kd10e-6_Drug20e-6_stimdur06.txt');

GIFN_1e5_BL = load('Simple_Full_gate_immobilization_noninactive_drug_Upstroke_rest_30s_Diffval1e5_kd10e-6_Drug20e-6_stimdur06.txt');
GIFN_1e3_BL = load('Simple_Full_gate_immobilization_noninactive_drug_Upstroke_rest_30s_Diffval1e3_kd10e-6_Drug20e-6_stimdur06.txt');
GIFN_1e1_BL = load('Simple_Full_gate_immobilization_noninactive_drug_Upstroke_rest_30s_Diffval1e1_kd10e-6_Drug20e-6_stimdur06.txt');


cd(oldFolder);

save('Upstroke_rest_30s_3_tau_b.mat', 'nodrug_BL',...
    'GRHHI_1e5_BL', 'GRHHI_1e3_BL', 'GRHHI_1e1_BL',...
    'GRHHN_1e5_BL', 'GRHHN_1e3_BL', 'GRHHN_1e1_BL',...
    'GIHHI_1e5_BL', 'GIHHI_1e3_BL', 'GIHHI_1e1_BL',...
    'GIHHN_1e5_BL', 'GIHHN_1e3_BL', 'GIHHN_1e1_BL',...
    'GRFI_1e5_BL', 'GRFI_1e3_BL', 'GRFI_1e1_BL',...
    'GRFN_1e5_BL', 'GRFN_1e3_BL', 'GRFN_1e1_BL',...
    'GIFI_1e5_BL', 'GIFI_1e3_BL', 'GIFI_1e1_BL',...
    'GIFN_1e5_BL', 'GIFN_1e3_BL', 'GIFN_1e1_BL');
