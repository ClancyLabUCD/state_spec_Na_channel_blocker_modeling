%Steffen Docken (Lewis Lab). All rights reserved.
%  Published in the Journal of Theoretical Biology under the title "Rate-dependent effects of state-specific sodium channel blockers in cardiac tissue: Insights from idealized models"
%2-27-17
%Code to load and save data files as a large .mat file
%This data was generated with a stimulus that lasted for .6ms
%dt = 0.01 and the BCL was 1000ms

oldFolder = cd('Ten_Tusscher_model_with_drug_output/Upstroke_vs_tau_b_BCL1000_12_7_17/');

nodrug = load('Simple_nodrug_drug_Upstroke_BCL1000_Stimdur06.txt');

GRHHI = load('Simple_HH_guarded_receptor_inactive_drug_Upstroke_BCL1000_koninit_1e5_kd10e-6_Drug20e-6_Stimdur06.txt');

GRHHN = load('Simple_HH_guarded_receptor_noninactive_drug_Upstroke_BCL1000_koninit_1e5_kd10e-6_Drug20e-6_Stimdur06.txt');

GIHHI = load('Simple_HH_gate_immobilization_inactive_drug_Upstroke_BCL1000_koninit_1e5_kd10e-6_Drug20e-6_Stimdur06.txt');

GIHHN = load('Simple_HH_gate_immobilization_noninactive_drug_Upstroke_BCL1000_koninit_1e5_kd10e-6_Drug20e-6_Stimdur06.txt');

GRFI = load('Simple_Full_guarded_receptor_inactive_drug_Upstroke_BCL1000_koninit_1e5_kd10e-6_Drug20e-6_Stimdur06.txt');

GRFN = load('Simple_Full_guarded_receptor_noninactive_drug_Upstroke_BCL1000_koninit_1e5_kd10e-6_Drug20e-6_Stimdur06.txt');

GIFI = load('Simple_Full_gate_immobilization_inactive_drug_Upstroke_BCL1000_koninit_1e5_kd10e-6_Drug20e-6_Stimdur06.txt');

GIFN = load('Simple_Full_gate_immobilization_noninactive_drug_Upstroke_BCL1000_koninit_1e5_kd10e-6_Drug20e-6_Stimdur06.txt');

cd(oldFolder);

save('Upstroke_vs_taub_BCL1000_data_12_7_17.mat', 'nodrug',...
    'GRHHI', 'GRHHN', 'GIHHI', 'GIHHN', 'GRFI', 'GRFN', 'GIFI', 'GIFN');
