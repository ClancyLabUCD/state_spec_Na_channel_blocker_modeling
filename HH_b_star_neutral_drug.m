%Steffen Docken (Lewis Lab). All rights reserved.
%  Published in the Journal of Theoretical Biology under the title "Rate-dependent effects of state-specific sodium channel blockers in cardiac tissue: Insights from idealized models"
%8-8-16
%Code to calculate steady state value for iterative map of b for the HH
%formulation of the guarded receptor and gate immobilization models with 
%the ten Tusscher model for the fast inactivation gate (h).
%This code also calculates the steady state value for the iterative map of
%a non-state specific binding drug.
%for neutral drug

function [h_ss,b_ss] = HH_b_star_neutral_drug(A, DI, V_A, V_DI, k_on, k_off_0, Drug,...
    inact_noninact, GIHH_GRHH_nons)

%inact_noninact = 0 means drug binds to the inactive state (1 means the
%non-inactiveated state).
%GIHH_GRHH_nons = 0 means drug binds via the gate immoilization model (which is 
%modeled using the HH formulation), 1 means the drug binds via the guarded receptor model (which is
%is modeled using the HH formulation), and 2 means the drug is not state
%specific.

R = 8314.472;	% J/mol*K
T = 310;% K
F = 96485.3415;

a_h = @(V)a_h_6_14_2016(V, T);
b_h = @(V)b_h_6_14_2016(V, T);

k_off = k_off_0;

tau_h_GRHH = @(V)(1/(a_h(V)+b_h(V))); %time constant of inactivation for the
%guarded receptor model in HH formulation
h_infty_GRHH = @(V)(a_h(V)/(a_h(V)+b_h(V)));%steady state fraction of non-
%inactivated channels in the guarded receptor model

tau_b_GRHH = @(D)(1/(D*k_on + k_off)); %time constant of drug binding 
%for a non-state specific drug
b_infty_GRHH = @(D)(D*k_on/(D*k_on + k_off)); %steady state fraction of 
%channels bound to drug for a non-state specific drug (or guarded receptor)

b_infty_GIHH_inact = @(V, D)( (1 - h_infty_GRHH(V))*D*k_on)/( (1 - ...
    h_infty_GRHH(V))*D*k_on + k_off); %steady state fraction of 
%channels bound to drug for a gate immobilization with inactivated state
%binding drug
        
tau_b_GIHH_inact = @(V,D)(1/( (1 - h_infty_GRHH(V))*D*k_on + k_off));%time constant for drug binding in gate immobilization with inactivated
%state binding model in HH formulation

b_infty_GIHH_noninact = @(V, D)( h_infty_GRHH(V)*D*k_on)/...
    (h_infty_GRHH(V)*D*k_on + k_off); %steady state fraction of 
%channels bound to drug for a gate immobilization with noninactivated state
%binding drug

tau_b_GIHH_noninact = @(V,D)(1/(h_infty_GRHH(V)*D*k_on + k_off));
%time constant for drug binding in gate immobilization with noninactivated
%state binding model in HH formulation

h_ss = (h_infty_GRHH(V_DI)-(h_infty_GRHH(V_DI)-h_infty_GRHH(V_A)*(1-exp(-A/tau_h_GRHH(V_A))))*...
    exp(-DI/tau_h_GRHH(V_DI)))/(1-exp(-A/tau_h_GRHH(V_A)-DI/tau_h_GRHH(V_DI)));
%steady state value of h before each stimulus.

if (GIHH_GRHH_nons == 0) %Gate immobilization model using HH 
    %formulation  %% Note that this assumes h is instantaneous and
    %immediately reaches its steady state value
    if (inact_noninact == 0) %inactive state
    
        b_ss = (b_infty_GIHH_inact(V_DI,Drug) - (b_infty_GIHH_inact(V_DI,Drug)...
            - (1 - exp(-A/tau_b_GIHH_inact(V_A,Drug)))*...
            b_infty_GIHH_inact(V_A,Drug))*exp(-DI/tau_b_GIHH_inact(V_DI,Drug)))/...
            (1 - exp(-A/tau_b_GIHH_inact(V_A,Drug) - ...
            DI/tau_b_GIHH_inact(V_DI,Drug)));
        
    elseif (inact_noninact == 1) %closed state block
         b_ss = (b_infty_GIHH_noninact(V_DI,Drug) - (b_infty_GIHH_noninact(V_DI,Drug)...
            - (1 - exp(-A/tau_b_GIHH_noninact(V_A,Drug)))*...
            b_infty_GIHH_noninact(V_A,Drug))*exp(-DI/tau_b_GIHH_noninact(V_DI,Drug)))/...
            (1 - exp(-A/tau_b_GIHH_noninact(V_A,Drug) - ...
            DI/tau_b_GIHH_noninact(V_DI,Drug))); 
    end

elseif (GIHH_GRHH_nons == 1) %guarded receptor model (HH formulation)
    if (inact_noninact == 0) %inactive state
    
        b_exp_A_inact = exp(-1/tau_b_GRHH(Drug)*((1-h_infty_GRHH(V_A))*A + ...
            (h_infty_GRHH(V_A)-h_ss)*tau_h_GRHH(V_A)*...
            (1-exp(-A/tau_h_GRHH(V_A)))));

        b_exp_DI_inact = exp(-1/tau_b_GRHH(Drug)*((1-h_infty_GRHH(V_DI))*DI +...
            (h_infty_GRHH(V_DI) - h_infty_GRHH(V_A) +(h_infty_GRHH(V_A) - h_ss)*...
            exp(-A/tau_h_GRHH(V_A)))*tau_h_GRHH(V_DI)*(1-exp(-DI/tau_h_GRHH(V_DI)))));

        b_ss = (b_infty_GRHH(Drug) - (b_infty_GRHH(Drug)-b_infty_GRHH(Drug)*...
            (1-b_exp_A_inact))*b_exp_DI_inact)/...
            (1-b_exp_A_inact*b_exp_DI_inact);
        
    elseif (inact_noninact == 1) %non-inactive state block
         b_exp_A_noninact = exp(-1/tau_b_GRHH(Drug)*(h_infty_GRHH(V_A)*A -...
             (h_infty_GRHH(V_A) - h_ss)*tau_h_GRHH(V_A)*...
             (1-exp(-A/tau_h_GRHH(V_A)))));
         
         b_exp_DI_noninact = exp(-1/tau_b_GRHH(Drug)*(h_infty_GRHH(V_DI)*...
             DI - (h_infty_GRHH(V_DI) - h_infty_GRHH(V_A) + (h_infty_GRHH(V_A) - ...
             h_ss)*exp(-A/tau_h_GRHH(V_A)))*tau_h_GRHH(V_DI)*(1 - ...
             exp(-DI/tau_h_GRHH(V_DI)))));
         
         b_ss = (b_infty_GRHH(Drug) - (b_infty_GRHH(Drug) - ...
             b_infty_GRHH(Drug)*(1 - b_exp_A_noninact))*b_exp_DI_noninact)/...
             (1 - b_exp_A_noninact*b_exp_DI_noninact);   
    end
    
elseif (GIHH_GRHH_nons == 2) %non-state-specific binding
    b_ss = (b_infty_GRHH(Drug) - (b_infty_GRHH(Drug) - (1 -...
        exp(-A/tau_b_GRHH(Drug)))*b_infty_GRHH(Drug))*exp(-DI/...
        tau_b_GRHH(Drug)))/(1 - exp(-A/tau_b_GRHH(Drug) - ...
        DI/tau_b_GRHH(Drug)));
    
end

        
