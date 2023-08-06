%Steffen Docken (Lewis Lab). All rights reserved.
%  Published in the Journal of Theoretical Biology under the title "Rate-dependent effects of state-specific sodium channel blockers in cardiac tissue: Insights from idealized models"
%8-23-16
%Code to calculate the analytic solution for the guarded receptor model
%(HH Formulation) with the ten Tusscher model for
%the fast inactivation gate (h).
%This code also calculates the analytic solution of a non-state specific 
%binding drug, and a gate immobilization drug using the HH formulation.
%This code assumes that V is constant

%V is the current transmembrane potential, dt is the size of the time step,
%h_0 and b_0 are the initial h and b values, k_on is the binding rate,
%k_off_0 is the unbinding rate at V = 0mV, Drug is the concentration of
%drug (in M), inact_noninact is 0 for inactivated state binding, and 1 for
%non-inactivated state binding, GIHH_GRHH_nons = 0 means drug binds via the gate immoilization model (which is 
%modeled using the HH formulation), 1 means the drug binds via the guarded receptor model (which is
%is modeled using the HH formulation), and 2 means the drug is not state
%specific.

function [h_new, b_new] = HH_analytic(V, dt, h_0, b_0, k_on, ...
    k_off_0, Drug, inact_noninact, GIHH_GRHH_nons)

R = 8314.472;	% J/mol*K
T = 310;% K
F = 96485.3415;

k_off = @(V)(k_off_0*exp(-0.7*V*F/(R*T)));

a_h = @(V)a_h_6_14_2016(V,T);
b_h = @(V)b_h_6_14_2016(V,T);

tau_h_GRHH = @(V)(1/(a_h(V)+b_h(V))); %time constant of inactivation for the
%guarded receptor model in HH formulation
h_infty_GRHH = @(V)(a_h(V)/(a_h(V)+b_h(V))); %steady state fraction of non-
%inactivated channels in the guarded receptor model

tau_b_GRHH = @(V,D)(1/(D*k_on + k_off(V))); %time constant of drug binding 
%for a non-state specific drug
b_infty_GRHH = @(V,D)(D*k_on/(D*k_on + k_off(V))); %steady state fraction of 
%channels bound to drug for a non-state specific drug (or guarded receptor)

b_infty_GIHH_inact = @(V, D)( b_infty_GRHH(V,D)*(1 - h_infty_GRHH(V)))/( 1 - ...
    h_infty_GRHH(V)*b_infty_GRHH(V,D));%steady state fraction of 
%channels bound to drug for a gate immobilization with inactivated state
%binding drug
        
tau_b_GIHH_inact = @(V,D)(tau_b_GRHH(V,D)./( 1 - ...
    h_infty_GRHH(V)*b_infty_GRHH(V,D)));
%time constant for drug binding in gate immobilization with inactivated
%state binding model in HH formulation

b_infty_GIHH_noninact = @(V, D)(  b_infty_GRHH(V,D)*h_infty_GRHH(V))/( 1 - ...
    (1-h_infty_GRHH(V))*b_infty_GRHH(V,D));%steady state fraction of 
%channels bound to drug for a gate immobilization with noninactivated state
%binding drug

tau_b_GIHH_noninact = @(V,D)(tau_b_GRHH(V,D)./( 1 - ...
    (1 - h_infty_GRHH(V))*b_infty_GRHH(V,D)));
%time constant for drug binding in gate immobilization with noninactivated
%state binding model in HH formulation

h_new = h_infty_GRHH(V) - (h_infty_GRHH(V) - h_0)*exp(-dt/...
    tau_h_GRHH(V));

if (GIHH_GRHH_nons == 0) %gate immobilization model (HH formulation)
    %h_new = h_infty_GRHH(V); %for this formulation of the gate immobilization
    %model, h is just at it's steady state value.
    
    if (inact_noninact == 0) %inactive state block
        
        b_new = b_infty_GIHH_inact(V, Drug) - (b_infty_GIHH_inact(V, Drug) -...
            b_0)*exp(-dt/tau_b_GIHH_inact(V, Drug));
        
    elseif (inact_noninact == 1) %closed state block
        
        b_new = b_infty_GIHH_noninact(V, Drug) - (b_infty_GIHH_noninact(V, Drug) -...
            b_0)*exp(-dt/tau_b_GIHH_noninact(V, Drug));
    end
    
elseif (GIHH_GRHH_nons == 1) %guarded receptor model (HH formulation)
    if (inact_noninact == 0) %inactive state

        b_new = b_infty_GRHH(V,Drug) - (b_infty_GRHH(V,Drug) - b_0)*...
            exp(-1/tau_b_GRHH(V,Drug)*((1-h_infty_GRHH(V))*dt + tau_h_GRHH(V)*...
            (h_infty_GRHH(V) - h_0)*(1 - exp(-dt/tau_h_GRHH(V)))));
        
    elseif (inact_noninact == 1) %closed state block
        b_new = b_infty_GRHH(V,Drug) - (b_infty_GRHH(V,Drug) - b_0)*...
            exp(-1/tau_b_GRHH(V,Drug)*(h_infty_GRHH(V)*dt - (h_infty_GRHH(V) - ...
            h_0)*tau_h_GRHH(V)*(1 - exp(-dt/tau_h_GRHH(V)))));
    end
elseif (GIHH_GRHH_nons ==2) %%non-state-specific binding
    b_new = b_infty_GRHH(V,Drug) - (b_infty_GRHH(V,Drug) - b_0)*...
        exp(-dt/tau_b_GRHH(V,Drug));
end
