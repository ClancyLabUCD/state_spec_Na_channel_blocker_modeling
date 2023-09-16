%Steffen Docken (Lewis Lab). All rights reserved.
%  Published in the Journal of Theoretical Biology under the title "Rate-dependent effects of state-specific sodium channel blockers in cardiac tissue: Insights from idealized models"
%1-18-17
%Code to calculate instantaneous steady state values and time constants of 
%different HH Formulation drug models

%V is the current transmembrane potential, h_0 is the current value of h,
%k_on is the binding rate,
%k_off_0 is the unbinding rate at V = 0mV, Drug is the concentration of
%drug (in M), inact_noninact is 0 for inactivated state binding, and 1 for
%non-inactivated state binding, GIHH_GRHH_nons = 0 means drug binds via the gate immoilization model (which is 
%modeled using the HH formulation), 1 means the drug binds via the guarded receptor model (which is
%is modeled using the HH formulation), and 2 means the drug is not state
%specific, Na_model is the Na-channel model being used

function [h_infty, b_infty, tau_h, tau_b] = HH_infty_tau(V, h_0, k_on,...
    k_off_0, Drug, inact_noninact, GIHH_GRHH_nons)

R = 8314.472;	% J/mol*K
T = 310;% K
F = 96485.3415;

k_off = @(V)(k_off_0*exp(-0.7*V*F/(R*T)));

a_h = @(V)a_h_6_14_2016(V, T);
b_h = @(V)b_h_6_14_2016(V, T);

h_infty = (a_h(V)/(a_h(V)+b_h(V)));

tau_h = (1/(a_h(V)+b_h(V))); %h dynamics are independent of model used

b_infty_nons = (Drug*k_on/(Drug*k_on + k_off(V)));

tau_b_nons = (1/(Drug*k_on + k_off(V))); %steady state and tau values for 
%nonstate-specific binding system

%HH formulations of the differential equations for the drug binding
%models can be reformulated as db/dt = 1/tau_b(b_infty - b) where tau_b and
%b_infty can be functions of h.  The resulting tau_b and b_infty are:

if (GIHH_GRHH_nons == 0) %Gate immobilization model
    if (inact_noninact == 0) %inactivated state binding
        
        b_infty = ((1 - h_0)*b_infty_nons)/(1 - ...
            h_0*b_infty_nons);
        
        tau_b = tau_b_nons/(1 - h_0*b_infty_nons);
        
    elseif (inact_noninact == 1) %noninactivated state binding
        
        b_infty = (b_infty_nons*h_0)/(1 - (1-h_0)*b_infty_nons);
        
        tau_b = tau_b_nons/(1 - (1-h_0)*b_infty_nons);
    end
    
elseif (GIHH_GRHH_nons == 1) %Guarded Receptor model
    if (inact_noninact == 0) %inactivated state binding
        
        b_infty = b_infty_nons;
        
        tau_b = tau_b_nons/(1-h_0);
        
    elseif (inact_noninact == 1)
        
        b_infty = b_infty_nons;
        
        tau_b = tau_b_nons/h_0;
        
    end
    
elseif (GIHH_GRHH_nons == 2) %nonstate specific binding
    
    b_infty = b_infty_nons;
    
    tau_b = tau_b_nons;
    
end
        
        
