%Steffen Docken (Lewis Lab). All rights reserved.
%  Published in the Journal of Theoretical Biology under the title "Rate-dependent effects of state-specific sodium channel blockers in cardiac tissue: Insights from idealized models"
%1-30-17
%Code to calculate instantaneous tau and steady state (or "target") values
%of different drug models and formulations when the dynamics of I_O are
%shifted into the form dI_O/dt =  1/tau(V,h,b) (I_O_infty(V,h,b) - I_O)

%V is transmembrane potential, h_0 is current value of h, b_0 is current
%value of b, Drug is concentration of drug, k_on is the drug binding rate,
%k_off_0 is the k_off value at V = 0, Na_Model is the model used for the Na
%inactivation gate, 
%gi_gr = 0 means the drug binds via the gate immobilization model (1 means
%the guarded receptor model).
%inact_noninact = 0 means drug binds to the inactive state (1 means the
%non-inactiveated state).
%HH_F = 0 means the HH Formulation (1 means the full formulation)


function tau_I_O_approx = I_O_approx_dynamics(V, h_0,...
    Drug, k_on, k_off_0, inact_noninact, gi_gr, HH_F)

R = 8314.472;% J/mol*K
T = 310;% K
F = 96485.3415;

%keyboard

k_off = k_off_0*exp(-0.7*V*F/(R*T));

K_D = k_off/k_on;

%Model from 2-2-2016
a_h = @(V)a_h_6_14_2016(V, T);
b_h = @(V)b_h_6_14_2016(V, T);

tau_h = 1/(a_h(V) + b_h(V));
tau_b = 1/(Drug*k_on + k_off);%time constants of h and drug binding (in the 
%absence of state-specificity)

h_infty = a_h(V)/(a_h(V) + b_h(V));
b_infty = Drug*k_on/(Drug*k_on + k_off); %steady state of h and b (in the 
%absence of state-specificity)

if (gi_gr == 0)
    if (inact_noninact == 0) %Gate Immobilization inactivated state binding 
        if (HH_F == 0) %HH Formulation
            tau_I_O_approx = (tau_h*tau_b*(1+K_D/Drug))/...
                ((1-h_0)*tau_h + tau_b*(1+K_D/Drug));
            
            
        elseif (HH_F == 1) %Full Formulation
            tau_I_O_approx = tau_h;
            
        end
        
    elseif (inact_noninact == 1) %Gate Immobilization non-inactivated state binding
        if (HH_F == 0) %HH Formulation
            tau_I_O_approx = (tau_h*tau_b*(1+K_D/Drug))/...
                (h_0*tau_h + tau_b*(1+K_D/Drug));
            
        elseif (HH_F == 1) %Full Formulation
            tau_I_O_approx = (tau_h*tau_b*(1+K_D/Drug))/...
                (tau_h + tau_b*(1+K_D/Drug));
            
        end
        
    end
    
elseif (gi_gr == 1)
    if (inact_noninact == 0) %Guarded Receptor Inactivated state binding
        if (HH_F == 0) %HH Formulation
            tau_I_O_approx = tau_h*tau_b/(tau_b + (1 - h_0)*tau_h);
        
        elseif (HH_F == 1) %Full Formulation
            tau_I_O_approx = tau_h;
            
        end
        
    elseif (inact_noninact == 1) %Guarded Receptor Noninactivated state binding
        if (HH_F == 0) %HH Formulation
            tau_I_O_approx = tau_h*tau_b/(tau_b + h_0*tau_h);
            
        elseif (HH_F == 1)
            tau_I_O_approx = tau_h*tau_b/(tau_b + tau_h);
            
        end
    end
end
    
