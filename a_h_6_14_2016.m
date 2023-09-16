%Steffen Docken (Lewis Lab). All rights reserved.
%  Published in the Journal of Theoretical Biology under the title "Rate-dependent effects of state-specific sodium channel blockers in cardiac tissue: Insights from idealized models"
%6-7-17
%Function for alpha_h for Simple Na-channel model.  The Q10 factor alters
%the rates from what they are at 310K

function a_h = a_h_6_14_2016(V, Temp)

Q10 = 3; %Q10 factor for change of temp from body temp

a_h = 6.169e-05*exp(V/-9.328)*Q10^((Temp - 310)/10);
