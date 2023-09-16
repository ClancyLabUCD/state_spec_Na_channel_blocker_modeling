%Steffen Docken (Lewis Lab). All rights reserved.
%  Published in the Journal of Theoretical Biology under the title "Rate-dependent effects of state-specific sodium channel blockers in cardiac tissue: Insights from idealized models"
%9-3-16
%Function for beta_h for Simple Na-channel model at 37.  The Q10 factor alters
%the rates from what they are at 310K

function b_h = b_h_6_14_2016(V, Temp)

Q10 = 3; %Q10 factor for change of temp from body temp

b_h = 14.15*exp(V/14.91)*Q10^((Temp - 310)/10);
