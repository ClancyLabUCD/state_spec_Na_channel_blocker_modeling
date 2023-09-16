%Steffen Docken (Lewis Lab). All rights reserved.
%  Published in the Journal of Theoretical Biology under the title "Rate-dependent effects of state-specific sodium channel blockers in cardiac tissue: Insights from idealized models"
%8-1-16
%Function to calculate matrix with eigenvalues raised as 
%exponentials (D_mat), matrix of eigen vectors (E), inverse of 
%matrix of eigenvectors (E_inv), and vector of steady state values of 
%variables for gate immobilization and guarded receptor models with Ten 
%Tusscher model

%Derivation for these functions is in notes from around 6-30-16.  

function [D_mat, E, E_inv, X_ss] = Model_Matrices(V, k_on, k_off_0, D, ...
    inact_noninact, gi_gr)

%inact_noninact = 0 means drug binds to the inactive state (1 means the
%non-inactiveated state).
%gi_gr = 0 means the drug binds via the gate immobilization model (1 means
%the guarded receptor model).

R = 8314.472;	% J/mol*K
T = 310;% K
F = 96485.3415;

k_off = k_off_0*exp(-0.7*V*F/(R*T));

if (inact_noninact == 0) %inactive state
    if (gi_gr == 0) %gate immobilization model
        
        a_h = a_h_6_14_2016(V,T);
        b_h = b_h_6_14_2016(V,T);
        
        a = (-(a_h + b_h));
        b = (-a_h);
        c = (-D*k_on);
        d = (-(D*k_on + k_off));

        l_1 = (a + d + sqrt((a + d)^2 -4*(a*d - ...
            b*c)))/2;

        l_2 = (a + d - sqrt((a + d)^2 -4*(a*d - ...
            b*c)))/2;

        D_mat = [exp(l_1), 0;
            0, exp(l_2)];

        E = [1, 1;
            (l_1 - a)/b, (l_2 - a)/b];

        E_inv = [(l_2 - a)/(l_2 - l_1), ...
            -1/((l_2 - l_1)/b);
            (-l_1 + a)/(l_2 - l_1), ...
            1/((l_2 - l_1)/b)];

        X_ss = 1/(a*d - b*c)*[d*b - b*c;
            -b*c + a*c];
        
    elseif (gi_gr == 1) %guarded receptor model
    
        a_h = a_h_6_14_2016(V,T);
        b_h = b_h_6_14_2016(V,T);
        
        a = -(a_h + b_h);
        b = k_off;
        c = -(D*k_on + k_off);
        d = -a_h;
        f = D*k_on;

        l_1 = (c + a + sqrt(a^2 - 2*a*c + c^2 + 4*c*d))/2;

        l_2 = (c + a - sqrt(a^2 - 2*a*c + c^2 + 4*c*d))/2;

        l_3 = a;

        D_mat = [exp(l_1), 0, 0;
            0, exp(l_2), 0;
            0, 0, exp(l_3)];

        E = [0, 0, -c/b;
            c/(l_1-c), c/(l_2 - c), 0;
            1, 1, 1];

        E_inv = [(l_1 - c)*b/(c*sqrt(a^2 - 2*a*c + c^2 + 4*c*d)), ...
            (l_1 - c)*(c - l_2)/(c*sqrt(a^2 - 2*a*c + c^2 + 4*c*d)), ...
            (l_1 - c)/(sqrt(a^2 - 2*a*c + c^2 + 4*c*d));
            (c - l_2)*b/(c*sqrt(a^2 - 2*a*c + c^2 + 4*c*d)), ...
            (l_1 - c)*(l_2 - c)/(c*sqrt(a^2 - 2*a*c + c^2 + 4*c*d)),...
            (c - l_2)/(sqrt(a^2 - 2*a*c + c^2 + 4*c*d));
            -b/c, 0, 0];

        M_inv = [1/a, 0, 0;
            -b/(c*(a - d)), a/(c*(a - d)), -1/(a-d);
            b*d/(a*c*(a-d)), -d/(c*(a-d)), 1/(a-d)]; %inverse of matrix for differential equation

        X_ss = M_inv*[d; -f; d];
        %keyboard

        
    end
        
elseif (inact_noninact == 1) %non-inactive state binding
    if (gi_gr == 0) %gate immobilization model
        a_h = a_h_6_14_2016(V,T);
        b_h = b_h_6_14_2016(V,T);
        
        a = -(a_h + b_h + D*k_on);
        b = k_off-a_h;
        c = D*k_on;
        d = -k_off;
        f = -a_h;

        l_1 = (a + d + sqrt((a + d)^2 -4*(a*d - ...
            b*c)))/2;

        l_2 = (a + d - sqrt((a + d)^2 -4*(a*d - ...
            b*c)))/2;

        D_mat = [exp(l_1), 0;
            0, exp(l_2)];
        
        E = [1, 1;
            (l_1 - a)/b, (l_2 - a)/b];
    
        E_inv = [(l_2 - a)/(l_2 - l_1), ...
            -1/((l_2 - l_1)/b);
            (-l_1 + a)/(l_2 - l_1), ...
            1/((l_2 - l_1)/b)];
        
        X_ss = 1/(a*d - b*c)*[f*d;
            -f*c];
        
    elseif (gi_gr == 1) %guarded receptor model
        a_h = a_h_6_14_2016(V,T);
        b_h = b_h_6_14_2016(V,T);
        
        a = -(a_h + b_h);
        b = k_off;
        c = D*k_on + k_off;
        d = -a_h;

        l_1 = a;

        l_2 = (a - c + sqrt(a^2 - 2*a*c + c^2 + 4*c*d))/2;

        l_3 = (a - c - sqrt(a^2 - 2*a*c + c^2 + 4*c*d))/2;

        f = sqrt(a^2 - 2*a*c + c^2 + 4*c*d); %a factor used a lot in E_inv

        D_mat = [exp(l_1), 0, 0;
            0, exp(l_2), 0;
            0, 0, exp(l_3)];

        E = [c/b, 0, 0;
            0, c/l_2, c/l_3;
            1, 1, 1];

        E_inv = [b/c, 0, 0;
            -l_2*b/(c*f), -l_2*l_3/(c*f), l_2/f;
            l_3*b/(c*f), l_2*l_3/(c*f), -l_3/f];

        M_inv = [1/a, 0, 0;
            -b/(c*d), -(a-c)/(c*d), 1/d;
            b/(a*c), 1/c, 0]; %inverse of matrix for differential equation

        X_ss = -M_inv*[-d; 0; -d];
        %keyboard
        
    end
end
