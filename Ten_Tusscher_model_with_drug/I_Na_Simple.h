//
//  I_Na_Simple.h
//  Ten_Tusscher_model_with_drug
//
//  Created by Steffen Docken on 12/5/17.
//  Copyright © 2017 Steffen Docken (Lewis Lab). All rights reserved.
//  Published in the Journal of Theoretical Biology under the title "Rate-dependent effects of state-specific sodium channel blockers in cardiac tissue: Insights from idealized models"
//

#ifndef I_Na_Simple_h
#define I_Na_Simple_h

#include <math.h>
#include "Global_variables.h"
#include "Numerical_Methods.h"

//My simple model of the sodium current.  These parameters were found on 6-14-16. h parameters come from taking the values from Moreno Methods paper and then fitting h_infty (as a function of 2 parameters) to SS Availability and then holding these parameters constant while tau_h is fitted to tau_50.  Then, holding the h parameters constant, the parameters for m were fit to SS Activation (here the full Na-channel model was used, not just m) and tau_m data. Parameters were then adjusted to 310K
void Calculate_I_Na_Simple(Cell_param *Cell_ptr, double t, double tTemp, double kon){
    
    double a_m, b_m, a_h, b_h;
    double tau_m, tau_h, m_inf, h_inf;
    const double G_Na=17.25; //nS/pF
    
    double E_Na=(R*T/F)*log(Na_out/Cell_ptr->Na_in);
    
    a_m = 45.43*exp((Cell_ptr->V)/13.78);
    b_m = 0.6628*exp((Cell_ptr->V)/-23.25);
    
    a_h = 6.169e-05*exp((Cell_ptr->V)/-9.328);
    b_h = 14.15*exp((Cell_ptr->V)/14.91); //m and h params adjusted to 310K.
    
    m_inf = a_m/(a_m + b_m);
    h_inf = a_h/(a_h + b_h);
    
    tau_m = 1.0/(a_m + b_m);
    tau_h = 1.0/(a_h + b_h);
    
    double k_on = kon;
    double k_off = kon*kd0_drug*exp(-dd_drug*(Cell_ptr->V)*F/(R*T));
    //V-dependent  kon*kd0_drug;//non V-dependent    
    
    double tau_b = 1.0/(Drug*k_on + k_off);
    double b_inf = Drug*k_on/(Drug*k_on + k_off);
    
    Cell_ptr->m = Rush_Larsen(Cell_ptr->m, m_inf, tau_m);
    
    switch (Bind_Scheme) {
        case 1: {//This is for the HH formulation of guarded receptor binding scheme with inactive state binding.
            Cell_ptr->h = Rush_Larsen(Cell_ptr->h, h_inf, tau_h);
            Cell_ptr->b = Rush_Larsen(Cell_ptr->b, b_inf, tau_b/(1.0-Cell_ptr->h));
            
            Cell_ptr->I_Na = G_Na*(Cell_ptr->V-E_Na)*(Cell_ptr->m)*(Cell_ptr->m)*(Cell_ptr->m)*(Cell_ptr->h)*(1-(Cell_ptr->b));
            
            break;
        }
        case 2: {//This is for the HH formulation of guarded receptor binding scheme with non-inactivated state state binding.
            Cell_ptr->h = Rush_Larsen(Cell_ptr->h, h_inf, tau_h);
            Cell_ptr->b = Rush_Larsen(Cell_ptr->b, b_inf, tau_b/(Cell_ptr->h));
            
            Cell_ptr->I_Na = G_Na*(Cell_ptr->V-E_Na)*(Cell_ptr->m)*(Cell_ptr->m)*(Cell_ptr->m)*(Cell_ptr->h)*(1-(Cell_ptr->b));
            
            break;
        }
        case 3: {//This is for the HH formulation of gate immobilization binding scheme with inact state state binding.
            double h_0 = Cell_ptr->h; //value of h from previous time step
            double b_inf_GISI = (1.0-h_0)*b_inf/(1.0 - h_0*b_inf); //steady state of b for HH Formulation of Gate Immobilization model with inactivated state binding (for the current value of h)
            double tau_b_GISI = tau_b/(1.0-h_0*b_inf); //time constant of b for HH Formulation of Gate Immobilization model with inactivated state binding (for the current value of h)
            
            Cell_ptr->h = Rush_Larsen(Cell_ptr->h, h_inf, tau_h);
            Cell_ptr->b = Rush_Larsen(Cell_ptr->b, b_inf_GISI, tau_b_GISI);
            
            Cell_ptr->I_Na = G_Na*(Cell_ptr->V-E_Na)*(Cell_ptr->m)*(Cell_ptr->m)*(Cell_ptr->m)*(Cell_ptr->h)*(1.0-(Cell_ptr->b));
            
            break;
        }
        case 4: {//This is for the HH formulation of gate immobilization binding scheme with non-inact state state binding.
            double h_0 = Cell_ptr->h; //value of h from previous time step
            double b_inf_GISN = h_0*b_inf/(1.0 - (1.0-h_0)*b_inf); //steady state of b for HH Formulation of Gate Immobilization model with non-inactivated state binding (for the current value of h)
            double tau_b_GISN = tau_b/(1.0-(1.0-h_0)*b_inf); //time constant of b for HH Formulation of Gate Immobilization model with non-inactivated state binding (for the current value of h)
            
            Cell_ptr->h = Rush_Larsen(Cell_ptr->h, h_inf, tau_h);
            Cell_ptr->b = Rush_Larsen(Cell_ptr->b, b_inf_GISN, tau_b_GISN);
            
            Cell_ptr->I_Na = G_Na*(Cell_ptr->V-E_Na)*(Cell_ptr->m)*(Cell_ptr->m)*(Cell_ptr->m)*(Cell_ptr->h)*(1.0-(Cell_ptr->b));
            
            break;
        }
        case 5: {//Full Guarded Receptor Inactive State
            double h_0 = Cell_ptr->h; //value of h (fraction of non-inactivated channels) from previous time step
            double b_0 = Cell_ptr->b; //value of b (fraction of channels bound to drug) from previous time step
            double i_O_0 = Cell_ptr->i_O; //value of i_O (fraction of channels that are both non-inactivatd and non-drug bound) from previous time step
            
            double b_inf_GRFI = (1.0-h_0)*b_inf + h_0 - i_O_0; //steady state of b for Full Formulation of Guarded Receptor model with inactivated state binding (for the current value of h)
            double i_O_inf = (1.0 - b_0)*h_inf;
            
            Cell_ptr->h = Rush_Larsen(Cell_ptr->h, h_inf, tau_h);
            Cell_ptr->b = Rush_Larsen(Cell_ptr->b, b_inf_GRFI, tau_b);
            Cell_ptr->i_O = Rush_Larsen(Cell_ptr->i_O, i_O_inf, tau_h);
            
            Cell_ptr->I_Na = G_Na*(Cell_ptr->V-E_Na)*(Cell_ptr->m)*(Cell_ptr->m)*(Cell_ptr->m)*(Cell_ptr->i_O);
            
            break;
        }
        case 6: {//Full Guarded Receptor Non-inactive State
            double h_0 = Cell_ptr->h; //value of h (fraction of non-inactivated channels) from previous time step
            double b_0 = Cell_ptr->b; //value of b (fraction of channels bound to drug) from previous time step
            double i_C_0 = 1.0 - Cell_ptr->b - Cell_ptr->i_O; //value of i_C (fraction of channels that are both inactivatd and non-drug bound) from previous time step
            
            double b_inf_GRFN = (1.0 - (1.0 - b_inf)*h_0 - i_C_0);  //steady state of b for Full Formulation of Guarded Receptor model with non-inactivated state binding (for the current value of h)
            double i_C_inf = (1.0 - h_inf)*(1.0 - b_0);
            
            Cell_ptr->h = Rush_Larsen(Cell_ptr->h, h_inf, tau_h);
            Cell_ptr->b = Rush_Larsen(Cell_ptr->b, b_inf_GRFN, tau_b);
            double i_C_new = Rush_Larsen(i_C_0, i_C_inf, tau_h);
            
            Cell_ptr->i_O = 1.0 - Cell_ptr->b - i_C_new;
            
            Cell_ptr->I_Na = G_Na*(Cell_ptr->V-E_Na)*(Cell_ptr->m)*(Cell_ptr->m)*(Cell_ptr->m)*(Cell_ptr->i_O);
            
            break;
        }
        case 7: {//This is for the Full gate immobilization binding scheme with inactive state binding
            
            double vars_final[2]; //Array that will hold the new h and b values
            Gate_Immob(Cell_ptr->h, h_inf, tau_h, Cell_ptr->b, b_inf, tau_b, vars_final);
            
            Cell_ptr->h = vars_final[0];
            Cell_ptr->b = vars_final[1];//Saving the new h and b values
            
            Cell_ptr->I_Na = G_Na*(Cell_ptr->V-E_Na)*(Cell_ptr->m)*(Cell_ptr->m)*(Cell_ptr->m)*(Cell_ptr->h);
            
            break;
        }
        case 8: {//This is for the Full gate immobilization binding scheme with non-inactivated state binding
            
            double h_bar = 1.0 - (Cell_ptr->h) - (Cell_ptr->b); //Defining h_bar (fraction of fast inactivated channels)
            double h_bar_inf = 1.0 - h_inf; //Defining h_bar_inf
            
            double vars_final[2]; //Array that will hold the new h and b values
            Gate_Immob(h_bar, h_bar_inf, tau_h, Cell_ptr->b, b_inf, tau_b, vars_final);
            
            Cell_ptr->h = 1.0 - vars_final[0] - vars_final[1]; //Updating h (1 - h_bar - b)
            Cell_ptr->b = vars_final[1];//Updating b
            
            Cell_ptr->I_Na = G_Na*(Cell_ptr->V-E_Na)*(Cell_ptr->m)*(Cell_ptr->m)*(Cell_ptr->m)*(Cell_ptr->h);
            
            break;
        }
        case 9: { //This is for nonstate specific drug binding.
            Cell_ptr->h = Rush_Larsen(Cell_ptr->h, h_inf, tau_h);
            Cell_ptr->b = Rush_Larsen(Cell_ptr->b, b_inf, tau_b);
            
            Cell_ptr->I_Na = G_Na*(Cell_ptr->V-E_Na)*(Cell_ptr->m)*(Cell_ptr->m)*(Cell_ptr->m)*(Cell_ptr->h)*(1-(Cell_ptr->b));
            
            break;
        }
        default: { //This is for the default, with no drug binding.
            Cell_ptr->h = Rush_Larsen(Cell_ptr->h, h_inf, tau_h);
            
            Cell_ptr->I_Na = G_Na*(Cell_ptr->V-E_Na)*(Cell_ptr->m)*(Cell_ptr->m)*(Cell_ptr->m)*(Cell_ptr->h);
            
            break;
        }
    }
    
}

#endif /* I_Na_Simple_h */
