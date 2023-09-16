//
//  I_Na_TT.h
//  1D_modified_Ten_Tusscher_with_drug_2017_Summer
//
//  Created by Steffen Docken on 7/11/17.
//  Copyright © 2017 Steffen Docken (Lewis Lab). All rights reserved.
//  Published in the Journal of Theoretical Biology under the title "Rate-dependent effects of state-specific sodium channel blockers in cardiac tissue: Insights from idealized models"
//

////7-11-17 This code was copied from Ten_Tusscher_model_with_drug_2017_Winter

#ifndef I_Na_TT_h
#define I_Na_TT_h

#include <math.h>
#include "Global_variables.h"
#include "Numerical_Methods.h"

//The ten Tusscher Na-channel model
void Calculate_I_Na_TT(Cell_param *Cell_ptr, double t, double kon){//From Ten Tusscher et al. 2004 paper
    
    double a_m, b_m, a_h, b_h, a_j, b_j;
    double tau_m, tau_h, tau_j, m_inf, h_inf, j_inf;
    const double G_Na=14.838; //nS/pF
    
    double E_Na=(R*T/F)*log(Na_out/Cell_ptr->Na_in);
    
    a_m = 1./(1.0 + exp(-60.0-Cell_ptr->V)/5.0);
    b_m = 0.1/(1.0+exp((Cell_ptr->V+35.0)/5.0))+0.1/(1.0+exp((Cell_ptr->V-50.0)/200.0));
    tau_m = a_m*b_m;
    m_inf = 1.0/pow((1.0+exp((-56.86-Cell_ptr->V)/9.03)),2.0);
    
    if(Cell_ptr->V>=-40.0){
        a_h=0.0;
        b_h=0.77/(0.13*(1.0+exp((Cell_ptr->V+10.66)/-11.1)));
        
        a_j=0.0;
        b_j=0.6*exp(0.057*Cell_ptr->V)/(1.0+exp(-0.1*(Cell_ptr->V+32.0)));
    }
    
    else {
        a_h=0.057*exp((80.0+Cell_ptr->V)/-6.8);
        b_h=2.7*exp(0.079*Cell_ptr->V)+3.1E5*exp(0.3485*Cell_ptr->V);
        
        a_j=((-2.5428E4)*exp(0.2444*Cell_ptr->V)-(6.978E-6)*exp(-0.04391*Cell_ptr->V))*(Cell_ptr->V+37.78)/(1+exp(0.311*(Cell_ptr->V+79.23)));
        b_j=0.02424*exp(-0.01052*Cell_ptr->V)/(1.0+exp(-0.1378*(Cell_ptr->V+40.14)));
    }
    
    tau_h = 1.0/(a_h+b_h);
    tau_j = 1.0/(a_j+b_j);
    
    h_inf = 1.0/pow((1.0+exp((Cell_ptr->V+71.55)/7.43)),2.0);
    j_inf = h_inf;
    
    double k_on = kon;
    double k_off = kon*kd0_drug*exp(-dd_drug*(Cell_ptr->V)*F/(R*T));
    
    double tau_b = 1.0/(Drug*k_on + k_off);
    double b_inf = Drug*k_on/(Drug*k_on + k_off);
    
    Cell_ptr->m = Rush_Larsen(Cell_ptr->m, m_inf, tau_m);
    
    Cell_ptr->j = Rush_Larsen(Cell_ptr->j, j_inf, tau_j);
    
    switch (Bind_Scheme) {
        case 1: {//This is for the HH formulation of guarded receptor binding scheme with inactive state binding.
            Cell_ptr->h = Rush_Larsen(Cell_ptr->h, h_inf, tau_h);
            Cell_ptr->b = Rush_Larsen(Cell_ptr->b, b_inf, tau_b/(1.0-Cell_ptr->h));
            
            Cell_ptr->I_Na = G_Na*(Cell_ptr->V-E_Na)*(Cell_ptr->m)*(Cell_ptr->m)*(Cell_ptr->m)*(Cell_ptr->h)*(Cell_ptr->j)*(1-(Cell_ptr->b));
            
            break;
        }
        case 2: {//This is for the HH formulation of guarded receptor binding scheme with non-inactivated state state binding.
            Cell_ptr->h = Rush_Larsen(Cell_ptr->h, h_inf, tau_h);
            Cell_ptr->b = Rush_Larsen(Cell_ptr->b, b_inf, tau_b/(Cell_ptr->h));
            
            Cell_ptr->I_Na = G_Na*(Cell_ptr->V-E_Na)*(Cell_ptr->m)*(Cell_ptr->m)*(Cell_ptr->m)*(Cell_ptr->h)*(Cell_ptr->j)*(1-(Cell_ptr->b));
            
            break;
        }
        case 3: {//This is for the HH formulation of gate immobilization binding scheme with inact state state binding.
            double h_0 = Cell_ptr->h; //value of h from previous time step
            double b_inf_GISI = (1.0-h_0)*b_inf/(1.0 - h_0*b_inf); //steady state of b for HH Formulation of Gate Immobilization model with inactivated state binding (for the current value of h)
            double tau_b_GISI = tau_b/(1.0-h_0*b_inf); //time constant of b for HH Formulation of Gate Immobilization model with inactivated state binding (for the current value of h)
            
            Cell_ptr->h = Rush_Larsen(Cell_ptr->h, h_inf, tau_h);
            Cell_ptr->b = Rush_Larsen(Cell_ptr->b, b_inf_GISI, tau_b_GISI);
            
            Cell_ptr->I_Na = G_Na*(Cell_ptr->V-E_Na)*(Cell_ptr->m)*(Cell_ptr->m)*(Cell_ptr->m)*(Cell_ptr->h)*(1.0-(Cell_ptr->b))*(Cell_ptr->j);
            
            break;
        }
        case 4: {//This is for the HH formulation of gate immobilization binding scheme with non-inact state state binding.
            double h_0 = Cell_ptr->h; //value of h from previous time step
            double b_inf_GISN = h_0*b_inf/(1.0 - (1.0-h_0)*b_inf); //steady state of b for HH Formulation of Gate Immobilization model with non-inactivated state binding (for the current value of h)
            double tau_b_GISN = tau_b/(1.0-(1.0-h_0)*b_inf); //time constant of b for HH Formulation of Gate Immobilization model with non-inactivated state binding (for the current value of h)
            
            Cell_ptr->h = Rush_Larsen(Cell_ptr->h, h_inf, tau_h);
            Cell_ptr->b = Rush_Larsen(Cell_ptr->b, b_inf_GISN, tau_b_GISN);
            
            Cell_ptr->I_Na = G_Na*(Cell_ptr->V-E_Na)*(Cell_ptr->m)*(Cell_ptr->m)*(Cell_ptr->m)*(Cell_ptr->h)*(1.0-(Cell_ptr->b))*(Cell_ptr->j);
            
            break;
        }
        case 5: {//Full Guarded Receptor Inactive State
            double h_0 = Cell_ptr->h; //value of h from previous time step
            double b_0 = Cell_ptr->b; //value of b from previous time step
            double i_O_0 = Cell_ptr->i_O; //value of i_O from previous time step
            
            double b_inf_GRFI = (1.0-h_0)*b_inf + h_0 - i_O_0; //steady state of b for Full Formulation of Guarded Receptor model with inactivated state binding (for the current value of h)
            double i_O_inf = (1.0 - b_0)*h_inf;
            
            Cell_ptr->h = Rush_Larsen(Cell_ptr->h, h_inf, tau_h);
            Cell_ptr->b = Rush_Larsen(Cell_ptr->b, b_inf_GRFI, tau_b);
            Cell_ptr->i_O = Rush_Larsen(Cell_ptr->i_O, i_O_inf, tau_h);
            
            Cell_ptr->I_Na = G_Na*(Cell_ptr->V-E_Na)*(Cell_ptr->m)*(Cell_ptr->m)*(Cell_ptr->m)*(Cell_ptr->i_O)*(Cell_ptr->j);
            
            break;
        }
        case 6: {//Full Guarded Receptor Non-inactive State
            double h_0 = Cell_ptr->h; //value of h from previous time step
            double b_0 = Cell_ptr->b; //value of b from previous time step
            double i_C_0 = 1.0 - Cell_ptr->b - Cell_ptr->i_O; //value of i_O from previous time step
            
            double b_inf_GRFN = (1.0 - (1.0 - b_inf)*h_0 - i_C_0);  //steady state of b for Full Formulation of Guarded Receptor model with non-inactivated state binding (for the current value of h)
            double i_C_inf = (1.0 - h_inf)*(1.0 - b_0);
            
            Cell_ptr->h = Rush_Larsen(Cell_ptr->h, h_inf, tau_h);
            Cell_ptr->b = Rush_Larsen(Cell_ptr->b, b_inf_GRFN, tau_b);
            double i_C_new = Rush_Larsen(i_C_0, i_C_inf, tau_h);
            
            Cell_ptr->i_O = 1.0 - Cell_ptr->b - i_C_new;
            
            Cell_ptr->I_Na = G_Na*(Cell_ptr->V-E_Na)*(Cell_ptr->m)*(Cell_ptr->m)*(Cell_ptr->m)*(Cell_ptr->i_O)*(Cell_ptr->j);
            
            break;
        }
        case 7: {//This is for the Full gate immobilization binding scheme with inactive state binding
            
            double vars_final[2]; //Array that will hold the new h and b values
            Gate_Immob(Cell_ptr->h, h_inf, tau_h, Cell_ptr->b, b_inf, tau_b, vars_final);
            
            Cell_ptr->h = vars_final[0];
            Cell_ptr->b = vars_final[1];//Saving the new h and b values
            
            Cell_ptr->I_Na = G_Na*(Cell_ptr->V-E_Na)*(Cell_ptr->m)*(Cell_ptr->m)*(Cell_ptr->m)*(Cell_ptr->h)*(Cell_ptr->j);
            
            break;
        }
        case 8: {//This is for the Full gate immobilization binding scheme with non-inactivated state binding
            
            double h_bar = 1.0 - (Cell_ptr->h) - (Cell_ptr->b); //Defining h_bar (fraction of fast inactivated channels)
            double h_bar_inf = 1.0 - h_inf; //Defining h_bar_inf
            
            double vars_final[2]; //Array that will hold the new h and b values
            Gate_Immob(h_bar, h_bar_inf, tau_h, Cell_ptr->b, b_inf, tau_b, vars_final);
            
            Cell_ptr->h = 1.0 - vars_final[0] - vars_final[1]; //Updating h (1 - h_bar - b)
            Cell_ptr->b = vars_final[1];//Updating b
            
            Cell_ptr->I_Na = G_Na*(Cell_ptr->V-E_Na)*(Cell_ptr->m)*(Cell_ptr->m)*(Cell_ptr->m)*(Cell_ptr->h)*(Cell_ptr->j);
            
            break;
        }
        case 9: { //This is for nonstate specific drug binding.
            Cell_ptr->h = Rush_Larsen(Cell_ptr->h, h_inf, tau_h);
            Cell_ptr->b = Rush_Larsen(Cell_ptr->b, b_inf, tau_b);
            
            Cell_ptr->I_Na = G_Na*(Cell_ptr->V-E_Na)*(Cell_ptr->m)*(Cell_ptr->m)*(Cell_ptr->m)*(Cell_ptr->h)*(Cell_ptr->j)*(1-(Cell_ptr->b));
            
            break;
        }
        default: { //This is for the default, with no drug binding.
            Cell_ptr->h = Rush_Larsen(Cell_ptr->h, h_inf, tau_h);
            
            Cell_ptr->I_Na = G_Na*(Cell_ptr->V-E_Na)*(Cell_ptr->m)*(Cell_ptr->m)*(Cell_ptr->m)*(Cell_ptr->h)*(Cell_ptr->j);
            
            break;
        }
    }
    
}

#endif /* I_Na_TT_h */
