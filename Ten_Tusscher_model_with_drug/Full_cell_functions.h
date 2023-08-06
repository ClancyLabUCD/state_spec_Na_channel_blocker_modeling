//
//  Full_cell_functions.h
//  Ten_Tusscher_model_with_drug
//
//  Created by Steffen Docken on 12/5/17.
//  Copyright © 2017 Steffen Docken (Lewis Lab). All rights reserved.
//  Published in the Journal of Theoretical Biology under the title "Rate-dependent effects of state-specific sodium channel blockers in cardiac tissue: Insights from idealized models"
//

////This code for the Ten Tusscher model was originally copied from code attached to "Jonathan D. Moreno, Z. Iris Zhu, Pei-Chi Yang, John R. Bankston, Mao-Tsuen Jeng, Chaoyi Kang, Lian- guo Wang, Jason D. Bayer, David J. Christini, Natalia A. Trayanova, Crystal M. Ripplinger, Robert S. Kass, and Colleen E. Clancy. A computational model to predict the effects of class 1 anti-arrhythmic drugs on ventricular rhythms. Science Translational Medicine, 3, 2011."  Code in this file was copied from the Functions.h file of the Single_Cell folder of the Moreno source code.  Changes were then made to the code.
////The Functions.h file was split into Non_Na_currents.h, Ion_conc_dynamics.h, and Full_cell_functions.h.

#ifndef Full_cell_functions_h
#define Full_cell_functions_h

#include <math.h>
#include "Global_variables.h"

void Calculate_Reset (Cell_param *Cell_ptr, double t, double tTemp){
    Cell_ptr->I_Na_peak = 0;
    Cell_ptr->I_Na_L_peak = 0;
    Cell_ptr->t_90_old = Cell_ptr->t_90;
    Cell_ptr->flag2 = 0;
    Cell_ptr->flag_EAD=0;
    Cell_ptr->flag_EAD2=0;
    Cell_ptr->peak_slope = 0;
    Cell_ptr->V_min = -88.654973;
    Cell_ptr->t_min = 0;
    Cell_ptr->V_thr = -88.654973;
    Cell_ptr->t_thr = 0;
    Cell_ptr->V_max = -88.654973;
    Cell_ptr->t_max = t;
    Cell_ptr->t_EAD= t;
    Cell_ptr->t_EAD2=t;
    Cell_ptr->t_90=t;
    Cell_ptr->V_EAD= -88.654973;
    Cell_ptr->V_EAD2=-88.654973;
    Cell_ptr->V_90 = -88.654973;
    Cell_ptr->dV_old = 0;
    
    Cell_ptr->flagI_total = 0;//resetting variables that store cell characteristics (like peak upstroke velocity or APD90) and variables used to calculate these characteristics.
    
}

void Calculate_Update (Cell_param *Cell_ptr, double t, double tTemp){
    Cell_ptr->dV_old = Cell_ptr->dV;
    Cell_ptr->dV = -1*(Cell_ptr->I_total)*dt;
    Cell_ptr->V = Cell_ptr->V + Cell_ptr->dV;
}

void Calculate_I_total(Cell_param *Cell_ptr, double t, double tTemp){
    Cell_ptr->I_total_old = Cell_ptr->I_total;
    Cell_ptr->dI_total_old = Cell_ptr->dI_total;
    
    if ((t>waitTime) && (tTemp <=I_duration)) {Cell_ptr->I_stim=stimulus;}
    else {Cell_ptr->I_stim=0;}
    
    Cell_ptr->I_Na_ion_total = Cell_ptr->I_Na + Cell_ptr->I_Na_b + 3*Cell_ptr->I_Na_K + 3*Cell_ptr->I_Na_Ca;
    Cell_ptr->I_K_ion_total = Cell_ptr->I_Kr + Cell_ptr->I_Ks + Cell_ptr->I_K1 + Cell_ptr->I_Kp - 2*Cell_ptr->I_Na_K + Cell_ptr->I_to;
    Cell_ptr->I_Ca_ion_total = Cell_ptr->I_Ca_L + Cell_ptr->I_p_Ca + Cell_ptr->I_Ca_b - 2*Cell_ptr->I_Na_Ca;
    Cell_ptr->I_total = Cell_ptr->I_stim + Cell_ptr->I_Na_ion_total + Cell_ptr->I_K_ion_total + Cell_ptr->I_Ca_ion_total;
    
    Cell_ptr->dI_total = (Cell_ptr->I_total - Cell_ptr->I_total_old)/dt;
}

void Calculate_Points(Cell_param *Cell_ptr, double t, double tTemp){
    if (fabs(Cell_ptr->I_Na) > fabs(Cell_ptr->I_Na_peak) ){
        Cell_ptr->I_Na_peak = Cell_ptr->I_Na;
        Cell_ptr->t_I_Na_peak = t;}
    
    if (fabs(Cell_ptr->I_Na_L) > fabs(Cell_ptr->I_Na_L_peak) ){
        Cell_ptr->I_Na_L_peak = Cell_ptr->I_Na_L;
        Cell_ptr->t_I_Na_L_peak = t;}
    
    //Minimum voltage catcher
    if ((t>=waitTime) && (tTemp<dt)){
        Cell_ptr->t_min = t;
        Cell_ptr->V_min = Cell_ptr->V;
    }
    
    //Threshold voltage point catcher needed for t_min
    if ((t>=waitTime) && ((Cell_ptr->dV/dt) > Cell_ptr->peak_slope)) {
        Cell_ptr->peak_slope = Cell_ptr->dV/dt;
        Cell_ptr->V_thr = Cell_ptr->V;
        Cell_ptr->t_thr = t;
    }
    
    //Maximum voltage point catcher
    if ((t>=waitTime) && (Cell_ptr->dV_old >0.) && (Cell_ptr->dV < 0.) && (Cell_ptr->V >10.)) {
        Cell_ptr->V_max = Cell_ptr->V; //this is incorrect.  V_max is reassigned the peak plateau value with this formulation
        Cell_ptr->t_max = t;
    }
    
    //I_total minimum point
    if ((t>=waitTime) && (Cell_ptr->dI_total_old < 0.) && (Cell_ptr->dI_total> 0.) && (tTemp >200.) && (Cell_ptr->flagI_total ==0)){
        Cell_ptr->flagI_total = 1;
        Cell_ptr->t_I_total = t;
        Cell_ptr->I_total_pt = Cell_ptr->I_total;
        Cell_ptr->V_I_total = Cell_ptr->V;
    }
    
    //Local minima
    if ((t>=waitTime) && (Cell_ptr->dV_old < 0.) && (Cell_ptr->dV > 0.) && (Cell_ptr->V > -40.) && (tTemp > 200.) && (Cell_ptr->flag_EAD==0)) {
        Cell_ptr->flag_EAD=1;
        Cell_ptr->V_EAD = Cell_ptr->V;
        Cell_ptr->t_EAD = t;
    }
    
    //Local maxima
    if ((t > Cell_ptr->t_EAD) && (Cell_ptr->dV_old > 0.) && (Cell_ptr->dV < 0.) && (tTemp > 200.) && (Cell_ptr->flag_EAD2==0)) {
        Cell_ptr->flag_EAD2=1;
        Cell_ptr->V_EAD2 = Cell_ptr->V;
        Cell_ptr->t_EAD2 = t;
    }
    
    //V_90 voltage point catcher, cycle number and APD_90 calculation
    Cell_ptr->V_90 = Cell_ptr->V_max - .90*(Cell_ptr->V_max - Cell_ptr->V_min);
    
    //Diastolic Interval
    Cell_ptr->DI = Cell_ptr->t_thr - Cell_ptr->t_90_old;
    
    //Height (or length) of the EAD
    Cell_ptr->L_EAD=Cell_ptr->V_EAD2 - Cell_ptr->V_EAD;
    
    if ((t>=waitTime) && (fabs(Cell_ptr->V - Cell_ptr->V_90) < (0.015)) && (Cell_ptr->dV < 0.) && (tTemp > 30.) && (t > Cell_ptr->t_max) && (Cell_ptr->flag2 ==0)) {
        Cell_ptr->flag2 = 1;
        Cell_ptr->t_90 = t;
        Cell_ptr->APD_90 = Cell_ptr->t_90 - Cell_ptr->t_thr;
    }
}

#endif /* Full_cell_functions_h */
