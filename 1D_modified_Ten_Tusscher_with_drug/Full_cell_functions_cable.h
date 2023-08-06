//
//  Full_cell_functions_cable.h
//  1D_modified_Ten_Tusscher_with_drug_2017_Summer
//
//  Created by Steffen Docken on 8/18/17.
//  Copyright © 2017 Steffen Docken (Lewis Lab). All rights reserved.
//  Published in the Journal of Theoretical Biology under the title "Rate-dependent effects of state-specific sodium channel blockers in cardiac tissue: Insights from idealized models"
//

#ifndef Full_cell_functions_cable_h
#define Full_cell_functions_cable_h

#include <math.h>
#include "Global_variables.h"

//This code is edited from Functions.h in Jonathan Moreno's 100 Cell cable code in the ScienceTM Code.

void Calculate_Reset (Cell_param *Cell_ptr, double t){//resetting variables that store cell characteristics (like peak upstroke velocity or APD90) and variables used to calculate these characteristics.
    Cell_ptr->t_APD90_old = Cell_ptr->t_APD90;
    Cell_ptr->t_APD90 = 0.0;
    Cell_ptr->peak_slope = 0.0;
    Cell_ptr->CV_1cell = 0.0;
    Cell_ptr->CV_2cell = 0.0;
    Cell_ptr->CV_10cell = 0.0;
    Cell_ptr->V_min = Cell_ptr->V;
    Cell_ptr->V_max = -40;
    Cell_ptr->t_thr_old = Cell_ptr->t_thr;
    
    Cell_ptr->V_thr = 0.0;
    Cell_ptr->DI = 0.0;
    Cell_ptr->T = 0.0;
    Cell_ptr->AP_flag = 0; //resetting flag for the fact that V went above -40mV
    
}

void Calculate_Update (Cell_param *Cell_ptr, double t, int n){
    Cell_ptr[n].dV_old = Cell_ptr[n].dV;
    
    if (n==0){
        Cell_ptr[n].I_axial = (Cell_ptr[n+1].V - Cell_ptr[n].V)*Diffusion/(dx*dx); //This takes into account the zero flux boundary condition.
    }
    else if (n>=1 && n<(nCell-1)){Cell_ptr[n].I_axial = (Cell_ptr[n+1].V - 2*Cell_ptr[n].V + Cell_ptr[n-1].V)*Diffusion/(dx*dx);} //This is the centered difference approximation of the second derivative of V multiplied by the Diffusion constant.
    else if (n==(nCell-1)){
        Cell_ptr[n].I_axial = (- Cell_ptr[n].V + Cell_ptr[n-1].V)*Diffusion/(dx*dx); //This takes into account the zero flux boundary condition.
    } /* Diffusion is the diffusion constant, which is 1/(rho*S*C_m), where C_m (uF/cm^2) is the membrane capacitance, rho (Ohm cm) is the cellular resistivity, and S (cm^-1) is the surface to volume ratio. Therefore, the units of Diffusion are cm^2/ms, and the units of I_axial are mV/ms */
    
    Cell_ptr[n].dV = (Cell_ptr[n].I_axial - Cell_ptr[n].I_total)*dt; /* The partial differential equation for V is dV/dt = - I_ion/C_m + 1/(rho*S*C_m)*d^2V/dx^2, where C_m is the membrane capacitance, rho is the cellular resistivity, and S is the surface to volume ratio.  In the ten Tusscher et al. model, dividing by C_m has already been factored into the maximal conductance parameters for the various currents.  Therefore, the units for I_total are (nS/pF)mV or equivalently mV/ms. I_axial also has units of mV/ms, because 1/(rho*S*C_m) (written as the Diffusion constant here) has units of cm^2/ms, as explained above.  For more details on the derivation of these formulas, look at the book "Mathematical Physiology 1: Cellular Physiology" Second Edition by Keener and Sneyd, the article "Instabilities of a Propagating Pulse in a Ring of Excitable Media" 1993 by Courtemanche, Glass, and Keener, and my notes from 9-22-17. */
    Cell_ptr[n].V_new = Cell_ptr[n].V + Cell_ptr[n].dV;
}

void Calculate_I_total(Cell_param *Cell_ptr, double t, int n){
    Cell_ptr->I_total_old = Cell_ptr->I_total;
    Cell_ptr->dI_total_old = Cell_ptr->dI_total;
    
    if ((flag_S2 == 1) && (t <= I_duration_S2 - dt/2.0) && (n == 259)) {Cell_ptr->I_stim=stimulus;} //only middle cell is stimulated
    else {Cell_ptr->I_stim=0;} //Setting initial cells to suprathreshold, so don't need this.
    
    Cell_ptr->I_Na_ion_total = Cell_ptr->I_Na + Cell_ptr->I_Na_b + 3*Cell_ptr->I_Na_K + 3*Cell_ptr->I_Na_Ca;
    /*Cell_ptr->I_Na_ion_total = Cell_ptr->I_Na + Cell_ptr->I_Na_L + Cell_ptr->I_Na_b + 3*Cell_ptr->I_Na_K + 3*Cell_ptr->I_Na_Ca;*/ //I_Na_L used if heart failure model included.
    Cell_ptr->I_K_ion_total = Cell_ptr->I_Kr + Cell_ptr->I_Ks + Cell_ptr->I_K1 + Cell_ptr->I_Kp - 2*Cell_ptr->I_Na_K + Cell_ptr->I_to;
    Cell_ptr->I_Ca_ion_total = Cell_ptr->I_Ca_L + Cell_ptr->I_p_Ca + Cell_ptr->I_Ca_b - 2*Cell_ptr->I_Na_Ca;
    Cell_ptr->I_total = Cell_ptr->I_stim + Cell_ptr->I_Na_ion_total + Cell_ptr->I_K_ion_total + Cell_ptr->I_Ca_ion_total;
    
    Cell_ptr->dI_total = (Cell_ptr->I_total - Cell_ptr->I_total_old)/dt;
}

void Calculate_Points(Cell_param *Cell_ptr, double t, int n){
    
    if ((Cell_ptr[n].AP_flag == 0) && (Cell_ptr[n].V > -40.0)) {
        Cell_ptr[n].AP_flag = 1; //updating flag that V went above -40mV
    }
    
    //Threshold voltage point catcher needed for t_min
    if ((Cell_ptr[n].AP_flag == 1) && (Cell_ptr[n].dV/dt) > Cell_ptr[n].peak_slope) {
        
        Cell_ptr[n].peak_slope = Cell_ptr[n].dV/dt; //recording peak slope
        Cell_ptr[n].V_thr = Cell_ptr[n].V; //recording V of peak slope
        Cell_ptr[n].t_thr = t; //recording time of peak slope
        
        if (n == 1){
            Cell_ptr[n].CV_1cell = dx/(Cell_ptr[n].t_thr - Cell_ptr[n-1].t_thr); //recording CV based on time it takes AP to travel 1 cell (not recorded for the first cell, because there is no cell beforehand to compare to).
        }else if ((n > 1) && (n<10)){
            Cell_ptr[n].CV_1cell = dx/(Cell_ptr[n].t_thr - Cell_ptr[n-1].t_thr); //recording CV based on time it takes AP to travel 1 cell (not recorded for the first cell, because there is no cell beforehand to compare to).
            
            Cell_ptr[n-1].CV_2cell = 2.0*dx/(Cell_ptr[n].t_thr - Cell_ptr[n-2].t_thr); //recording CV based on time it takes AP to travel 2 cells (centered difference), but this calculates CV at the previous cell, so it is saved there (not recorded for the first (or last) cell, because there is no cell beforehand (or after) to compare to).
        }else if (n >=10){
            Cell_ptr[n].CV_1cell = dx/(Cell_ptr[n].t_thr - Cell_ptr[n-1].t_thr); //recording CV based on time it takes AP to travel 1 cell (not recorded for the first cell, because there is no cell beforehand to compare to).
            
            Cell_ptr[n-1].CV_2cell = 2.0*dx/(Cell_ptr[n].t_thr - Cell_ptr[n-2].t_thr); //recording CV based on time it takes AP to travel 2 cells (centered difference), but this calculates CV at the previous cell, so it is saved there (not recorded for the first (or last) cell, because there is no cell beforehand (or after) to compare to).
            
            Cell_ptr[n-5].CV_10cell = 10.0*dx/(Cell_ptr[n].t_thr - Cell_ptr[n-10].t_thr);//recording CV based on time it takes AP to travel 10 cells (centered difference), but this calculates CV at the 5th previous cell, so it is saved there (not recorded for the first (or last) cell, because there is no cell beforehand (or after) to compare to).
        }
        
        //Diastolic Interval
        Cell_ptr[n].DI = Cell_ptr[n].t_thr - Cell_ptr[n].t_APD90;
        //Period
        Cell_ptr[n].T = Cell_ptr[n].t_thr - Cell_ptr[n].t_thr_old;
    }
    
    if ((Cell_ptr[n].dV_old < 0.0)&&(Cell_ptr[n].dV > 0.0)&&(Cell_ptr[n].V < Cell_ptr[n].V_min)) {
        
        Cell_ptr[n].V_min = Cell_ptr[n].V; //recording minimum V during diastolic interval
        Cell_ptr[n].t_min = t; //recording time of minimum V during diastolic interval
        
    }
    
    if ((Cell_ptr[n].AP_flag == 1) && (Cell_ptr[n].dV_old > 0.0)&&(Cell_ptr[n].dV < 0.0)&&(Cell_ptr[n].V > Cell_ptr[n].V_max)) {
        
        Cell_ptr[n].V_max = Cell_ptr[n].V;  //recording maximum V
        Cell_ptr[n].t_max = t; //recording time of max V
        
        Cell_ptr[n].V_90 = Cell_ptr[n].V_min + 0.1*(Cell_ptr[n].V_max - Cell_ptr[n].V_min); //calculating potential at which the cell is 90% repolarized for APD90
        
        Cell_ptr[n].APD90_flag = 1; //stating APD90 can now be calculated.
    }
    
    if ((Cell_ptr[n].V < Cell_ptr[n].V_90)&&(Cell_ptr[n].APD90_flag == 1)) {
        
        Cell_ptr[n].t_APD90 = t; //saving time of APD90
        Cell_ptr[n].APD_90 = t - Cell_ptr[n].t_thr; //Calculating APD90
        
        Cell_ptr[n].APD90_flag = 0; //ensures APD90 is not recalculated
        Cell_ptr[n].cycle_num = Cell_ptr[n].cycle_num + 1; //increases the counter for the period the cell is on by one.
        
    }
    
}


#endif /* Full_cell_functions_cable_h */
