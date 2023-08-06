//
//  load_IC.h
//  Ten_Tusscher_model_with_drug
//
//  Created by Steffen Docken on 12/5/17.
//  Copyright © 2017 Steffen Docken (Lewis Lab). All rights reserved.
//  Published in the Journal of Theoretical Biology under the title "Rate-dependent effects of state-specific sodium channel blockers in cardiac tissue: Insights from idealized models"
//

#ifndef load_IC_h
#define load_IC_h

#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>
//#include <omp.h>
#include <unistd.h>
#include "Global_variables.h"

//This program sets up the initial conditions for the Cell at the beginning of each protocol
void load_IC(Cell_param *Cell_ptr){
    if (sim_type ==1)    {
        //Loading in initial conditions for each cell in the array
        Cell_ptr->V=-86.2;
        //Cell_ptr->V_new=-86.2; //Not used
        Cell_ptr->dV=0;
        Cell_ptr->Na_in = 7.67;
        Cell_ptr->K_in = 138.3;
        Cell_ptr->Ca_in = 0.00007;
        Cell_ptr->Ca_sr = 1.3;
        Cell_ptr->Ca_ss = 0.00007;
        Cell_ptr->Ca_in_buffer = 0;
        Cell_ptr->Ca_ss_buffer = 0;
        Cell_ptr->Ca_sr_buffer = 0;
        Cell_ptr->m = 0.0 ;
        Cell_ptr->h = 0.75;
        if (Na_model == 0) {
            Cell_ptr->j = 0.75;
        }else if (Na_model == 1){
            Cell_ptr->j = 1.0;
        }
        Cell_ptr->b = 0.0;
        Cell_ptr->i_O = 0.75;
        
        //Cell_ptr->mL =  0.00111859;
        //Cell_ptr->hL =  0.339310414;//not used if there is no heart failure included
        
        Cell_ptr->d = 0;
        Cell_ptr->f = 1;
        Cell_ptr->f2 = 1;
        Cell_ptr->f_Ca = 1;
        Cell_ptr->r = 0;
        Cell_ptr->s = 1;
        Cell_ptr->xr1 = 0;
        Cell_ptr->xr2 = 1;
        Cell_ptr->xs = 0;
        Cell_ptr->OO = 0;
        Cell_ptr->R_bar = 1;
        
        Cell_ptr->I_Na = 0;
        Cell_ptr->I_Na_L = 0;
        Cell_ptr->I_Ca_L = 0;
        Cell_ptr->I_Kr = 0;
        Cell_ptr->I_Ks = 0;
        Cell_ptr->I_K1 = 0;
        Cell_ptr->I_Kp = 0;
        Cell_ptr->I_to = 0;
        Cell_ptr->I_Na_Ca = 0;
        Cell_ptr->I_Na_K = 0;
        Cell_ptr->I_p_Ca = 0;
        Cell_ptr->I_Ca_b = 0;
        Cell_ptr->I_Na_b = 0;
        Cell_ptr->I_stim = 0;
        Cell_ptr->I_tr = 0;
        Cell_ptr->I_leak = 0;
        Cell_ptr->I_up = 0;
        Cell_ptr->I_rel = 0;
        
        Cell_ptr->I_Na_ion_total = 0;
        Cell_ptr->I_K_ion_total = 0;
        Cell_ptr->I_Ca_ion_total = 0;
        Cell_ptr->I_total = 0.0;
        Cell_ptr->I_axial=0.0;
        
        Cell_ptr->peak_slope = 0;
        Cell_ptr->V_min = -88.654973;
        Cell_ptr->t_min = 0;
        Cell_ptr->V_thr = -88.654973;
        Cell_ptr->t_thr = 0;
        Cell_ptr->V_max = -88.654973;
        Cell_ptr->t_max = 0;
        Cell_ptr->t_EAD= 0;
        Cell_ptr->V_EAD= -88.654973;
        Cell_ptr->t_EAD2= 0;
        Cell_ptr->V_EAD2= -88.654973;
        Cell_ptr->V_90 = -88.654973;
        Cell_ptr->t_90 = 0;
        Cell_ptr->t_90_old = 0;
        Cell_ptr->dV_old = 0;
        Cell_ptr->flag2 = 0;
        Cell_ptr->flag_EAD=0;
        Cell_ptr->flag_EAD2=0;
        
        
    } //Initial Conditions; sim_type ==1
    
}

#endif /* load_IC_h */
