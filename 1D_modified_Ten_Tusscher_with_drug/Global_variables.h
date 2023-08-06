//
//  Global_variables.h
//  1D_modified_Ten_Tusscher_with_drug_2017_Summer
//
//  Created by Steffen Docken on 7/12/17.
//  Copyright © 2017 Steffen Docken (Lewis Lab). All rights reserved.
//  Published in the Journal of Theoretical Biology under the title "Rate-dependent effects of state-specific sodium channel blockers in cardiac tissue: Insights from idealized models"
//

////7-12-17 This code was originally copied from Ten_Tusscher_model_with_drug_2017_Winter.  Some alterations have been made for 1-D simulations

////This code for the Ten Tusscher model was originally copied from code attached to "Jonathan D. Moreno, Z. Iris Zhu, Pei-Chi Yang, John R. Bankston, Mao-Tsuen Jeng, Chaoyi Kang, Lian- guo Wang, Jason D. Bayer, David J. Christini, Natalia A. Trayanova, Crystal M. Ripplinger, Robert S. Kass, and Colleen E. Clancy. A computational model to predict the effects of class 1 anti-arrhythmic drugs on ventricular rhythms. Science Translational Medicine, 3, 2011."  Code in this file was copied from the Global_variables.h file of the Single_Cell folder of the Moreno source code.  Changes were then made to the code.

#ifndef Global_variables_h
#define Global_variables_h

#include <math.h>
/*********************************************************************************************************
 Universal Constants
 *********************************************************************************************************/

//Ion Valences and Universal Constants
const double R = 8314.472;				// J/mol*K
const double T = 310;					// K
const double F = 96485.3415;            // C/mol
const double Cm = 2.0;					// uF/ cm^2
const double CAP = 0.185;				//Cellular capacitance  From Ten Tusscher code
const double rho = 162.;					// ohm*cm
const double z_Na = 1.;
const double z_Ca = 2.;
const double z_K = 1.;
const double dx= 0.01;					//in cm

//const int nCell = 110;                  //setting number of cells in 1.1 cm cable.
const int nCell_cable_VW = 510;                //setting number of cells in 5.1 cm cable.
const int nCell = nCell_cable_VW;

//Cell Geometry
const long double pi = 3.141592;
const double S_cg = 0.2;				//Surface to volume ratio (um^-1)
const double Diffusion = 0.00154;       //Diffusion coefficient for V from ten Tusscher et al. 2004

//Intracellular volumes
const double V_cyto=0.016404;
const double V_sr=0.001094;
const double V_ss=0.00005468; //Not sure what these units are, since they don't match with what is in the paper, but do match what is in the code.

//Extraceullular concentrations (mM)
const double Na_out = 140;
const double Ca_out = 2.0;
const double K_out = 5.4;

//Drug parameters
double Diffusion_drug = 1.0e2; //This is a default value, but most programs change this value
const double kd0_drug = 10.0e-6; //Kd at 0mV (M)
const double dd_drug = 0.7;
const double Drug = 20.0e-6; //Concentration of drug (M^-1)

//Beating parameters
const double I_duration_S2 = 2.0; //length of extra stimulus
const double stimulus = -400; //Strength of extra stimulus
int flag_S2 = 0; // this is the flag used to create an early stimulus when needed for VW simulations.
const double dt = 0.001;
const double waitTime = 10000; //10 seconds  //30*60000; //30 minutes for I.C.				//If sim_type =0, waitTime ==30000; if sim_type =1, waitTime==0.01;

const int Na_model = 1; //Determines whether TT (0) or the Simple model (1) is used for I_Na
//const int Bind_Scheme = 0; // Determines the drug binding scheme used. (0) No drug, (1) HH Guarded Receptor Inactive State, (2) HH Guarded Receptor Non-inactivated State, (3) HH Gate Immobilization Inactive State, (4) HH Gate Immobilization Non-inactivated State, (5) Full Guarded Receptor Inactive State, (6) Full Guarded Receptor Non-inactivated State, (7) Full Gate Immobilization Inactive State, (8) Full Gate Immobilization Non-inactivated State, (9) Nonstate specific binding0
int Bind_Scheme = 0; //For programs that allow Bind_Scheme to change


/*********************************************************************************************************
 Structures
 *********************************************************************************************************/
struct Cell_param{ //Moreno Code uses typedef to define this.
    double V, dV, V_new;
    double Na_in, K_in, Ca_in, Ca_sr, Ca_ss, Ca_in_buffer, Ca_ss_buffer, Ca_sr_buffer;
    double I_Na, I_Na_L, I_Ca_L, I_Kr, I_Ks, I_K1, I_Kp, I_to, I_Na_Ca, I_Na_K, I_p_Ca, I_Ca_b, I_Na_b, I_stim;
    double I_Na_ion_total, I_Ca_ion_total, I_K_ion_total, I_total, I_axial;
    double I_tr, I_leak, I_up, I_rel;
    int Cell_type;
    double m, h, j, b, i_O, d, f, f2, f_Ca, r, s, xr1, xr2, xs, OO, R_bar; //b was added (which is the fraction of channels bound to drug), and i_O (which is the fraction of channels channels in the non-inactivated, non-drug bound state in the Full Guarded Receptor model)
    
    //double mL, hL; //These variables are for a heart failure model not used here
    
    double peak_slope, t_min, V_min, t_thr, V_thr, t_thr_old, t_max, V_max, I_Na_peak, t_I_Na_peak, I_Na_L_peak, t_I_Na_L_peak;
    double t_APD90, V_90, dV_old;
    int APD90_flag, cycle_num, AP_flag;
    double APD_90, DI;
    double t_APD90_old;
    double CV_1cell, CV_2cell, CV_10cell, T; //CV_1cell and CV_2cell calculate CV based on time it takes AP to travel 1 or 2 cells (to see if there is a difference in the estimation)
    
    double I_total_old, dI_total, dI_total_old;
};

struct simState {
    double t, tTemp;
    int beat;
    Cell_param cellData[nCell];
    
};// this is a structure type for an array of Cell_params.

#endif /* Global_variables_h */
