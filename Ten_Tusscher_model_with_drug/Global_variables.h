//
//  Global_variables.h
//  Ten_Tusscher_model_with_drug
//
//  Created by Steffen Docken on 12/5/17.
//  Copyright © 2017 Steffen Docken (Lewis Lab). All rights reserved.
//  Published in the Journal of Theoretical Biology under the title "Rate-dependent effects of state-specific sodium channel blockers in cardiac tissue: Insights from idealized models"
//

////This code for the Ten Tusscher model was originally copied from code attached to "Jonathan D. Moreno, Z. Iris Zhu, Pei-Chi Yang, John R. Bankston, Mao-Tsuen Jeng, Chaoyi Kang, Lian- guo Wang, Jason D. Bayer, David J. Christini, Natalia A. Trayanova, Crystal M. Ripplinger, Robert S. Kass, and Colleen E. Clancy. A computational model to predict the effects of class 1 anti-arrhythmic drugs on ventricular rhythms. Science Translational Medicine, 3, 2011."  Code in this file was copied from the Global_variables.h file of the Single_Cell folder of the Moreno source code.  Changes were then made to the code.

#ifndef Global_variables_h
#define Global_variables_h

#include <math.h>
/*********************************************************************************************************
 Universal Constants
 *********************************************************************************************************/

//Ion Valences and Universal Constants
const double R = 8314.472;                // J/mol*K
const double T = 310;                    // K
const double F = 96485.3415;            // C/mol
const double Cm = 2.0;                    // uF/ cm^2
const double CAP = 0.185;                //Cellular capacitance  From Ten Tusscher code
const double rho = 162.;                    // ohm*cm
const double z_Na = 1.;
const double z_Ca = 2.;
const double z_K = 1.;
const double dx= 0.01;                    //.01 //Can go from 0.01 - 0.02 cm

//Cell Geometry
const long double pi = 3.141592;
const double S_cg = 0.2;                //Surface to volume ratio (um^-1)
const double mutant = 0.0;

//Intracellular volumes
const double V_cyto=0.016404;
const double V_sr=0.001094;
const double V_ss=0.00005468; //Not sure what these units are, since they don't match with what is in the paper, but do match what is in the code.

//Extraceullular concentrations (mM)
const double Na_out = 140;
const double Ca_out = 2;
const double K_out = 5.4;

//Drug parameters
double Diffusion_drug = 1.0e2; //This is a default value, but most programs change this value
const double kd0_drug = 10.0e-6; //Kd at 0mV (M)
const double dd_drug = 0.7;
const double Drug = 20.0e-6; // 5.0e-6; //Concentration of drug (M^-1)

//Beating parameters
const double I_duration = 0.6;
const double stimulus = -80;
double dt = 0.005;//0.01; //decreasing dt from 0.01 to 0.005 causes about a 1% change in upstroke velocity for both Ten Tusscher Na and the Simple Na model, so this is at convergence.
const double waitTime = 10000; //10 seconds

const int sim_type = 1;        //0 is for loading initial conditions from a long wait time from a file (which is not in this folder), 1 is for initials conditions that must be allowed to settle to steady state

const int Na_model = 1; //Determines whether TT (0) or our Simple model (1) is used for I_Na
//const int Bind_Scheme = 0; // Determines the drug binding scheme used. (0) No drug, (1) HH Guarded Receptor Inactive State, (2) HH Guarded Receptor Non-inactivated State, (3) HH Gate Immobilization Inactive State, (4) HH Gate Immobilization Non-inactivated State, (5) Full Guarded Receptor Inactive State, (6) Full Guarded Receptor Non-inactivated State, (7) Full Gate Immobilization Inactive State, (8) Full Gate Immobilization Non-inactivated State, (9) Nonstate specific binding0
int Bind_Scheme = 0; //For programs that allow Bind_Scheme to change

/*********************************************************************************************************
 Structures
 *********************************************************************************************************/
struct Cell_param{
    double V, dV;
    double Na_in, K_in, Ca_in, Ca_sr, Ca_ss, Ca_in_buffer, Ca_ss_buffer, Ca_sr_buffer;
    double I_Na, I_Na_L, I_Ca_L, I_Kr, I_Ks, I_K1, I_Kp, I_to, I_Na_Ca, I_Na_K, I_p_Ca, I_Ca_b, I_Na_b, I_stim;
    double I_Na_ion_total, I_Ca_ion_total, I_K_ion_total, I_total, I_axial;
    double I_tr, I_leak, I_up, I_rel;
    int Cell_type;
    double m, h, j, b, i_O, d, f, f2, f_Ca, r, s, xr1, xr2, xs, OO, R_bar; //b was added (which is the fraction of channels bound to drug), and i_O (which is the fraction of channels channels in the non-inactivated, non-drug bound state in the Full Guarded Receptor model)
    
    //double mL, hL; //These variables are for a heart failure model not used here
    
    double peak_slope, t_min, V_min, t_thr, V_thr, t_max, V_max, t_EAD, V_EAD, I_Na_peak, t_I_Na_peak, I_Na_L_peak, t_I_Na_L_peak;
    double t_EAD2, V_EAD2, t_90, V_90, dV_old, flag2;
    double APD_90, DI, L_EAD;
    double t_90_old, flag_EAD, flag_EAD2;
    double CV, CV_EAD;
    
    double I_total_old, dI_total, dI_total_old, t_I_total, I_total_pt, V_I_total, flagI_total;
};

#endif /* Global_variables_h */
