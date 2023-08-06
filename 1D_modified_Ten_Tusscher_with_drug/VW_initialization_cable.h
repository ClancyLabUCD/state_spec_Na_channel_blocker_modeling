//
//  VW_initialization_cable.h
//  1D_modified_Ten_Tusscher_with_drug_2017_Summer
//
//  Created by Steffen Docken on 9/14/17.
//  Copyright © 2017 Steffen Docken (Lewis Lab). All rights reserved.
//  Published in the Journal of Theoretical Biology under the title "Rate-dependent effects of state-specific sodium channel blockers in cardiac tissue: Insights from idealized models"
//

#ifndef VW_initialization_cable_h
#define VW_initialization_cable_h

//FOR RUNNING THIS CODE, CHANGE nCell TO 510 FOR A 5.1 CM LENGTH OF TISSUE

#define use_omp //This is commented out to exclude using omp.h (not do parallelization)

#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>

#ifdef use_omp
#include <omp.h>  //Includes omp.h if parallelization will be used
#else
#include <ctime>  //Includes ctime for timing program if omp.h is not used
#endif

#include <unistd.h>
#include "Global_variables.h"
#include "Non_Na_currents.h"
#include "Ion_conc_dynamics.h"
#include "Na_currents.h"
#include "Full_cell_functions_cable.h"

void VW_initialization_cable(){
    
#ifdef use_omp
    double start = omp_get_wtime();
#else
    clock_t start_s = clock();
#endif   //initial time for either parallelized or non-parallelized code.
    
    simState Cell_struct;
    simState *C_ptr = &Cell_struct;
    int n;
    
    double BCL;
    
    double kon = 1.0e1; //in M^-1ms^-1
    
    std::string kon_val = "1e1";
    
    std::string Model_type;
    switch (Na_model) {
        case 1:{
            Model_type = "Simple";
            break;
        }
            
        default:{
            Model_type = "Ten_Tusscher";
            break;
        }
    }
    
    std::string dt_val, dx_val, BCL_val;
    
    if (dt == 0.001) {
        dt_val = "1eneg3";
    }else if (dt == 0.0005){
        dt_val = "5eneg4";
    }//setting value of dt to be used in saving the results files
    
    if (dx == 0.01) {
        dx_val = "1eneg2";
    }else if (dx == 0.005){
        dx_val = "5eneg3";
    }//setting value of dx to be used in saving the results files
    
    double BCL_vec[3] = {1.0e3, 6.75e2, 4.0e2}; //list of BCLs to be used
    
    std::ofstream result_file_time, result_file_480, result_file_500, result_file;  // These files will contain the total time of the simulation, cell characteristics after 480 stimuli, cell characteristics after 500 stimuli, and the variables for each cell following the full 500 stimuli.
    FILE *f_final_vars;  //pointer to file for saving initial conditions for next simulation
    
    //Initilizations
    std::string bindscheme;
    
    for (int jj = 0; jj < 3; jj++) {
        
        BCL = BCL_vec[jj]; //setting the BCL for 500 pacings
        
        if (jj == 0){
            BCL_val = "1000";
        }else if(jj == 1){
            BCL_val = "675";
        }else if(jj == 2){
            BCL_val = "400";
        }//setting the BCL to be used in saving the results files
        
        for(Bind_Scheme = 0; Bind_Scheme < 5; Bind_Scheme++){
            C_ptr->t = 0.0;  //initializing t
            
            
            switch (Bind_Scheme) {
                case 1:{
                    bindscheme = "HH_guarded_receptor_inactive";
                    break;
                }
                case 2:{
                    bindscheme = "HH_guarded_receptor_noninactive";
                    break;
                }
                case 3:{
                    bindscheme = "HH_gate_immobilization_inactive";
                    break;
                }
                case 4:{
                    bindscheme = "HH_gate_immobilization_noninactive";
                    break;
                }
                case 5:{
                    bindscheme = "Full_guarded_receptor_inactive";
                    break;
                }
                case 6:{
                    bindscheme = "Full_guarded_receptor_noninactive";
                    break;
                }
                case 7:{
                    bindscheme = "Full_gate_immobilization_inactive";
                    break;
                }
                case 8:{
                    bindscheme = "Full_gate_immobilization_noninactive";
                    break;
                }
                case 9:{
                    bindscheme = "nonstate_specific";
                    break;
                }
                    
                default:{
                    bindscheme = "nodrug";
                    break;
                }
            }//setting the binding scheme to be used in saving the results files
            
            if (Bind_Scheme == 0) {
                
                f_final_vars = fopen((Model_type+'_'+bindscheme+'_'+"VW_initialization_BCL_"+BCL_val+"ms_dt" + dt_val +"_dx"+dx_val+"_state_vars.init_cond").c_str(), "w");
                
                result_file_480.open((Model_type+'_'+bindscheme+'_'+"VW_initialization_BCL_"+BCL_val+"ms_dt" + dt_val +"_dx"+dx_val+"_params_480.txt").c_str());
                
                result_file_500.open((Model_type+'_'+bindscheme+'_'+"VW_initialization_BCL_"+BCL_val+"ms_dt" + dt_val +"_dx"+dx_val+"_params_500.txt").c_str());
                
                result_file.open((Model_type+'_'+bindscheme+'_'+"VW_initialization_BCL_"+BCL_val+"ms_dt" + dt_val +"_dx"+dx_val+"_params_cell260.txt").c_str());
                
            } else {
                
                f_final_vars = fopen((Model_type+'_'+bindscheme+'_'+"VW_initialization_BCL_"+BCL_val+"ms_dt" + dt_val +"_dx"+dx_val+"_kon_" + kon_val+ "_kd10e-6_Drug20e-6_state_vars.init_cond").c_str(), "w");
                
                result_file_480.open((Model_type+'_'+bindscheme+'_'+"VW_initialization_BCL_"+BCL_val+"ms_dt" + dt_val +"_dx"+dx_val+"_kon_" + kon_val+ "_kd10e-6_Drug20e-6_params_480.txt").c_str());
                
                result_file_500.open((Model_type+'_'+bindscheme+'_'+"VW_initialization_BCL_"+BCL_val+"ms_dt" + dt_val +"_dx"+dx_val+"_kon_" + kon_val+ "_kd10e-6_Drug20e-6_params_500.txt").c_str());
                
                result_file.open((Model_type+'_'+bindscheme+'_'+"VW_initialization_BCL_"+BCL_val+"ms_dt" + dt_val +"_dx"+dx_val+"_kon_" + kon_val+ "_kd10e-6_Drug20e-6_params_cell260.txt").c_str());
            }
            
            //Loading in initial conditions for each cell in the array
            for (n=0; n<nCell; n+=1) {
                
                C_ptr->cellData[n].Cell_type = 3; // setting cell type to Epi
                C_ptr->cellData[n].V=-86.2;
                C_ptr->cellData[n].V_new=-86.2;
                C_ptr->cellData[n].dV=0;
                C_ptr->cellData[n].Na_in = 7.67;
                C_ptr->cellData[n].K_in = 138.3;
                C_ptr->cellData[n].Ca_in = 0.00007;
                C_ptr->cellData[n].Ca_sr = 1.3;
                C_ptr->cellData[n].Ca_ss = 0.00007;
                C_ptr->cellData[n].Ca_in_buffer = 0;
                C_ptr->cellData[n].Ca_ss_buffer = 0;
                C_ptr->cellData[n].Ca_sr_buffer = 0;
                C_ptr->cellData[n].m = 0.0;
                C_ptr->cellData[n].h = 0.75;
                C_ptr->cellData[n].j = 1.0;
                C_ptr->cellData[n].b = 0.0;
                
                /*C_ptr->cellData[n].mL =  0.00111859;
                 C_ptr->cellData[n].hL =  0.339310414;*/
                
                C_ptr->cellData[n].d = 0 ;
                C_ptr->cellData[n].f = 1;
                C_ptr->cellData[n].f2 = 1;
                C_ptr->cellData[n].f_Ca = 1;
                C_ptr->cellData[n].r = 0;
                C_ptr->cellData[n].s = 1;
                C_ptr->cellData[n].xr1 = 0;
                C_ptr->cellData[n].xr2 = 1;
                C_ptr->cellData[n].xs = 0;
                C_ptr->cellData[n].OO = 0;
                C_ptr->cellData[n].R_bar = 1;
                
                C_ptr->cellData[n].I_Na = 0;
                C_ptr->cellData[n].I_Na_L = 0;
                C_ptr->cellData[n].I_Ca_L = 0;
                C_ptr->cellData[n].I_Kr = 0;
                C_ptr->cellData[n].I_Ks = 0;
                C_ptr->cellData[n].I_K1 = 0;
                C_ptr->cellData[n].I_Kp = 0;
                C_ptr->cellData[n].I_to = 0;
                C_ptr->cellData[n].I_Na_Ca = 0;
                C_ptr->cellData[n].I_Na_K = 0;
                C_ptr->cellData[n].I_p_Ca = 0;
                C_ptr->cellData[n].I_Ca_b = 0;
                C_ptr->cellData[n].I_Na_b = 0;
                C_ptr->cellData[n].I_stim = 0;
                C_ptr->cellData[n].I_tr = 0;
                C_ptr->cellData[n].I_leak = 0;
                C_ptr->cellData[n].I_up = 0;
                C_ptr->cellData[n].I_rel = 0;
                
                C_ptr->cellData[n].I_Na_ion_total = 0;
                C_ptr->cellData[n].I_K_ion_total = 0;
                C_ptr->cellData[n].I_Ca_ion_total = 0;
                C_ptr->cellData[n].I_total = 0.0;
                C_ptr->cellData[n].I_axial=0.0;
                
                C_ptr->cellData[n].peak_slope = 0;
                C_ptr->cellData[n].V_min = -88.654973;
                C_ptr->cellData[n].t_min = 0;
                C_ptr->cellData[n].V_thr = -88.654973;
                C_ptr->cellData[n].t_thr = 0;
                C_ptr->cellData[n].V_max = -88.654973;
                C_ptr->cellData[n].t_max = 0;
                C_ptr->cellData[n].V_90 = -88.654973;
                C_ptr->cellData[n].t_APD90 = 0;
                C_ptr->cellData[n].t_APD90_old = 0;
                C_ptr->cellData[n].dV_old = 0;
                C_ptr->cellData[n].APD90_flag = 0;
                
                C_ptr->cellData[n].cycle_num = 0;
                C_ptr->cellData[n].AP_flag = 0;
                
            }	// End of for n = 0 to nCell, loading in of initial conditions
            
            for (int ii = 1; ii <= 500; ii++) { //pacing 500 times
                C_ptr->beat = ii;
                
                
                //Resets cell specific parameters every beat such as APD, DI, V90 etc.
                for (n = 0; n<nCell; n=n+1) {
                    Calculate_Reset (&(C_ptr->cellData[n]), C_ptr->t);
                    
                    if ((double) n < 0.1/dx)  {
                        C_ptr->cellData[n].V = 0.0; //establishing "stimulus" to start next propagated wave
                    }
                }
                
                for (C_ptr->tTemp = 0.0;C_ptr->tTemp < BCL; (C_ptr->tTemp) = (C_ptr->tTemp)+dt){
                    
#ifdef use_omp
#pragma omp parallel for
#endif
                    for (n=0; n<nCell; n=n+1){
                        
                        Calculate_Points (C_ptr->cellData, C_ptr->t, n);  //Calculates t_min, V_min, t_max, V_max etc. etc.;
                        
                        //Updating current calculations for each cell
                        Calculate_I_Na(&(C_ptr->cellData[n]),C_ptr->tTemp , kon);
                        //Calculate_I_Na_L(&(C_ptr->cellData[n]), t );
                        Calculate_I_Ca_L(&(C_ptr->cellData[n]),C_ptr->tTemp );
                        Calculate_I_Kr(&(C_ptr->cellData[n]),C_ptr->tTemp );
                        Calculate_I_Ks(&(C_ptr->cellData[n]),C_ptr->tTemp );
                        Calculate_I_K1(&(C_ptr->cellData[n]),C_ptr->tTemp );
                        Calculate_I_Kp(&(C_ptr->cellData[n]),C_ptr->tTemp );
                        Calculate_I_to(&(C_ptr->cellData[n]),C_ptr->tTemp );
                        Calculate_I_Na_Ca(&(C_ptr->cellData[n]),C_ptr->tTemp );
                        Calculate_I_Na_K(&(C_ptr->cellData[n]),C_ptr->tTemp);
                        Calculate_I_p_Ca(&(C_ptr->cellData[n]),C_ptr->tTemp );
                        Calculate_I_Ca_b(&(C_ptr->cellData[n]),C_ptr->tTemp );
                        Calculate_I_Na_b(&(C_ptr->cellData[n]),C_ptr->tTemp );
                        Calculate_I_total(&(C_ptr->cellData[n]),C_ptr->tTemp, n);
                        
                        //Updating ionic concentrations
                        Calculate_Na_in(&(C_ptr->cellData[n]),C_ptr->tTemp );
                        Calculate_K_in(&(C_ptr->cellData[n]),C_ptr->tTemp );
                        Calculate_I_tr(&(C_ptr->cellData[n]),C_ptr->tTemp);
                        Calculate_I_leak(&(C_ptr->cellData[n]),C_ptr->tTemp);
                        Calculate_I_up(&(C_ptr->cellData[n]),C_ptr->tTemp);
                        Calculate_I_rel(&(C_ptr->cellData[n]),C_ptr->tTemp);
                        Calculate_Ca_sr(&(C_ptr->cellData[n]),C_ptr->tTemp );
                        Calculate_Ca_ss(&(C_ptr->cellData[n]),C_ptr->tTemp );
                        Calculate_Ca_in(&(C_ptr->cellData[n]),C_ptr->tTemp );
                        
                        // Update calculates axial currents and updates the voltage
                        Calculate_Update (C_ptr->cellData,C_ptr->tTemp, n);
                        
                        
                    } // End of for n
                    
                    //Update the new voltages for each cell AFTER the parallelized "for n" loop.
                    
                    for (n=0; n<nCell; n=n+1){C_ptr->cellData[n].V = C_ptr->cellData[n].V_new;}
                    
                    C_ptr->t = C_ptr->t+dt; //update time
                    
                    if (fmod(C_ptr->t,1000)<dt) {std::cout << "t=" << C_ptr->t << ", Binding Scheme is " << Bind_Scheme << ", drug is " << Drug <<  std::endl;} //outputs progress
                    
                }
                
                /* Results section */
                
                if (ii == 480) {
                    for (n = 0; n<nCell; n++) {
                        result_file_480 << n << ", " << C_ptr->cellData[n].CV_1cell << ", " << C_ptr->cellData[n].CV_2cell << ", " <<  C_ptr->cellData[n].CV_10cell << ", " << C_ptr->cellData[n].peak_slope << ", " << C_ptr->cellData[n].APD_90 << ", " << C_ptr->cellData[n].V_thr << ", " << C_ptr->cellData[n].cycle_num << std::endl; //This records the cell number, CV_1cell, CV_2cell, CV_10cell, peak upstroke velocity, and APD90 for each cell at the end of the 480th stimulus so that I can check that this is close to steady state
                    }
                    
                }else if (ii == 499){
                    result_file << ii << ", " << BCL << ", " << C_ptr->cellData[259].CV_1cell << ", " << C_ptr->cellData[259].CV_2cell << ", " << C_ptr->cellData[259].CV_10cell << ", " << C_ptr->cellData[259].peak_slope << ", " << C_ptr->cellData[259].APD_90 << ", " << C_ptr->cellData[259].V_thr << ", " << C_ptr->cellData[259].cycle_num  << std::endl; //This records the stimulus number, BCL, CV_1cell, CV_2cell, peak upstroke velocity, and APD90
                    
                }else if (ii == 500){
                    
                    result_file << ii << ", " << BCL << ", " << C_ptr->cellData[259].CV_1cell << ", " << C_ptr->cellData[259].CV_2cell << ", " << C_ptr->cellData[259].CV_10cell << ", " << C_ptr->cellData[259].peak_slope << ", " << C_ptr->cellData[259].APD_90 << ", " << C_ptr->cellData[259].V_thr << ", " << C_ptr->cellData[259].cycle_num  << std::endl; //This records the stimulus number, BCL, CV_1cell, CV_2cell, peak upstroke velocity, and APD90
                    
                    for (n = 0; n<nCell; n++) {
                        
                        result_file_500 << n << ", " << C_ptr->cellData[n].CV_1cell << ", " << C_ptr->cellData[n].CV_2cell << ", " <<  C_ptr->cellData[n].CV_10cell << ", " << C_ptr->cellData[n].peak_slope << ", " << C_ptr->cellData[n].APD_90 << ", " << C_ptr->cellData[n].V_thr << ", " << C_ptr->cellData[n].cycle_num << std::endl; //This records the cell number, CV_1cell, CV_2cell, peak upstroke velocity, and APD90 for each cell at the end of the 500th stimulus so that I can check that this is close to steady state
                    }
                }
                
                
#ifdef use_omp
                double stop_s1 = omp_get_wtime();
                std::cout << "time: " << stop_s1- start << " seconds" << std::endl;
#else
                clock_t stop_s1 = clock();
                std::cout << "time: " << ((double)stop_s1 - (double)start_s)/double(CLOCKS_PER_SEC) << " seconds" << std::endl;
#endif
                
                
            }
            
            
            fwrite(C_ptr, sizeof(simState), 1, f_final_vars);
            fclose(f_final_vars);
            
            result_file.close();
            result_file_500.close();
            result_file_480.close();//closing these files
        }
        
    }
    
    result_file_time.open(("VW_initialization_dt" + dt_val +"_dx"+dx_val+"_timing.txt").c_str());
    
#ifdef use_omp
    double end = omp_get_wtime();
    result_file_time << "time: " << end - start << " seconds" << std::endl;
#else
    clock_t stop_s = clock();
    result_file_time << "time: " << ((double)stop_s - (double)start_s)/double(CLOCKS_PER_SEC) << " seconds" << std::endl;
#endif
    
    result_file_time.close();
    
    
}

#endif /* VW_initialization_cable_h */
