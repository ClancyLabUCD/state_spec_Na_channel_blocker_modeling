//
//  CV_vs_BCL_cable.h
//  1D_modified_Ten_Tusscher_with_drug_2017_Summer
//
//  Created by Steffen Docken on 8/18/17.
//  Copyright © 2017 Steffen Docken (Lewis Lab). All rights reserved.
//  Published in the Journal of Theoretical Biology under the title "Rate-dependent effects of state-specific sodium channel blockers in cardiac tissue: Insights from idealized models"
//

#ifndef CV_vs_BCL_cable_h
#define CV_vs_BCL_cable_h

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

void CV_vs_BCL_cable(){
    
#ifdef use_omp
    double start = omp_get_wtime();
#else
    clock_t start_s = clock();
#endif   //initial time for either parallelized or non-parallelized code.
    
    Cell_param Cell[nCell];
    int n;
    double t, tTemp;
    
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
    
    std::string dt_val, dx_val; //These strings will be used to write the file names for the various results files.
    
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
    
    //chdir ("/Users/steffendocken/Documents/Research/2017_Summer/Ten_Tusscher_1D_sims_output/sandbox/cable/CV_vs_BCL_9_4_17/"); //for my laptop
    std::ofstream result_file_Vs_1000ms, result_file_Vs_300ms, result_file_time, result_file_480, result_file_500, result_file;  // These files will contain the results.  The first will contain the Vs for all cells every .5 ms for the last 2 stimuli with a BCL of 1000ms, the second will contain the Vs for all cells every .5 ms for the last 2 stimuli with a BCL of 300ms, the third contains the total time it takes to run the simulation, the fourth and fifth contain cell characteristics for all cells after the 480th and 500th stimuli, respectively, to check that steady state has been reached, the final file contains peak upstroke velocity (and other cell characteristics) of the cell .6cm down the cable during the final two stimuli at each BCL.
    
    //Initilizations
    std::string bindscheme;
    
    for(Bind_Scheme = 0; Bind_Scheme < 5; Bind_Scheme++){
        t = 0.0;  //initializing t
        
        
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
            result_file_Vs_1000ms.open((Model_type+'_'+bindscheme+'_'+"CV_v_BCL_dt" + dt_val +"_dx"+dx_val+"_Vs_1000ms.txt").c_str());
            
            result_file_Vs_300ms.open((Model_type+'_'+bindscheme+'_'+"CV_v_BCL_dt" + dt_val +"_dx"+dx_val+"_Vs_300ms.txt").c_str());
            
            result_file_480.open((Model_type+'_'+bindscheme+'_'+"CV_v_BCL_dt" + dt_val +"_dx"+dx_val+"_params_480s.txt").c_str());
            
            result_file_500.open((Model_type+'_'+bindscheme+'_'+"CV_v_BCL_dt" + dt_val +"_dx"+dx_val+"_params_500s.txt").c_str());
            
            result_file.open((Model_type+'_'+bindscheme+'_'+"CV_v_BCL_dt" + dt_val +"_dx"+dx_val+"_params_cell60.txt").c_str());
            
        } else {
            result_file_Vs_1000ms.open((Model_type+'_'+bindscheme+'_'+"CV_v_BCL_dt" + dt_val +"_dx"+dx_val+"_kon_" + kon_val+ "_kd10e-6_Drug20e-6_Vs_1000ms.txt").c_str());
            
            result_file_Vs_300ms.open((Model_type+'_'+bindscheme+'_'+"CV_v_BCL_dt" + dt_val +"_dx"+dx_val+"_kon_" + kon_val+ "_kd10e-6_Drug20e-6_Vs_300ms.txt").c_str());
            
            result_file_480.open((Model_type+'_'+bindscheme+'_'+"CV_v_BCL_dt" + dt_val +"_dx"+dx_val+"_kon_" + kon_val+ "_kd10e-6_Drug20e-6_params_480s.txt").c_str());
            
            result_file_500.open((Model_type+'_'+bindscheme+'_'+"CV_v_BCL_dt" + dt_val +"_dx"+dx_val+"_kon_" + kon_val+ "_kd10e-6_Drug20e-6_params_500s.txt").c_str());
            
            result_file.open((Model_type+'_'+bindscheme+'_'+"CV_v_BCL_dt" + dt_val +"_dx"+dx_val+"_kon_" + kon_val+ "_kd10e-6_Drug20e-6_params_cell60.txt").c_str());
        }
        
        //Loading in initial conditions for each cell in the array.  These initial conditions are close to steady state.  The cable is allowed to rest for 10s and then stimulated 500 times to ensure the system approaches dynamic steady state.
        for (n=0; n<nCell; n+=1) {
                
            Cell[n].Cell_type = 3; // setting cell type to Epi
            Cell[n].V=-86.2;
            Cell[n].V_new=-86.2;
            Cell[n].dV=0;
            Cell[n].Na_in = 7.67;
            Cell[n].K_in = 138.3;
            Cell[n].Ca_in = 0.00007;
            Cell[n].Ca_sr = 1.3;
            Cell[n].Ca_ss = 0.00007;
            Cell[n].Ca_in_buffer = 0;
            Cell[n].Ca_ss_buffer = 0;
            Cell[n].Ca_sr_buffer = 0;
            Cell[n].m = 0.0;
            Cell[n].h = 0.75;
            Cell[n].j = 1.0;
            Cell[n].b = 0.0;
            
            /*Cell[n].mL =  0.00111859;
                Cell[n].hL =  0.339310414;*/
                
            Cell[n].d = 0 ;
            Cell[n].f = 1;
            Cell[n].f2 = 1;
            Cell[n].f_Ca = 1;
            Cell[n].r = 0;
            Cell[n].s = 1;
            Cell[n].xr1 = 0;
            Cell[n].xr2 = 1;
            Cell[n].xs = 0;
            Cell[n].OO = 0;
            Cell[n].R_bar = 1;
            
            Cell[n].I_Na = 0;
            Cell[n].I_Na_L = 0;
            Cell[n].I_Ca_L = 0;
            Cell[n].I_Kr = 0;
            Cell[n].I_Ks = 0;
            Cell[n].I_K1 = 0;
            Cell[n].I_Kp = 0;
            Cell[n].I_to = 0;
            Cell[n].I_Na_Ca = 0;
            Cell[n].I_Na_K = 0;
            Cell[n].I_p_Ca = 0;
            Cell[n].I_Ca_b = 0;
            Cell[n].I_Na_b = 0;
            Cell[n].I_stim = 0;
            Cell[n].I_tr = 0;
            Cell[n].I_leak = 0;
            Cell[n].I_up = 0;
            Cell[n].I_rel = 0;
            
            Cell[n].I_Na_ion_total = 0;
            Cell[n].I_K_ion_total = 0;
            Cell[n].I_Ca_ion_total = 0;
            Cell[n].I_total = 0.0;
            Cell[n].I_axial=0.0;
            
            Cell[n].peak_slope = 0;
            Cell[n].V_min = -88.654973;
            Cell[n].t_min = 0;
            Cell[n].V_thr = -88.654973;
            Cell[n].t_thr = 0;
            Cell[n].V_max = -88.654973;
            Cell[n].t_max = 0;
            Cell[n].V_90 = -88.654973;
            Cell[n].t_APD90 = 0;
            Cell[n].t_APD90_old = 0;
            Cell[n].dV_old = 0;
            Cell[n].APD90_flag = 0;
            
            Cell[n].cycle_num = 0;
            Cell[n].AP_flag = 0;
                
        }	// End of for n = 0 to nCell, loading in of initial conditions
        
        BCL = 1.0e3; //setting the BCL to 1000ms for first 500 pacings
        
        for (int ii = 1; ii <= 500; ii++) { //pacing 500 times with a BCL of 1000ms
            
            //Resets cell specific parameters every beat such as APD, DI, V90 etc.
            for (n = 0; n<nCell; n=n+1) {
                Calculate_Reset (&Cell[n], t);//resetting Cell properties that record AP characteristics
                
                if ((double) n < 0.1/dx)  {
                    Cell[n].V = 0.0; //establishing "stimulus" to start next propagated wave
                }
            }
            
            for (tTemp = 0.0; tTemp < BCL; tTemp=tTemp+dt){
                
#ifdef use_omp
#pragma omp parallel for
#endif
                for (n=0; n<nCell; n=n+1){
                    
                    Calculate_Points (Cell, t, n);  //Calculates t_min, V_min, t_max, V_max etc. etc.;
                    
                    //Updating current calculations for each cell
                    Calculate_I_Na(&Cell[n], tTemp , kon);
                    //Calculate_I_Na_L(&Cell[n], t );
                    Calculate_I_Ca_L(&Cell[n], tTemp );
                    Calculate_I_Kr(&Cell[n], tTemp );
                    Calculate_I_Ks(&Cell[n], tTemp );
                    Calculate_I_K1(&Cell[n], tTemp );
                    Calculate_I_Kp(&Cell[n], tTemp );
                    Calculate_I_to(&Cell[n], tTemp );
                    Calculate_I_Na_Ca(&Cell[n], tTemp );
                    Calculate_I_Na_K(&Cell[n], tTemp);
                    Calculate_I_p_Ca(&Cell[n], tTemp );
                    Calculate_I_Ca_b(&Cell[n], tTemp );
                    Calculate_I_Na_b(&Cell[n], tTemp );
                    Calculate_I_total(&Cell[n], tTemp, n);
                    
                    //Updating ionic concentrations
                    Calculate_Na_in(&Cell[n], tTemp );
                    Calculate_K_in(&Cell[n], tTemp );
                    Calculate_I_tr(&Cell[n], tTemp);
                    Calculate_I_leak(&Cell[n], tTemp);
                    Calculate_I_up(&Cell[n], tTemp);
                    Calculate_I_rel(&Cell[n], tTemp);
                    Calculate_Ca_sr(&Cell[n], tTemp );
                    Calculate_Ca_ss(&Cell[n], tTemp );
                    Calculate_Ca_in(&Cell[n], tTemp );
                    
                    // Update calculates axial currents and updates the voltage
                    Calculate_Update (Cell, tTemp, n); //Calculating new transmembrane potentials
                    
                    
                } // End of for n
                
                //Update the new voltages for each cell AFTER the parallelized "for n" loop.
                
                for (n=0; n<nCell; n=n+1){Cell[n].V = Cell[n].V_new;}
                
                t = t+dt; //update time
                
                if (fmod(t,1000)<dt) {std::cout << "t=" << t << ", Binding Scheme is " << Bind_Scheme << ", drug is " << Drug <<  std::endl;}
                
                if ((ii > 498) && (fmod(t,0.5)<dt)) {
                    for (n = 0; n < nCell; n++) {
                        result_file_Vs_1000ms << Cell[n].V << ", ";
                    }
                    result_file_Vs_1000ms << std::endl;
                }
                
            }
            
            /* Results section */
            
            if (ii == 480) {
                for (n = 0; n<nCell; n++) {
                    result_file_480 << n << ", " << Cell[n].CV_1cell << ", " << Cell[n].CV_2cell << ", " <<  Cell[n].CV_10cell << ", " << Cell[n].peak_slope << ", " << Cell[n].APD_90 << ", " << Cell[n].V_thr << ", " << Cell[n].cycle_num << std::endl; //This records the cell number, CV_1cell, CV_2cell, CV_10cell, peak upstroke velocity, and APD90 for each cell at the end of the 480th stimulus so that I can check that this is close to steady state
                }
                
            }else if (ii == 499){
                result_file << ii << ", " << BCL << ", " << Cell[59].CV_1cell << ", " << Cell[59].CV_2cell << ", " << Cell[59].CV_10cell << ", " << Cell[59].peak_slope << ", " << Cell[59].APD_90 << ", " << Cell[59].V_thr << ", " << Cell[59].cycle_num  << std::endl; //This records the stimulus number, BCL, CV_1cell, CV_2cell, peak upstroke velocity, and APD90
                
            }else if (ii == 500){
                
                result_file << ii << ", " << BCL << ", " << Cell[59].CV_1cell << ", " << Cell[59].CV_2cell << ", " << Cell[59].CV_10cell << ", " << Cell[59].peak_slope << ", " << Cell[59].APD_90 << ", " << Cell[59].V_thr << ", " << Cell[59].cycle_num  << std::endl; //This records the stimulus number, BCL, CV_1cell, CV_2cell, peak upstroke velocity, and APD90
                
                for (n = 0; n<nCell; n++) {
                    
                    result_file_500 << n << ", " << Cell[n].CV_1cell << ", " << Cell[n].CV_2cell << ", " <<  Cell[n].CV_10cell << ", " << Cell[n].peak_slope << ", " << Cell[n].APD_90 << ", " << Cell[n].V_thr << ", " << Cell[n].cycle_num << std::endl; //This records the cell number, CV_1cell, CV_2cell, peak upstroke velocity, and APD90 for each cell at the end of the 500th stimulus so that I can check that this is close to steady state
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
        
        for (BCL = 9.5e2; BCL > 2.51e2; BCL =  BCL - 50.0) {
            
            for (int ii = 1; ii <= 10; ii++) { //pacing 10 times for each new BCL
                
                //Resets cell specific parameters every beat such as APD, DI, V90 etc.
                for (n = 0; n<nCell; n=n+1) {
                    Calculate_Reset (&Cell[n], t);
                    
                    if ((double) n < 0.1/dx)  {
                        Cell[n].V = 0.0; //establishing "stimulus" to start next propagated wave
                    }
                }
                
                for (tTemp = 0.0; tTemp < BCL; tTemp=tTemp+dt){
                    
#ifdef use_omp
#pragma omp parallel for
#endif
                    for (n=0; n<nCell; n=n+1){
                        
                        Calculate_Points (Cell, t, n);  //Calculates t_min, V_min, t_max, V_max etc. etc.;
                        
                        //Updating current calculations for each cell
                        Calculate_I_Na(&Cell[n], tTemp , kon);
                        //Calculate_I_Na_L(&Cell[n], t );
                        Calculate_I_Ca_L(&Cell[n], tTemp );
                        Calculate_I_Kr(&Cell[n], tTemp );
                        Calculate_I_Ks(&Cell[n], tTemp );
                        Calculate_I_K1(&Cell[n], tTemp );
                        Calculate_I_Kp(&Cell[n], tTemp );
                        Calculate_I_to(&Cell[n], tTemp );
                        Calculate_I_Na_Ca(&Cell[n], tTemp );
                        Calculate_I_Na_K(&Cell[n], tTemp );
                        Calculate_I_p_Ca(&Cell[n], tTemp );
                        Calculate_I_Ca_b(&Cell[n], tTemp );
                        Calculate_I_Na_b(&Cell[n], tTemp );
                        Calculate_I_total(&Cell[n], tTemp, n);
                        
                        //Updating ionic concentrations
                        Calculate_Na_in(&Cell[n], tTemp );
                        Calculate_K_in(&Cell[n], tTemp );
                        Calculate_I_tr(&Cell[n], tTemp);
                        Calculate_I_leak(&Cell[n], tTemp);
                        Calculate_I_up(&Cell[n], tTemp);
                        Calculate_I_rel(&Cell[n], tTemp);
                        Calculate_Ca_sr(&Cell[n], tTemp );
                        Calculate_Ca_ss(&Cell[n], tTemp );
                        Calculate_Ca_in(&Cell[n], tTemp );
                        
                        // Update calculates axial currents and updates the voltage
                        Calculate_Update (Cell, tTemp, n); //Calculating new transmembrane potentials
                        
                        
                    } // End of for n
                    
                    //Update the new voltages for each cell AFTER the parallelized "for n" loop.
                    
                    for (n=0; n<nCell; n=n+1){Cell[n].V = Cell[n].V_new;}
                    
                    t = t+dt; //update time
                    
                    if (fmod(t,1000)<dt) {std::cout << "t=" << t << ", Binding Scheme is " << Bind_Scheme << ", drug is " << Drug <<  std::endl;}
                    
                    if ((BCL < 325.0) && (ii > 8) && (fmod(t,0.5)<dt)) {
                        for (n = 0; n < nCell; n++) {
                            result_file_Vs_300ms << Cell[n].V << ", ";
                        }
                        result_file_Vs_300ms << std::endl;
                    }
                    
                }
                
                //Results section
                if (ii > 8){
                    
                    result_file << ii << ", " << BCL << ", " << Cell[59].CV_1cell << ", " << Cell[59].CV_2cell << ", " << Cell[59].CV_10cell << ", " << Cell[59].peak_slope << ", " << Cell[59].APD_90 << ", " << Cell[59].V_thr << ", " << Cell[59].cycle_num  << std::endl; //This records the stimulus number, BCL, CV_1cell, CV_2cell, peak upstroke velocity, and APD90
                }
                
            } //end of loop through 10 stimuli
            
        } //end of loop through BCLs
        
        result_file.close();
        result_file_500.close();
        result_file_480.close();
        result_file_Vs_1000ms.close();
        result_file_Vs_300ms.close();//closing these files
        
    } //end of loop through binding schemes.
    
    result_file_time.open(("CV_v_BCL_dt" + dt_val +"_dx"+dx_val+"_timing.txt").c_str());
    
#ifdef use_omp
    double end = omp_get_wtime();
    result_file_time << "time: " << end - start << " seconds" << std::endl;
#else
    clock_t stop_s = clock();
    result_file_time << "time: " << ((double)stop_s - (double)start_s)/double(CLOCKS_PER_SEC) << " seconds" << std::endl;
#endif
    
    result_file_time.close();
    
    
}


#endif /* CV_vs_BCL_cable_h */
