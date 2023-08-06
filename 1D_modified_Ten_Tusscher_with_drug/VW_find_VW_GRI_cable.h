//
//  VW_find_VW_GRI_cable.h
//  1D_modified_Ten_Tusscher_with_drug_2017_Summer
//
//  Created by Steffen Docken on 9/21/17.
//  Copyright © 2017 Steffen Docken (Lewis Lab). All rights reserved.
//  Published in the Journal of Theoretical Biology under the title "Rate-dependent effects of state-specific sodium channel blockers in cardiac tissue: Insights from idealized models"
//

#ifndef VW_find_VW_GRI_cable_h
#define VW_find_VW_GRI_cable_h

// This code is the same as VW_find_VW.h, except it is specifically for GRI, and therefore the times for S2 are much later, to account for the much slower conduction time in GRI

// This code takes the initial conditions from VW_initialization_cable.h and varies the time of S2 to find the size of the VW


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

void VW_find_VW_GRI_cable(){
#ifdef use_omp
    double start = omp_get_wtime();
#else
    clock_t start_s = clock();
#endif   //initial time for either parallelized or non-parallelized code.
    
    simState Cell_struct;
    simState *C_ptr = &Cell_struct;
    simState Cell_struct_initial;
    simState *C_ptr_initial = &Cell_struct_initial;//will hold the initial conditions of the Cells
    int n;
    
    //initialization
    double BCL;
    const double Sim_len = 2.5e3; //length of time simulated
    
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
    
    std::string dt_val, dx_val, BCL_val;  //These strings will be used to write the file names for the various results files.
    
    if (dt == 0.001) {
        dt_val = "1eneg3";
    }else if (dt == 0.0005){
        dt_val = "5eneg4";
    }
    
    if (dx == 0.01) {
        dx_val = "1eneg2";
    }else if (dx == 0.005){
        dx_val = "5eneg3";
    }
    
    double BCL_vec[3] = {1.0e3, 6.75e2, 4.0e2}; //list of BCLs to be used
    
    std::ofstream result_file_previous_Vs, result_file_Vs, result_file_time, result_file_S1, result_file_S2, result_file_S2_times, error_file;  // The first two files contain all of the voltages for each cell every 2ms for the t_S2 just before and just after the borders of the Vulnerable Window, the third file is for the total time of the simulation, the next two files contain the various parameters for each cell following the S1 and S2 stimuli, the second to last file contains the times of S2 for one way and two way propagation, and the last file contains any errors that are output.
    FILE *f_initial_conditions;  //pointer to file for loading initial conditions
    
    error_file.open(("errors_VW_find_VW_GRI_dt" + dt_val +"_dx"+dx_val+".txt").c_str());  //creates and opens the file to store any errors that occur when loading the initial conditions
    
    
    //Initilizations
    double t_S2;
    const int V_num = (int) (Sim_len + 0.5)/5 + 1; //This is the number of time points at which Vs should saved
    std::string bindscheme;
    double V_mat[V_num][nCell+1] = {};
    double prev_V_mat[V_num][nCell+1] = {};//These arrays will hold the V values (and corresponding t as the first entry of each row) from the current and previous loops of the program
    double params_S1[nCell][8] = {};
    double params_S2[nCell][8] = {}; //These arrays will hold the parameters describing conduction following S1 and S2
    int flag_retro_prop1, flag_retro_prop2, flag_ante_prop1, flag_ante_prop2;
    int mat_index; //This int will increase to make sure the Vs recorded every 2ms aren't ever recorded in the same row of the V_mat's.
    
    for (int jj = 0; jj < 2; jj++) {  //Looping over the BCL values for which the VW will be checked.
        
        BCL = BCL_vec[jj]; //BCL that has been used for setting initial conditions
        
        if (jj == 0){
            BCL_val = "1000";
        }else if(jj == 1){
            BCL_val = "675";
        }else if(jj == 2){
            BCL_val = "400";
        }
        
        Bind_Scheme = 1; //fixing Bind_Scheme
        
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
        }
        
        //chdir ("/Users/steffendocken/Documents/Research/2017_Summer/Ten_Tusscher_1D_sims_output/cable/VW/VW_initialization_9_27_17/"); //for my laptop (changes the directory from which the iniitial conditions will be loaded.  When this code is run on linux system, all initial condition files must be in the same folder as where this program is saved.)
        
        if (Bind_Scheme == 0) {
            
            f_initial_conditions = fopen((Model_type+'_'+bindscheme+'_'+"VW_initialization_BCL_"+BCL_val+"ms_dt" + dt_val +"_dx"+dx_val+"_state_vars.init_cond").c_str(), "r"); //opens the file from which the initial conditions will be read
            
            result_file_S2_times.open((Model_type+'_'+bindscheme+'_'+"VW_find_VW_BCL_" + BCL_val+ "ms_dt" + dt_val +"_dx"+dx_val+"_VW_times.txt").c_str()); //creates and opens the file to store the times for the borders of the vulnerable window.  Organizationally, it makes more sense to put this file in the results folder, as opposed to the initial conditions folder (which is what would happen if chdir is not commented out).  However, when this is run on a linux system and chdir is commented out, all result files (including this one) will be saved in the same folder as this program, and organization of files can be done later.
            
        } else {
            
            f_initial_conditions = fopen((Model_type+'_'+bindscheme+'_'+"VW_initialization_BCL_"+BCL_val+"ms_dt" + dt_val +"_dx"+dx_val+"_kon_" + kon_val+ "_kd10e-6_Drug20e-6_state_vars.init_cond").c_str(), "r");  //opens the file from which the initial conditions will be read
            
            result_file_S2_times.open((Model_type+'_'+bindscheme+'_'+"VW_find_VW_BCL_" + BCL_val+ "ms_dt" + dt_val +"_dx"+dx_val+"_kon_" + kon_val+ "_kd10e-6_Drug20e-6_VW_times.txt").c_str());  //creates and opens the file to store the times for the borders of the vulnerable window.  Organizationally, it makes more sense to put this file in the results folder, as opposed to the initial conditions folder (which is what would happen if chdir is not commented out).  However, when this is run on a linux system and chdir is commented out, all result files (including this one) will be saved in the same folder as this program, and organization of files can be done later.
        }
        
        if (f_initial_conditions == NULL) {
            error_file << "Could not load initial condiditions for Model " << Model_type << " with " << bindscheme << " and BCL = " << BCL_val << std::endl;  //saving error if the initial conditions folder is not opened correctly.
        } else {
            
            fread(C_ptr_initial, sizeof(simState), 1, f_initial_conditions);
            fclose(f_initial_conditions); //establishing initial conditions
        }
        
        
        
        //chdir ("/Users/steffendocken/Documents/Research/2017_Summer/Ten_Tusscher_1D_sims_output/sandbox/cable/VW_test/"); //for my laptop (changes the directory from which the iniitial conditions will be loaded.  When this code is run on linux system, all initial condition files must be in the same folder as where this program is saved.)
        
        
        flag_retro_prop1 = 0;
        flag_retro_prop2 = 0;
        flag_ante_prop1 = 0;
        flag_ante_prop2 = 0;//These flags will be switched to 1 for when retro and antegrade propagation occur
        
        for (t_S2 = 715.0; t_S2 < 1000.0; t_S2 = t_S2 + 0.5) {
            flag_retro_prop1 = 0;
            flag_ante_prop1 = 0; //resetting flags that indicate that the first wave has passed
            
            Cell_struct = Cell_struct_initial; //copying initial state over
            
            flag_S2 = 0; //Declaring that an early stimulus should not be applied yet.
            mat_index = 0; //resetting the value to keep track of which row to add V values to
            
            //Resets cell specific parameters every beat such as APD, DI, V90 etc.
            for (n = 0; n<nCell; n=n+1) {
                Calculate_Reset (&(C_ptr->cellData[n]), C_ptr->t); //resetting Cell properties that record AP characteristics
                
                if ((double) n < 0.1/dx)  {
                    C_ptr->cellData[n].V = 0.0; //establishing "stimulus" to start next propagated wave
                }
            }
            
            C_ptr->beat = C_ptr->beat + 1; //updating the number of stimuli applied to the system.
            
            for (C_ptr->tTemp = 0.0; C_ptr->tTemp < t_S2 - dt/2; (C_ptr->tTemp) = (C_ptr->tTemp) +dt){  //Running the model upto the point S2 is applied.   The "-dt/2" is included to ensure the t_S2 time point is not included in this loop (even if dt and t_S2 with double precision are slightly off)
                
#ifdef use_omp
#pragma omp parallel for
#endif
                for (n=0; n<nCell; n=n+1){
                    
                    Calculate_Points (C_ptr->cellData, C_ptr->t, n);  //Calculates t_min, V_min, t_max, V_max etc. etc.;
                    
                    //Updating current calculations for each cell
                    Calculate_I_Na(&(C_ptr->cellData[n]), C_ptr->tTemp , kon);
                    //Calculate_I_Na_L(&(C_ptr->cellData[n]), t );
                    Calculate_I_Ca_L(&(C_ptr->cellData[n]), C_ptr->tTemp );
                    Calculate_I_Kr(&(C_ptr->cellData[n]), C_ptr->tTemp );
                    Calculate_I_Ks(&(C_ptr->cellData[n]), C_ptr->tTemp );
                    Calculate_I_K1(&(C_ptr->cellData[n]), C_ptr->tTemp );
                    Calculate_I_Kp(&(C_ptr->cellData[n]), C_ptr->tTemp );
                    Calculate_I_to(&(C_ptr->cellData[n]), C_ptr->tTemp );
                    Calculate_I_Na_Ca(&(C_ptr->cellData[n]), C_ptr->tTemp );
                    Calculate_I_Na_K(&(C_ptr->cellData[n]), C_ptr->tTemp);
                    Calculate_I_p_Ca(&(C_ptr->cellData[n]), C_ptr->tTemp );
                    Calculate_I_Ca_b(&(C_ptr->cellData[n]), C_ptr->tTemp );
                    Calculate_I_Na_b(&(C_ptr->cellData[n]), C_ptr->tTemp );
                    Calculate_I_total(&(C_ptr->cellData[n]), C_ptr->tTemp, n); //Calculating Currents
                    
                    //Updating ionic concentrations
                    Calculate_Na_in(&(C_ptr->cellData[n]), C_ptr->tTemp );
                    Calculate_K_in(&(C_ptr->cellData[n]), C_ptr->tTemp );
                    Calculate_I_tr(&(C_ptr->cellData[n]), C_ptr->tTemp);
                    Calculate_I_leak(&(C_ptr->cellData[n]), C_ptr->tTemp);
                    Calculate_I_up(&(C_ptr->cellData[n]), C_ptr->tTemp);
                    Calculate_I_rel(&(C_ptr->cellData[n]), C_ptr->tTemp);
                    Calculate_Ca_sr(&(C_ptr->cellData[n]), C_ptr->tTemp );
                    Calculate_Ca_ss(&(C_ptr->cellData[n]), C_ptr->tTemp );
                    Calculate_Ca_in(&(C_ptr->cellData[n]), C_ptr->tTemp ); //Updating Ion concentrations
                    
                    // Update calculates axial currents and updates the voltage
                    Calculate_Update (C_ptr->cellData, C_ptr->tTemp, n); //Calculating new transmembrane potentials
                    
                    
                } // End of for n
                
                //Update the new voltages for each cell AFTER the parallelized "for n" loop.
                
                for (n=0; n<nCell; n=n+1){C_ptr->cellData[n].V = C_ptr->cellData[n].V_new;}
                
                if (fmod(C_ptr->tTemp + dt/2,5)<dt) {//adding dt/2 to allow tTemp to be between dt/2 less than an even number and dt/2 more than an even number
                    prev_V_mat[mat_index][0] = V_mat[mat_index][0];
                    
                    V_mat[mat_index][0] = C_ptr->tTemp;
                    
                    for (n = 0; n < nCell; n++) {
                        prev_V_mat[mat_index][n+1] = V_mat[mat_index][n+1]; //Saving the previous loop's V values
                        
                        V_mat[mat_index][n+1] = C_ptr->cellData[n].V; //Saving the current loop's V values
                    }
                    mat_index++;
                } //storing Vs
                
                C_ptr->t = C_ptr->t+dt; //update time
                
            }
            
            for (n = 0; n<nCell; n++) {
                params_S1[n][0] = n;
                params_S1[n][1] = C_ptr->cellData[n].CV_1cell;
                params_S1[n][2] = C_ptr->cellData[n].CV_2cell;
                params_S1[n][3] = C_ptr->cellData[n].CV_10cell;
                params_S1[n][4] = C_ptr->cellData[n].peak_slope;
                params_S1[n][5] = C_ptr->cellData[n].APD_90;
                params_S1[n][6] = C_ptr->cellData[n].V_thr;
                params_S1[n][7] = C_ptr->cellData[n].cycle_num; //This records the cell number, CV_1cell, CV_2cell, CV_10cell, peak upstroke velocity, and APD90 for each cell at the end of the S1 stimulus so that I can check that this is close to the values after the 500th stimulus in initialization
            }
            
            flag_S2 = 1; //Declaring that an early stimulus should be applied
            
            //Resets cell specific parameters every beat such as APD, DI, V90 etc.
            for (n = 0; n<nCell; n=n+1) {
                Calculate_Reset (&(C_ptr->cellData[n]), C_ptr->t);
            }
            
            C_ptr->beat = C_ptr->beat + 1; //Updating the number of stimuli applied to the system
            
            for (C_ptr->tTemp = 0.0; C_ptr->tTemp < Sim_len + dt/2 - t_S2; C_ptr->tTemp=C_ptr->tTemp+dt){ //Simulating the cable from the time of the S2 stimulus through the end of the simulation time.   The "+dt/2" is included to ensure the Sim_len time point is included in this loop (even if dt and t_S2 with double precision are slightly off)
                
#ifdef use_omp
#pragma omp parallel for
#endif
                for (n=0; n<nCell; n=n+1){
                    
                    Calculate_Points (C_ptr->cellData, C_ptr->t, n);  //Calculates t_min, V_min, t_max, V_max etc. etc.;
                    
                    //Updating current calculations for each cell
                    Calculate_I_Na(&(C_ptr->cellData[n]), C_ptr->tTemp , kon);
                    //Calculate_I_Na_L(&(C_ptr->cellData[n]), t );
                    Calculate_I_Ca_L(&(C_ptr->cellData[n]), C_ptr->tTemp );
                    Calculate_I_Kr(&(C_ptr->cellData[n]), C_ptr->tTemp );
                    Calculate_I_Ks(&(C_ptr->cellData[n]), C_ptr->tTemp );
                    Calculate_I_K1(&(C_ptr->cellData[n]), C_ptr->tTemp );
                    Calculate_I_Kp(&(C_ptr->cellData[n]), C_ptr->tTemp );
                    Calculate_I_to(&(C_ptr->cellData[n]), C_ptr->tTemp );
                    Calculate_I_Na_Ca(&(C_ptr->cellData[n]), C_ptr->tTemp );
                    Calculate_I_Na_K(&(C_ptr->cellData[n]), C_ptr->tTemp);
                    Calculate_I_p_Ca(&(C_ptr->cellData[n]), C_ptr->tTemp );
                    Calculate_I_Ca_b(&(C_ptr->cellData[n]), C_ptr->tTemp );
                    Calculate_I_Na_b(&(C_ptr->cellData[n]), C_ptr->tTemp );
                    Calculate_I_total(&(C_ptr->cellData[n]), C_ptr->tTemp, n); //Calculating Currents
                    
                    //Updating ionic concentrations
                    Calculate_Na_in(&(C_ptr->cellData[n]), C_ptr->tTemp );
                    Calculate_K_in(&(C_ptr->cellData[n]), C_ptr->tTemp );
                    Calculate_I_tr(&(C_ptr->cellData[n]), C_ptr->tTemp);
                    Calculate_I_leak(&(C_ptr->cellData[n]), C_ptr->tTemp);
                    Calculate_I_up(&(C_ptr->cellData[n]), C_ptr->tTemp);
                    Calculate_I_rel(&(C_ptr->cellData[n]), C_ptr->tTemp);
                    Calculate_Ca_sr(&(C_ptr->cellData[n]), C_ptr->tTemp );
                    Calculate_Ca_ss(&(C_ptr->cellData[n]), C_ptr->tTemp );
                    Calculate_Ca_in(&(C_ptr->cellData[n]), C_ptr->tTemp ); //Updating Ion concentrations
                    
                    // Update calculates axial currents and updates the voltage
                    Calculate_Update (C_ptr->cellData, C_ptr->tTemp, n);  //Calculating new transmembrane potentials
                    
                    
                } // End of for n
                
                //Update the new voltages for each cell AFTER the parallelized "for n" loop.
                
                for (n=0; n<nCell; n=n+1){C_ptr->cellData[n].V = C_ptr->cellData[n].V_new;}
                
                if (fmod(C_ptr->tTemp + dt/2 + t_S2,5)<dt) {//adding dt/2 to allow tTemp to be between dt/2 less than an even number and dt/2 more than an even number
                    prev_V_mat[mat_index][0] = V_mat[mat_index][0];
                    
                    V_mat[mat_index][0] = C_ptr->tTemp + t_S2;//adding t_S2, since tTemp restarts at 0 when S2 is given
                    
                    for (n = 0; n < nCell; n++) {
                        prev_V_mat[mat_index][n+1] = V_mat[mat_index][n+1]; //Saving the previous loop's V values
                        
                        V_mat[mat_index][n+1] = C_ptr->cellData[n].V; //Saving the current loop's V values
                    }
                    mat_index++; //increasing mat_index for next loop.
                } //storing Vs
                
                C_ptr->t = C_ptr->t+dt; //update time
                
                if ((flag_retro_prop1 == 0) && (C_ptr->cellData[39].V < - 60.0)) {
                    flag_retro_prop1 = 1; //updating flag that indicates that the first forward propagating wave has passed the beginning of the cable
                    
                }
                
                if ((flag_retro_prop1 == 1) && (flag_retro_prop2 == 0) && (C_ptr->cellData[39].V > - 40.0)) {
                    flag_retro_prop2 = 1; //updating flag that says retrograde propagation has occured.
                    
                }
                
                if ((jj == 0) && (C_ptr->tTemp + t_S2 > 895.0)) {//Based on simulation results, for a BCL of 1000ms, the wave resulting from the S1 stimulus is at Cell 480 895ms after S1, so this ensures that the program does not check that V of Cell 480 is less than -60 before the wave from S1 has reached cell 480
                    if ((flag_retro_prop2 == 2) && (flag_ante_prop1 == 0) && (C_ptr->cellData[479].V < -60.0)) {
                        flag_ante_prop1 = 1; //updating flag that indicates that the first forward propagating wave has passed
                    }
                } else if ((jj == 1) && (C_ptr->tTemp + t_S2 > 1145.0)) {//Based on simulation results, for a BCL of 675ms, the wave resulting from the S1 stimulus is at Cell 480 1145ms after S1, so this ensures that the program does not check that V of Cell 480 is less than -60 before the wave from S1 has reached cell 480
                    if ((flag_retro_prop2 == 2) && (flag_ante_prop1 == 0) && (C_ptr->cellData[479].V < -60.0)) {
                        flag_ante_prop1 = 1; //updating flag that indicates that the first forward propagating wave has passed
                    }
                }
                
                if ((flag_ante_prop1 == 1) && (flag_ante_prop2 == 0) && (C_ptr->cellData[479].V > -40.0)) {
                    flag_ante_prop2 = 1; //updating flag that indicates that a second forward propagating wave has passed
                }
                
            }
            
            if (mat_index > V_num){
                error_file << "Too many V time points were captured for Model " << Model_type << " with " << bindscheme << " and BCL = " << BCL_val << std::endl;  //saving error if something went wrong with saving the Vs, which means the arrays might get overloaded and alter other data values.
            }
            
            for (n = 0; n<nCell; n++) {
                params_S2[n][0] = n;
                params_S2[n][1] = C_ptr->cellData[n].CV_1cell;
                params_S2[n][2] = C_ptr->cellData[n].CV_2cell;
                params_S2[n][3] = C_ptr->cellData[n].CV_10cell;
                params_S2[n][4] = C_ptr->cellData[n].peak_slope;
                params_S2[n][5] = C_ptr->cellData[n].APD_90;
                params_S2[n][6] = C_ptr->cellData[n].V_thr;
                params_S2[n][7] = C_ptr->cellData[n].cycle_num; //This records the cell number, CV_1cell, CV_2cell, CV_10cell, peak upstroke velocity, and APD90 for each cell at the end of the S2 stimulus
            }
            
            std::cout << "t=" << C_ptr->t << ", Binding Scheme is " << Bind_Scheme << ", drug is " << Drug << ", BCL = " << BCL_val << ", S2 = " << t_S2 << std::endl;  //outputting progress
            
            if (flag_retro_prop2 == 1) {
                flag_retro_prop2 = 2; //updating flag so don't save Vs twice
                
                if (Bind_Scheme == 0) {
                    
                    result_file_previous_Vs.open((Model_type+'_'+bindscheme+'_'+"VW_find_VW_BCL_" + BCL_val+ "ms_dt" + dt_val +"_dx"+dx_val+"_1way_prop_previous_Vs.txt").c_str());
                    
                    result_file_Vs.open((Model_type+'_'+bindscheme+'_'+"VW_find_VW_BCL_" + BCL_val+ "ms_dt" + dt_val +"_dx"+dx_val+"_1way_prop_Vs.txt").c_str());
                    
                    result_file_S1.open((Model_type+'_'+bindscheme+'_'+"VW_find_VW_BCL_" + BCL_val+ "ms_dt" + dt_val +"_dx"+dx_val+"_1way_prop_params_S1.txt").c_str());
                    
                    result_file_S2.open((Model_type+'_'+bindscheme+'_'+"VW_find_VW_BCL_" + BCL_val+ "ms_dt" + dt_val +"_dx"+dx_val+"_1way_prop_params_S2.txt").c_str());
                    
                } else {
                    
                    result_file_previous_Vs.open((Model_type+'_'+bindscheme+'_'+"VW_find_VW_BCL_" + BCL_val+ "ms_dt" + dt_val +"_dx"+dx_val+"_kon_" + kon_val+ "_kd10e-6_Drug20e-6_1way_prop_previous_Vs.txt").c_str());
                    
                    result_file_Vs.open((Model_type+'_'+bindscheme+'_'+"VW_find_VW_BCL_" + BCL_val+ "ms_dt" + dt_val +"_dx"+dx_val+"_kon_" + kon_val+ "_kd10e-6_Drug20e-6_1way_prop_Vs.txt").c_str());
                    
                    result_file_S1.open((Model_type+'_'+bindscheme+'_'+"VW_find_VW_BCL_" + BCL_val+ "ms_dt" + dt_val +"_dx"+dx_val+"_kon_" + kon_val+ "_kd10e-6_Drug20e-6_1way_prop_params_S1.txt").c_str());
                    
                    result_file_S2.open((Model_type+'_'+bindscheme+'_'+"VW_find_VW_BCL_" + BCL_val+ "ms_dt" + dt_val +"_dx"+dx_val+"_kon_" + kon_val+ "_kd10e-6_Drug20e-6_1way_prop_params_S2.txt").c_str());
                    
                } //opening results files
                
                for (int ii = 0; ii < (int) (Sim_len + 0.5)/5 + 1; ii ++){
                    for (n = 0; n < nCell + 1; n++) { //entering nCell + 1 values at each time point, becasue first entry is the time, and the other nCell are the V's of each cell
                        result_file_previous_Vs << prev_V_mat[ii][n] << ", ";
                        
                        result_file_Vs << V_mat[ii][n] << ", ";
                    }
                    
                    result_file_previous_Vs << std::endl;
                    result_file_Vs << std::endl;
                    
                } //saving V's from S2 just before VW border
                
                for (n = 0; n<nCell; n++) {
                    result_file_S1 << params_S1[n][0] << ", " << params_S1[n][1] << ", " << params_S1[n][2] << ", " <<  params_S1[n][3] << ", " << params_S1[n][4] << ", " << params_S1[n][5] << ", " << params_S1[n][6] << ", " << params_S1[n][7] << std::endl; //This records the cell number, CV_1cell, CV_2cell, CV_10cell, peak upstroke velocity, and APD90 for each cell at the end of the S1 stimulus so that I can check that it has not changed from the initial conditions
                    
                    result_file_S2 << params_S2[n][0] << ", " << params_S2[n][1] << ", " << params_S2[n][2] << ", " <<  params_S2[n][3] << ", " << params_S2[n][4] << ", " << params_S2[n][5] << ", " << params_S2[n][6] << ", " << params_S2[n][7] << std::endl; //This records the cell number, CV_1cell, CV_2cell, CV_10cell, peak upstroke velocity, and APD90 for each cell at the end of the S2 stimulus
                }
                
                result_file_S2_times << t_S2 << std::endl; //This saves the smallest values of S2 that gives 1 way propagation
                
                result_file_previous_Vs.close();
                result_file_Vs.close();
                result_file_S1.close();
                result_file_S2.close();
                
            }
            
            if (flag_ante_prop2 == 1) {
                flag_ante_prop2 = 2; //updating flag so don't save Vs twice
                
                if (Bind_Scheme == 0) {
                    
                    result_file_previous_Vs.open((Model_type+'_'+bindscheme+'_'+"VW_find_VW_BCL_" + BCL_val+ "ms_dt" + dt_val +"_dx"+dx_val+"_2way_prop_previous_Vs.txt").c_str());
                    
                    result_file_Vs.open((Model_type+'_'+bindscheme+'_'+"VW_find_VW_BCL_" + BCL_val+ "ms_dt" + dt_val +"_dx"+dx_val+"_2way_prop_Vs.txt").c_str());
                    
                    result_file_S1.open((Model_type+'_'+bindscheme+'_'+"VW_find_VW_BCL_" + BCL_val+ "ms_dt" + dt_val +"_dx"+dx_val+"_2way_prop_params_S1.txt").c_str());
                    
                    result_file_S2.open((Model_type+'_'+bindscheme+'_'+"VW_find_VW_BCL_" + BCL_val+ "ms_dt" + dt_val +"_dx"+dx_val+"_2way_prop_params_S2.txt").c_str());
                    
                } else {
                    
                    result_file_previous_Vs.open((Model_type+'_'+bindscheme+'_'+"VW_find_VW_BCL_" + BCL_val+ "ms_dt" + dt_val +"_dx"+dx_val+"_kon_" + kon_val+ "_kd10e-6_Drug20e-6_2way_prop_previous_Vs.txt").c_str());
                    
                    result_file_Vs.open((Model_type+'_'+bindscheme+'_'+"VW_find_VW_BCL_" + BCL_val+ "ms_dt" + dt_val +"_dx"+dx_val+"_kon_" + kon_val+ "_kd10e-6_Drug20e-6_2way_prop_Vs.txt").c_str());
                    
                    result_file_S1.open((Model_type+'_'+bindscheme+'_'+"VW_find_VW_BCL_" + BCL_val+ "ms_dt" + dt_val +"_dx"+dx_val+"_kon_" + kon_val+ "_kd10e-6_Drug20e-6_2way_prop_params_S1.txt").c_str());
                    
                    result_file_S2.open((Model_type+'_'+bindscheme+'_'+"VW_find_VW_BCL_" + BCL_val+ "ms_dt" + dt_val +"_dx"+dx_val+"_kon_" + kon_val+ "_kd10e-6_Drug20e-6_2way_prop_params_S2.txt").c_str());
                    
                } //opening results files
                
                for (int ii = 0; ii < (int) (Sim_len + 0.5)/5 + 1; ii ++){
                    for (n = 0; n < nCell + 1; n++) {//entering nCell + 1 values at each time point, becasue first entry is the time, and the other nCell are the V's of each cell
                        result_file_previous_Vs << prev_V_mat[ii][n] << ", ";
                        
                        result_file_Vs << V_mat[ii][n] << ", ";
                    }
                    
                    result_file_previous_Vs << std::endl;
                    result_file_Vs << std::endl;
                    
                } //saving V's from just after VW border
                
                for (n = 0; n<nCell; n++) {
                    result_file_S1 << params_S1[n][0] << ", " << params_S1[n][1] << ", " << params_S1[n][2] << ", " <<  params_S1[n][3] << ", " << params_S1[n][4] << ", " << params_S1[n][5] << ", " << params_S1[n][6] << ", " << params_S1[n][7] << std::endl; //This records the cell number, CV_1cell, CV_2cell, CV_10cell, peak upstroke velocity, and APD90 for each cell at the end of the S1 stimulus so that I can check that it has not changed from the initial conditions
                    
                    result_file_S2 << params_S2[n][0] << ", " << params_S2[n][1] << ", " << params_S2[n][2] << ", " <<  params_S2[n][3] << ", " << params_S2[n][4] << ", " << params_S2[n][5] << ", " << params_S2[n][6] << ", " << params_S2[n][7] << std::endl; //This records the cell number, CV_1cell, CV_2cell, CV_10cell, peak upstroke velocity, and APD90 for each cell at the end of the S2 stimulus
                }
                
                result_file_S2_times << t_S2 << std::endl; //This saves the smallest value of S2 that gives 2 way propagation
                
                result_file_previous_Vs.close();
                result_file_Vs.close();
                result_file_S1.close();
                result_file_S2.close();
                
                break;
                
            }
            
        }
        
        result_file_S2_times.close(); //closes the file saving the times of the VW borders
        
        
    }
    
    result_file_time.open(("VW_find_VW_GRI_dt" + dt_val +"_dx"+dx_val+"_timing.txt").c_str());
    
#ifdef use_omp
    double end = omp_get_wtime();
    result_file_time << "time: " << end - start << " seconds" << std::endl;
#else
    clock_t stop_s = clock();
    result_file_time << "time: " << ((double)stop_s - (double)start_s)/double(CLOCKS_PER_SEC) << " seconds" << std::endl;
#endif //recording the time it takes to run the simulation.
    
    result_file_time.close();
    error_file.close(); //closing the file that records the amount of time the simulation takes and the file which records any errors that occur while loading initial conditions.
    
}

#endif /* VW_find_VW_GRI_cable_h */
