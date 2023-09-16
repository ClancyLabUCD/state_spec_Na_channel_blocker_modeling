//
//  APDs_all_HH_form.h
//  Ten_Tusscher_model_with_drug
//
//  Created by Steffen Docken on 2/26/18.
//  Copyright © 2018 Steffen Docken (Lewis Lab). All rights reserved.
//  Published in the Journal of Theoretical Biology under the title "Rate-dependent effects of state-specific sodium channel blockers in cardiac tissue: Insights from idealized models"
//

#ifndef APDs_all_HH_form_h
#define APDs_all_HH_form_h

#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>
//#include <omp.h>
#include <unistd.h>
#include "Global_variables.h"
#include "Non_Na_currents.h"
#include "Ion_conc_dynamics.h"
#include "Na_currents.h"
#include "Full_cell_functions.h"
#include "load_IC.h"

//This program brings the cell to steady state and then paces it 500 times
void APDs_all_HH_form(){
    Cell_param Cell;
    int S1_cycle;
    
    double t, tTemp, interval;
    double BCL = 750.0; //Basic cycle length
    
    int BCL_int = BCL + 0.5; //This is converting BCL to an int (the +.5 is to make sure the truncation gives the right value) to then be changed to a string
    std::string BCL_val = std::to_string(BCL_int);
    
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
    

    
    chdir ("../Ten_Tusscher_model_with_drug_output/APDs_all_HH_forms_2_26_18");
    
    
    std::ofstream result_file, result_file0, result_file2, result_file3, result_file4, result_file5;
    
    double Diff_val = 1.0e2; // Diffusion values used for drug binding.
    std::string Diffval = "1e2";
    
    for(Bind_Scheme = 0; Bind_Scheme < 5; Bind_Scheme++){//Starting at Diff_counter ensured that the no drug model is not run twice.
        
        
        //const char *stateFileName = "Epi_Initial_10min.j";    //DONT TOUCH THIS FILE NAME
        //const char *stateFileName2 = "TEST.j";    //DONT TOUCH THIS FILE NAME
        
        //Initilizations
        std::string bindscheme;
        
        int counter = 0;
        
        
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
        
        Diffusion_drug = Diff_val;// Setting Diffusion for drug binding.
            
        std::string time_fac; //This will store the time factor in string form to be used in naming the file
            
        if (Bind_Scheme == 0) {
            result_file.open((Model_type+'_'+bindscheme+'_'+"drug_APs_BCL"+BCL_val+"_Stimdur06.txt").c_str());    // This file contains all of the voltages and other variables for each cell (the master-esque file) over the last 5 beats
        result_file0.open((Model_type+'_'+bindscheme+'_'+"drug_APs_BCL"+BCL_val+"_Stimdur06_param.txt").c_str()   );    // These files hold data for individual cells to compare V_min, V_max, V_90 etc. to compare between cells
        result_file2.open((Model_type+'_'+bindscheme+'_'+"drug_APs_BCL"+BCL_val+"_Stimdur06_state.txt").c_str()   );    // This file contains all of the individual Na channel states
        result_file3.open((Model_type+'_'+bindscheme+'_'+"drug_APs_BCL"+BCL_val+"_Stimdur06_IC_check.txt").c_str()   );    // This file contains all the cell variables at two points to check if steady state was reached.
        result_file4.open((Model_type+'_'+bindscheme+'_'+"drug_APs_BCL"+BCL_val+"_Stimdur06_initial_APs.txt").c_str()   );    // This file contains all the variables in the result_file, but for the first 5 APs
        result_file5.open((Model_type+'_'+bindscheme+'_'+"drug_APs_BCL"+BCL_val+"_Stimdur06_upstroke.txt").c_str()   );    // This file contains all the variables in the result_file, but for the first 5 APs
            
        } else {
            result_file.open((Model_type+'_'+bindscheme+'_'+"drug_APs_Diffval_"+Diffval+"_kd10e-6_Drug20e-6_BCL" +BCL_val+"_Stimdur06.txt").c_str());   // This file contains all of the voltages for each cell (the master-esque file)
            result_file0.open((Model_type+'_'+bindscheme+'_'+"drug_APs_Diffval_"+Diffval+"_kd10e-6_Drug20e-6_BCL" +BCL_val+"_Stimdur06_param.txt").c_str()   );    // These files hold data for individual cells to compare V_min, V_max, V_90 etc. to compare between cells
            result_file2.open((Model_type+'_'+bindscheme+'_'+"drug_APs_Diffval_"+Diffval+"_kd10e-6_Drug20e-6_BCL" +BCL_val+"_Stimdur06_state.txt").c_str()   );    // This file contains all of the individual Na channel states
            result_file3.open((Model_type+'_'+bindscheme+'_'+"drug_APs_Diffval_"+Diffval+"_kd10e-6_Drug20e-6_BCL" +BCL_val+"_Stimdur06_IC_check.txt").c_str()   );    // This file contains all the cell variables at two points to check if steady state was reached.
            result_file4.open((Model_type+'_'+bindscheme+'_'+"drug_APs_Diffval_"+Diffval+"_kd10e-6_Drug20e-6_BCL" +BCL_val+"_Stimdur06_initial_APs.txt").c_str()   );    // This file contains all the variables in the result_file, but for the first 50 APs
            result_file5.open((Model_type+'_'+bindscheme+'_'+"drug_APs_Diffval_"+Diffval+"_kd10e-6_Drug20e-6_BCL" +BCL_val+"_Stimdur06_upstroke.txt").c_str()   );    // This file contains all the variables in the result_file, but for the first 50 APs
        }
            
        result_file.precision(16);
            
        Cell.Cell_type = 3;        // Cell type = (1) for endo; (2) for M,  (3) epi
            
        t = 0.0;
        tTemp = 0.0; //initializing t and tTemp
            
        load_IC(&Cell);
            
        /*********************************************************************************************************************************
            Begin Time Loop Here
            *********************************************************************************************************************************/
            
        for (S1_cycle = 0; S1_cycle <=500; S1_cycle = S1_cycle + 1){                //Enter the number of beats that you want to simulate
                
            if (S1_cycle ==0){interval=waitTime;}
                
            else {interval=BCL;}                                                    //Enter the BCL of the simulation here
                
            //Resets cell specific parameters every beat such as APD, DI, V90 etc.
            Calculate_Reset (&Cell, t, tTemp);
            for (tTemp = 0; tTemp <= interval; tTemp=tTemp+dt){
                Calculate_Points (&Cell, t, tTemp); //Check points before updating to the next time point
                    
                //Updating current calculations for each cell
                Calculate_I_Na(&Cell, t, tTemp, Diffusion_drug); //This is where the Na channel model is chosen
                //Calculate_I_Na_L(&Cell, t, tTemp ); //Used if heartfailure model included
                Calculate_I_Ca_L(&Cell, t, tTemp );
                Calculate_I_Kr(&Cell, t, tTemp );
                Calculate_I_Ks(&Cell, t, tTemp );
                Calculate_I_K1(&Cell, t, tTemp );
                Calculate_I_Kp(&Cell, t, tTemp );
                Calculate_I_to(&Cell, t, tTemp );
                Calculate_I_Na_Ca(&Cell, t, tTemp );
                Calculate_I_Na_K(&Cell, t, tTemp );
                Calculate_I_p_Ca(&Cell, t, tTemp );
                Calculate_I_Ca_b(&Cell, t, tTemp );
                Calculate_I_Na_b(&Cell, t, tTemp );
                Calculate_I_total(&Cell, t, tTemp);
                    
                //Updating ionic concentrations
                Calculate_Na_in(&Cell, t, tTemp );
                Calculate_K_in(&Cell, t, tTemp );
                Calculate_I_tr(&Cell, t, tTemp);
                Calculate_I_leak(&Cell, t, tTemp);
                Calculate_I_up(&Cell, t, tTemp);
                Calculate_I_rel(&Cell, t, tTemp);
                Calculate_Ca_sr(&Cell, t, tTemp );
                Calculate_Ca_ss(&Cell, t, tTemp );
                Calculate_Ca_in(&Cell, t, tTemp );
                    
                //Calculates t_min, V_min, t_max, V_max etc. etc.; Update calculates axial currents and updates the voltage
                Calculate_Update (&Cell, t, tTemp);
                    
                    
                t = t+dt; //update time
                    
                /*if (sim_type ==0){                                        //This will rewrite state file name for initial conditions (a very long hold time: 30,000ms)
                    if (S1_cycle ==0 && tTemp >= (interval-dt)  ){
                    FILE *fp = fopen (stateFileName, "w" );
                     fwrite (&Cell, sizeof(Cell), 1, fp);
                     fclose (fp);
                     exit(0);}
                     }
                     
                     if (sim_type ==1){                                        //This will rewrite state file name2 for initial conditions with drug
                     if (S1_cycle ==5000 && tTemp >= (interval-dt)  ){
                     FILE *fp = fopen (stateFileName2, "w" );
                     fwrite (&Cell, sizeof(Cell), 1, fp);
                     fclose (fp);
                     exit(0);}
                }*/
                    
                /*******************************************************************************************************************************************************
                     Beginning of result files
            *******************************************************************************************************************************************************/
                    
                    
                    
                //This file is for individual cells that calculates all of the min, max, V_90 parameters.
                    
                if (counter%10==0 && S1_cycle >=495 ){result_file << t << ", " << Cell.V << ", " << Cell.I_Na<< ", " << Cell.I_Ca_L << ", " << Cell.I_Na_L << ", "<< Cell.I_up << ", "<< Cell.I_leak << ", "  << Cell.I_Ks << ", " << Cell.Ca_in << ", " << Cell.Ca_ss << ", " << Cell.m << ", " << Cell.h << ", " << Cell.j << ", " << Cell.b  << ", " << Cell.i_O << std::endl;}
                    
                if (counter%1000==0 && S1_cycle >=0){result_file2 << t << ", " << Cell.m << ", " << Cell.h << ", " << Cell.j << ", " << Cell.b << ", " << Cell.i_O << std::endl;}
                    
                if (Cell.flag2==1) {
                    result_file0 << Cell.t_min << ", " << Cell.V_min << ", " << Cell.t_thr<< ", " << Cell.V_thr << ", " << Cell.t_max << ", " << Cell.V_max << ", "<< Cell.t_90 << ", " << Cell.V_90 << ", "<< Cell.t_EAD << ", " << Cell.V_EAD << ", "
                        << Cell.t_EAD2 << ", "<<Cell.V_EAD2<< ", " << Cell.L_EAD << ", " << Cell.APD_90 << ", " << Cell.DI << ", " << S1_cycle<< ", "<< Cell.peak_slope<< ", " << Cell.t_I_Na_peak << ", " << Cell.I_Na_peak << ", " << Cell.t_I_Na_L_peak << ", " << Cell.I_Na_L_peak <<  std::endl;
                        
                    Cell.flag2 = 2;
                }
                    
                if (((interval - tTemp)<dt) && S1_cycle > 0) {
                    result_file5 << S1_cycle << ", " << Cell.peak_slope << std::endl;
                }
                    
                if (((t>=0.99*waitTime) && (t<0.99*waitTime+dt)) || ((t>=waitTime-1.0) && (t<waitTime-1.0+dt))) {//Last data collection is at waitTime - 1 to ensure the stimulus current has not been added yet.
                    result_file3 << Cell.V << ", " << Cell.dV << ", " << Cell.Na_in<< ", " << Cell.K_in << ", " << Cell.Ca_in << ", " << Cell.Ca_sr << ", "<< Cell.Ca_ss << ", " << Cell.Ca_in_buffer << ", "<< Cell.Ca_ss_buffer << ", " << Cell.Ca_sr_buffer << ", "
                        << Cell.I_Na << ", "<<Cell.I_Na_L<< ", " << Cell.I_Ca_L << ", " << Cell.I_Kr << ", " << Cell.I_Ks << ", " << Cell.I_K1<< ", "<< Cell.I_Kp << ", " << Cell.I_to << ", " << Cell.I_Na_Ca << ", " << Cell.I_Na_K << ", " << Cell.I_p_Ca << ", " << Cell.I_Ca_b << ", " << Cell.I_Na_b << ", " << Cell.I_Na_ion_total << ", " << Cell.I_Ca_ion_total << ", " << Cell.I_K_ion_total << ", " << Cell.I_total << ", " << Cell.I_tr << ", " << Cell.I_leak << ", " << Cell.I_up << ", " << Cell.I_rel << ", " << Cell.m << ", " << Cell.h << ", " << Cell.j << ", " << Cell.d << ", " << Cell.f << ", " << Cell.f2 << ", " << Cell.f_Ca << ", " << Cell.r << ", " << Cell.s << ", " << Cell.xr1 << ", " << Cell.xr2 << ", " << Cell.xs << ", " << Cell.OO << ", " << Cell.R_bar << ", " << Cell.b << std::endl;
                        
                }
                    
                if ((counter%10==0) && (S1_cycle <=50) && (S1_cycle > 0)){result_file4 << t << ", " << Cell.V << ", " << Cell.I_Na<< ", " << Cell.I_Ca_L << ", " << Cell.I_Na_L << ", "<< Cell.I_up << ", "<< Cell.I_leak << ", "  << Cell.I_Ks << ", " << Cell.Ca_in << ", " << Cell.Ca_ss << ", " << Cell.m << ", " << Cell.h << ", " << Cell.j << ", " << Cell.b << std::endl;}// End of individual result files
                    
                    
                /*******************************************************************************************************************************************************
                     End of result files
                    *******************************************************************************************************************************************************/
                    
                counter++;
                    
                if (fmod(t,5000)<dt) {std::cout << "S1_cycle=" <<S1_cycle<< ", " << "t=" << t << ", " "interval =" << interval << ", " << "Drug = "<< Drug << ", " << "Mutant = "<< mutant << std::endl;}
                    
            }//End of For tTemp ==
            tTemp = 0;
                
        }//End of for S1
        S1_cycle = 1;    //Resets S1 cycle number
        
            
        result_file.close();
        result_file0.close();
        result_file2.close();
        result_file3.close();
        result_file4.close();
        result_file5.close();
    }
    
}

#endif /* APDs_all_HH_form_h */
