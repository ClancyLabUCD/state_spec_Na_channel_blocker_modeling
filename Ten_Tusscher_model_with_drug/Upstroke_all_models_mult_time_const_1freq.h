//
//  Upstroke_all_models_mult_time_const_1freq.h
//  Ten_Tusscher_model_with_drug
//
//  Created by Steffen Docken on 12/6/17.
//  Copyright © 2017 Steffen Docken (Lewis Lab). All rights reserved.
//  Published in the Journal of Theoretical Biology under the title "Rate-dependent effects of state-specific sodium channel blockers in cardiac tissue: Insights from idealized models"
//

#ifndef Upstroke_all_models_mult_time_const_1freq_h
#define Upstroke_all_models_mult_time_const_1freq_h

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

//This program brings the cell to steady state and then paces it 500 times
void Upstroke_all_models_mult_time_const_1freq(){
    Cell_param Cell; //initializing the structure that will hold model variables
    int S1_cycle; //initializing counter for number of stimuli
    double BCL = 1000.0;
    
    int BCL_int = BCL + 0.5; //This is converting BCL to an int (the +.5 is to make sure the truncation gives the right value) to then be changed to a string
    std::string BCL_val = std::to_string(BCL_int);
    
    double t, tTemp, interval;//initializing variable to hold total simulated time, simulated time within a for loop, and the BCL value

    std::string bindscheme; //The string that will hold the name of the Binding Scheme being simulated
    
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
    }//Model_type is set in Global_variables.h
    
    chdir ("/Users/z3525459/Library/CloudStorage/OneDrive-UNSW/Cardiac_research/Drug_interaction_paper_2017/JTB_2022_paper_C++_codes/Ten_Tusscher_model_with_drug_output/Upstroke_vs_tau_b_BCL1000_12_7_17"); //setting the folder to save results
    
    std::ofstream result_file;
    
    double Diffusion_init = 1.0e5; //Different Diffusion values used for drug binding.
    
    double time_factor[22] = {100.0, 50.0, 30.0, 10.0, 5.0, 3.0, 1.0, 0.5, 0.3, 0.1, 0.05, 0.03, 0.01, 5.0e-3, 3.0e-3, 1.0e-3, 5.0e-4, 3.0e-4, 1.0e-4, 5.0e-4, 3.0e-4, 1.0e-4};  //Factor Diffusion will be multiplied by for different loops
    
    
    for(Bind_Scheme = 0; Bind_Scheme < 9; Bind_Scheme++){//Looping through all binding schemes
        
        t = 0.0;
        tTemp = 0.0; //initializing t and tTemp
        
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
        
        if (Bind_Scheme == 0) {
            result_file.open(Model_type+'_'+bindscheme+'_'+"drug_Upstroke_BCL"+BCL_val+"_Stimdur06.txt");
        } else {
            result_file.open(Model_type+'_'+bindscheme+'_'+"drug_Upstroke_BCL"+BCL_val+"_koninit_1e5_kd10e-6_Drug20e-6_Stimdur06.txt");
        }
        
        for (int Diff_counter = 0; Diff_counter < 22; Diff_counter++) {
            Diffusion_drug = Diffusion_init*time_factor[Diff_counter];// Setting Diffusion for drug binding.
            
            if (Bind_Scheme == 0){
                Diff_counter = 22; //Ends loop after one iteration if no drug is present.
            }
            
            
            Cell.Cell_type = 3;        // Cell type = (1) for endo; (2) for M,  (3) epi
            
            load_IC(&Cell); //loads the initial conditions
            
            /*********************************************************************************************************************************
             Begin Time Loop Here
             *********************************************************************************************************************************/
            
            for (S1_cycle = 0; S1_cycle <=500; S1_cycle = S1_cycle + 1){                //Enter the number of beats that you want to simulate
                
                if (S1_cycle ==0){interval=waitTime;}//allows the model to equilibrate for the waitTime set in Global_variables.h
                
                else {interval=BCL;}                                                    //sets the BL
                
                //Resets cell specific parameters every beat such as APD, DI, V90 etc.
                Calculate_Reset (&Cell, t, tTemp);
                for (tTemp = 0; tTemp <= interval; tTemp=tTemp+dt){
                    Calculate_Points (&Cell, t, tTemp); //Check points before updating to the next time point
                    
                    //Updating current calculations for each cell
                    Calculate_I_Na(&Cell, t, tTemp, Diffusion_drug); //This is where the Na channel model is chosen
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
                    
                    Calculate_Update (&Cell, t, tTemp);//Calculate new V
                    
                    
                    t = t+dt; //update time
                    
                    
                    /*******************************************************************************************************************************************************
                     Beginning of result files
                     *******************************************************************************************************************************************************/
                    
                    
                    if (((interval - tTemp)<dt) && S1_cycle > 498) {
                        result_file << Diffusion_drug << "," << BCL << ", " << S1_cycle << ", " << Cell.peak_slope << ", " << Cell.APD_90 << std::endl;
                    }
                    
                    
                    /*******************************************************************************************************************************************************
                     End of result files
                     *******************************************************************************************************************************************************/
                    
                    counter++;
                    
                    if (fmod(t,5000)<dt) {std::cout << "S1_cycle=" <<S1_cycle<< ", " << "t=" << t << ", " "interval =" << interval << ", " << "Drug = "<< Drug << ", " << "Mutant = "<< mutant << std::endl;}
                    
                }//End of For tTemp ==
                
            }//End of for S1
            
            
            
        }
        
        result_file.close();
    }
}

#endif /* Upstroke_all_models_mult_time_const_1freq_h */
