//
//  Upstroke_rest_30s_all_models_mult_time_const.h
//  Ten_Tusscher_model_with_drug
//
//  Created by Steffen Docken on 12/7/17.
//  Copyright © 2017 Steffen Docken (Lewis Lab). All rights reserved.
//  Published in the Journal of Theoretical Biology under the title "Rate-dependent effects of state-specific sodium channel blockers in cardiac tissue: Insights from idealized models"
//

#ifndef Upstroke_rest_30s_all_models_mult_time_const_h
#define Upstroke_rest_30s_all_models_mult_time_const_h

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

//This program brings the cell to steady state and then paces it twice at a BCL of 30s to approximate the peak upstroke velocity with an infinitely long BCL
void Upstroke_rest_30s_all_models_mult_time_const(){
    Cell_param Cell; //initializing the structure that will hold model variables
    int S1_cycle; //initializing counter for number of stimuli
    double rest = 3.0e4; // Length of BCL
    
    double t, tTemp, interval;//initializing variable to hold total simulated time, simulated time within a for loop, and the BCL value
    
    std::string bindscheme;//The string that will hold the name of the Binding Scheme being simulated
    
    int counter = 0;//counts the number of time steps that have been taken
    
    std::string Diffval; //This will store the time factor in string form to be used in naming the file
    
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
    
    chdir("/Users/z3525459/Library/CloudStorage/OneDrive-UNSW/Cardiac_research/Drug_interaction_paper_2017/JTB_2022_paper_C++_codes/Ten_Tusscher_model_with_drug_output/Upstroke_vs_Period_3_tau_b_12_5_17/"); //V-dependent results //("/Users/z3525459/Library/CloudStorage/OneDrive-UNSW/Cardiac_research/Drug_interaction_paper_2017/JTB_2022_paper_C++_codes/Ten_Tusscher_model_with_drug_output/Upstroke_vs_Period_nonVdepend_2_28_18/");// Non-V-dependent results 
    
    std::ofstream result_file;

    double Diffusion_vec[3] = {1.0e5, 1.0e3, 1.0e1}; //Different Diffusion values used for drug binding.
    
    
    for(Bind_Scheme = 0; Bind_Scheme < 9; Bind_Scheme++){//Looping through all binding schemes
        
        //Initilizations for file names
        
        
        counter = 0;
        
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
        
        for (int Diff_counter = 0; Diff_counter < 3; Diff_counter++) {//Loops through the various drug binding rates
            Diffusion_drug = Diffusion_vec[Diff_counter];// Setting Diffusion for drug binding.
            
            switch (Diff_counter) {
                case 0:{
                    Diffval = "1e5";
                    break;
                }
                case 1:{
                    Diffval = "1e3";
                    break;
                }
                case 2:{
                    Diffval = "1e1";
                    break;
                }
            }//recording the diffusion value for the file name
            
            if (Bind_Scheme == 0) {
                result_file.open(Model_type+'_'+bindscheme+'_'+"drug_Upstroke_rest_30s_stimdur06.txt");
                
                Diff_counter = 4;
            } else {
                result_file.open(Model_type+'_'+bindscheme+'_'+"drug_Upstroke_rest_30s_Diffval"+Diffval+"_kd10e-6_Drug20e-6_stimdur06.txt");
            }
            
            
            Cell.Cell_type = 3;        // Cell type = (1) for endo; (2) for M,  (3) epi
            
            t = 0.0;
            tTemp = 0.0; //initializing t and tTemp
            
            load_IC(&Cell);
            
            /*********************************************************************************************************************************
             Begin Time Loop Here
             *********************************************************************************************************************************/
            
            for (S1_cycle = 0; S1_cycle <=3; S1_cycle = S1_cycle + 1){                //Enter the number of beats that you want to simulate
                
                if (S1_cycle ==0){interval=waitTime;}//allows the model to equilibrate for the waitTime set in Global_variables.h
                
                else {interval=rest;}
                //sets the BCL
                
                Calculate_Reset (&Cell, t, tTemp);//Resets cell specific parameters every beat such as APD, DI, V90 etc.
                for (tTemp = 0; tTemp <= interval; tTemp=tTemp+dt){
                    Calculate_Points (&Cell, t, tTemp); //Calculates t_min, V_min, t_max, V_max etc. etc.; Update calculates axial currents and updates the voltage
                    
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
                    
                    
                    if (((interval - tTemp)<dt)&&(S1_cycle > 1)) {
                        result_file << interval << ", " << S1_cycle << ", " << Cell.peak_slope << ", " << Cell.APD_90 << std::endl;//recording the BCL, number of stimuli, peak Upstroke, and APD90 for each stimuli
                        
                    }
                    
                    
                    /*******************************************************************************************************************************************************
                     End of result files
                     *******************************************************************************************************************************************************/
                    
                    counter++;
                    
                    if (fmod(t,5000)<dt) {std::cout << "S1_cycle=" <<S1_cycle<< ", " << "t=" << t << ", " "interval =" << interval << ", " << "Drug = "<< Drug << ", " << "Mutant = "<< mutant << std::endl;}
                    
                }//End of For tTemp ==
                
            }//End of for S1
            
            
            
            result_file.close();
            
        }
    }
}

#endif /* Upstroke_rest_30s_all_models_mult_time_const_h */
