//
//  Na_currents.h
//  1D_modified_Ten_Tusscher_with_drug_2017_Summer
//
//  Created by Steffen Docken on 7/11/17.
//  Copyright © 2017 Steffen Docken (Lewis Lab). All rights reserved.
//  Published in the Journal of Theoretical Biology under the title "Rate-dependent effects of state-specific sodium channel blockers in cardiac tissue: Insights from idealized models"
//

////7-11-17 This code was copied from Ten_Tusscher_model_with_drug_2017_Winter

////The Ten Tusscher model code was copied from code given to me by Pei-Chi Yang

#ifndef Na_currents_h
#define Na_currents_h

#include <math.h>
#include "Global_variables.h"
#include "I_Na_TT.h"
#include "I_Na_Simple.h"


//This function uses the correct I_Na model depending on the value of Na_model
void Calculate_I_Na(Cell_param *Cell_ptr, double t, double kon){
    if (Na_model == 0) {
        Calculate_I_Na_TT(Cell_ptr, t, kon);
    } else if (Na_model == 1){
        Calculate_I_Na_Simple(Cell_ptr, t, kon);
    }
}

#endif /* Na_currents_h */
