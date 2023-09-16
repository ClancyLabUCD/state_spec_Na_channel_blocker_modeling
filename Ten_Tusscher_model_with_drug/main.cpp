//
//  main.cpp
//  Ten_Tusscher_model_with_drug
//
//  Created by Steffen Docken on 12/5/17.
//  Copyright © 2017 Steffen Docken (Lewis Lab). All rights reserved.
//  Published in the Journal of Theoretical Biology under the title "Rate-dependent effects of state-specific sodium channel blockers in cardiac tissue: Insights from idealized models"
//

#include <iostream>
#include <fstream>
#include <math.h>
#include "Upstroke_vs_period_all_models_mult_time_const.h"
#include "Upstroke_all_models_mult_time_const_1freq.h"
#include "Upstroke_rest_30s_all_models_mult_time_const.h"
#include "APDs_all_HH_form.h"

int main() {
    Upstroke_vs_period_all_models_mult_time_const();
    //Upstroke_all_models_mult_time_const_1freq();
    //Upstroke_rest_30s_all_models_mult_time_const();
    //APDs_all_HH_form();
    
    return 0;
}
