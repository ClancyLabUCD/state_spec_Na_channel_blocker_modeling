//
//  main.cpp
//  1D_modified_Ten_Tusscher_with_drug_2017_Summer
//
//  Created by Steffen Docken on 7/11/17.
//  Copyright © 2017 Steffen Docken (Lewis Lab). All rights reserved.
//  Published in the Journal of Theoretical Biology under the title "Rate-dependent effects of state-specific sodium channel blockers in cardiac tissue: Insights from idealized models"
//

#include <iostream>
#include <fstream>
#include <math.h>

#include "VW_initialization_cable.h"
#include "VW_find_VW_high_resolution_cable.h"
#include "VW_find_VW_GRI_high_resolution_cable.h"
#include "CV_vs_BCL_cable.h"

int main() {
    //VW_initialization_cable();
    //VW_find_VW_high_resolution_cable();
    //VW_find_VW_GRI_high_resolution_cable();
    CV_vs_BCL_cable();
    return 0;
}
