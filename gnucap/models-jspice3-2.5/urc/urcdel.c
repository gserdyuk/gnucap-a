/***************************************************************************
JSPICE3 adaptation of Spice3f2 - Copyright (c) Stephen R. Whiteley 1992
Copyright 1990 Regents of the University of California.  All rights reserved.
Authors: 1985 Thomas L. Quarles
         1993 Stephen R. Whiteley
****************************************************************************/

#include "spice.h"
#include <stdio.h>
#include "urcdefs.h"
#include "sperror.h"
#include "util.h"


int
URCdelete(inModel,name,inst)
    GENmodel *inModel;
    IFuid name;
    GENinstance **inst;
{
    URCmodel *model = (URCmodel *)inModel;
    URCinstance **fast = (URCinstance**)inst;
    URCinstance **prev = NULL;
    URCinstance *here;

    for( ; model ; model = model->URCnextModel) {
        prev = &(model->URCinstances);
        for(here = *prev; here ; here = *prev) {
            if(here->URCname == name || (fast && here==*fast) ) {
                *prev= here->URCnextInstance;
                FREE(here);
                return(OK);
            }
            prev = &(here->URCnextInstance);
        }
    }
    return(E_NODEV);
}
