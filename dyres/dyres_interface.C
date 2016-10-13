#include "dyres_interface.h"
#include "settings.h"
#include "interface.h"

void dyres::init()
{
    //Initialise some DYRES settings
    g_param_.g_param_ = opts.g_param;
    nnlo_.order_ = opts.order;            //order (0=LO, 1=NLO, 2=NNLO)
    opts_.fixedorder_  = opts.fixedorder; //fixed order/resummation switch
    qtsub_.xqtcut_= opts.xqtcut;          //Cut on qt/Q
    qtsub_.qtcut_= opts.qtcut;            //Cut on qt
    //move here the flaq.eq.0 initialisation part of resumm() in main2 instead of using this initialisation flag
}
