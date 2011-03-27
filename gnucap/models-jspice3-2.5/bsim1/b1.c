/***************************************************************************
JSPICE3 adaptation of Spice3f2 - Copyright (c) Stephen R. Whiteley 1992
Copyright 1990 Regents of the University of California.  All rights reserved.
Authors: 1985 Hong J. Park, Thomas L. Quarles 
         1993 Stephen R. Whiteley
****************************************************************************/

#include "spice.h"
#include <stdio.h>
#include "bsim1def.h"
#include "uflags.h"

static IFparm B1pTable[] = { /* parameters */
 IOP( "l",   BSIM1_L,      IF_REAL,    "Length"),
 IOP( "w",   BSIM1_W,      IF_REAL,    "Width"),
 IOP( "ad",  BSIM1_AD,     IF_REAL,    "Drain area"),
 IOP( "as",  BSIM1_AS,     IF_REAL,    "Source area"),
 IOP( "pd",  BSIM1_PD,     IF_REAL,    "Drain perimeter"),
 IOP( "ps",  BSIM1_PS,     IF_REAL,    "Source perimeter"),
 IOP( "nrd", BSIM1_NRD,    IF_REAL,    "Number of squares in drain"),
 IOP( "nrs", BSIM1_NRS,    IF_REAL,    "Number of squares in source"),
 IOP( "off", BSIM1_OFF,    IF_FLAG,    "Device is initially off"),
 IOP( "vds", BSIM1_IC_VDS, IF_REAL,    "Initial D-S voltage"),
 IOP( "vgs", BSIM1_IC_VGS, IF_REAL,    "Initial G-S voltage"),
 IOP( "vbs", BSIM1_IC_VBS, IF_REAL,    "Initial B-S voltage"),
 IP(  "ic",  BSIM1_IC,     IF_VECTOR , "Vector of DS,GS,BS initial voltages")
};

static IFparm B1mPTable[] = { /* model parameters */
 IOP( "vfb",   BSIM1_MOD_VFB0,      IF_REAL,"Flat band voltage"),
 IOP( "lvfb",  BSIM1_MOD_VFBL,      IF_REAL,
                    "Length dependence of vfb"),
 IOP( "wvfb",  BSIM1_MOD_VFBW,      IF_REAL,
                    "Width dependence of vfb"),
 IOP( "phi",   BSIM1_MOD_PHI0,      IF_REAL,
        "Strong inversion surface potential "),
 IOP( "lphi",  BSIM1_MOD_PHIL,      IF_REAL,
        "Length dependence of phi"),
 IOP( "wphi",  BSIM1_MOD_PHIW,      IF_REAL,
        "Width dependence of phi"),
 IOP( "k1",    BSIM1_MOD_K10,       IF_REAL,
        "Bulk effect coefficient 1"),
 IOP( "lk1",   BSIM1_MOD_K1L,       IF_REAL,
        "Length dependence of k1"),
 IOP( "wk1",   BSIM1_MOD_K1W,       IF_REAL,
        "Width dependence of k1"),
 IOP( "k2",    BSIM1_MOD_K20,       IF_REAL,
        "Bulk effect coefficient 2"),
 IOP( "lk2",   BSIM1_MOD_K2L,       IF_REAL,
        "Length dependence of k2"),
 IOP( "wk2",   BSIM1_MOD_K2W,       IF_REAL,
        "Width dependence of k2"),
 IOP( "eta",   BSIM1_MOD_ETA0,      IF_REAL,
        "VDS dependence of threshold voltage"),
 IOP( "leta",  BSIM1_MOD_ETAL,      IF_REAL,
        "Length dependence of eta"),
 IOP( "weta",  BSIM1_MOD_ETAW,      IF_REAL,
        "Width dependence of eta"),
 IOP( "x2e",   BSIM1_MOD_ETAB0,     IF_REAL,
        "VBS dependence of eta"),
 IOP( "lx2e",  BSIM1_MOD_ETABL,     IF_REAL,
        "Length dependence of x2e"),
 IOP( "wx2e",  BSIM1_MOD_ETABW,     IF_REAL,
        "Width dependence of x2e"),
 IOP( "x3e",   BSIM1_MOD_ETAD0,     IF_REAL,
        "VDS dependence of eta"),
 IOP( "lx3e",  BSIM1_MOD_ETADL,     IF_REAL,
        "Length dependence of x3e"),
 IOP( "wx3e",  BSIM1_MOD_ETADW,     IF_REAL,
        "Width dependence of x3e"),
 IOP( "dl",    BSIM1_MOD_DELTAL,    IF_REAL,
        "Channel length reduction in um"),
 IOP( "dw",    BSIM1_MOD_DELTAW,    IF_REAL,
        "Channel width reduction in um"),
 IOP( "muz",   BSIM1_MOD_MOBZERO,   IF_REAL,
        "Zero field mobility at VDS=0 VGS=VTH"),
 IOP( "x2mz",  BSIM1_MOD_MOBZEROB0, IF_REAL,
        "VBS dependence of muz"),
 IOP( "lx2mz", BSIM1_MOD_MOBZEROBL, IF_REAL,
        "Length dependence of x2mz"),
 IOP( "wx2mz", BSIM1_MOD_MOBZEROBW, IF_REAL,
        "Width dependence of x2mz"),
 IOP( "mus",   BSIM1_MOD_MOBVDD0,   IF_REAL,
        "Mobility at VDS=VDD VGS=VTH, channel length modulation"),
 IOP( "lmus",  BSIM1_MOD_MOBVDDL,   IF_REAL,
        "Length dependence of mus"),
 IOP( "wmus",  BSIM1_MOD_MOBVDDW,   IF_REAL,
        "Width dependence of mus"),
 IOP( "x2ms",  BSIM1_MOD_MOBVDDB0,  IF_REAL,
            "VBS dependence of mus"),
 IOP( "lx2ms", BSIM1_MOD_MOBVDDBL,  IF_REAL,
        "Length dependence of x2ms"),
 IOP( "wx2ms", BSIM1_MOD_MOBVDDBW,  IF_REAL,
        "Width dependence of x2ms"),
 IOP( "x3ms",  BSIM1_MOD_MOBVDDD0,  IF_REAL,
        "VDS dependence of mus"),
 IOP( "lx3ms", BSIM1_MOD_MOBVDDDL,  IF_REAL,
        "Length dependence of x3ms"),
 IOP( "wx3ms", BSIM1_MOD_MOBVDDDW,  IF_REAL,
        "Width dependence of x3ms"),
 IOP( "u0",    BSIM1_MOD_UGS0,      IF_REAL,
        "VGS dependence of mobility"),
 IOP( "lu0",   BSIM1_MOD_UGSL,      IF_REAL,
                "Length dependence of u0"),
 IOP( "wu0",   BSIM1_MOD_UGSW,      IF_REAL,
        "Width dependence of u0"),
 IOP( "x2u0",  BSIM1_MOD_UGSB0,     IF_REAL,
        "VBS dependence of u0"),
 IOP( "lx2u0", BSIM1_MOD_UGSBL,     IF_REAL,
        "Length dependence of x2u0"),
 IOP( "wx2u0", BSIM1_MOD_UGSBW,     IF_REAL,
        "Width dependence of x2u0"),
 IOP( "u1",    BSIM1_MOD_UDS0,      IF_REAL,
        "VDS depence of mobility, velocity saturation"),
 IOP( "lu1",   BSIM1_MOD_UDSL,      IF_REAL,
        "Length dependence of u1"),
 IOP( "wu1",   BSIM1_MOD_UDSW,      IF_REAL,
        "Width dependence of u1"),
 IOP( "x2u1",  BSIM1_MOD_UDSB0,     IF_REAL,
        "VBS depence of u1"),
 IOP( "lx2u1", BSIM1_MOD_UDSBL,     IF_REAL,
        "Length depence of x2u1"),
 IOP( "wx2u1", BSIM1_MOD_UDSBW,     IF_REAL,
        "Width depence of x2u1"),
 IOP( "x3u1",  BSIM1_MOD_UDSD0,     IF_REAL,
        "VDS depence of u1"),
 IOP( "lx3u1", BSIM1_MOD_UDSDL,     IF_REAL,
        "Length dependence of x3u1"),
 IOP( "wx3u1", BSIM1_MOD_UDSDW,     IF_REAL,
        "Width depence of x3u1"),
 IOP( "n0",    BSIM1_MOD_N00,       IF_REAL,
            "Subthreshold slope"),
 IOP( "ln0",   BSIM1_MOD_N0L,       IF_REAL,
        "Length dependence of n0"),
 IOP( "wn0",   BSIM1_MOD_N0W,       IF_REAL,
        "Width dependence of n0"),
 IOP( "nb",    BSIM1_MOD_NB0,       IF_REAL,
        "VBS dependence of subthreshold slope"),
 IOP( "lnb",   BSIM1_MOD_NBL,       IF_REAL,
        "Length dependence of nb"),
 IOP( "wnb",   BSIM1_MOD_NBW,       IF_REAL,
        "Width dependence of nb"),
 IOP( "nd",    BSIM1_MOD_ND0,       IF_REAL,
        "VDS dependence of subthreshold slope"),
 IOP( "lnd",   BSIM1_MOD_NDL,       IF_REAL,
        "Length dependence of nd"),
 IOP( "wnd",   BSIM1_MOD_NDW,       IF_REAL,
        "Width dependence of nd"),
 IOP( "tox",   BSIM1_MOD_TOX,       IF_REAL,
        "Gate oxide thickness in um"),
 IOP( "temp",  BSIM1_MOD_TEMP,      IF_REAL,
        "Temperature in degree Celcius"),
 IOP( "vdd",   BSIM1_MOD_VDD,       IF_REAL,
                "Supply voltage to specify mus"),
 IOPA("cgso",  BSIM1_MOD_CGSO,      IF_REAL,
          "Gate source overlap capacitance per unit channel width(m)"),
 IOPA("cgdo",  BSIM1_MOD_CGDO,      IF_REAL,
          "Gate drain overlap capacitance per unit channel width(m)"),
 IOPA("cgbo",  BSIM1_MOD_CGBO,      IF_REAL,
          "Gate bulk overlap capacitance per unit channel length(m)"),
 IOP( "xpart", BSIM1_MOD_XPART,     IF_REAL,
      "Flag for channel charge partitioning"),
 IOP( "rsh",   BSIM1_MOD_RSH,       IF_REAL,
      "Source drain diffusion sheet resistance in ohm per square"),
 IOP( "js",    BSIM1_MOD_JS,        IF_REAL,
      "Source drain junction saturation current per unit area"),
 IOP( "pb",    BSIM1_MOD_PB,        IF_REAL,
      "Source drain junction built in potential"),
 IOPA("mj",    BSIM1_MOD_MJ,        IF_REAL,
       "Source drain bottom junction capacitance grading coefficient"),
 IOPA("pbsw",  BSIM1_MOD_PBSW,      IF_REAL,
       "Source drain side junction capacitance built in potential"),
 IOPA("mjsw",  BSIM1_MOD_MJSW,      IF_REAL,
       "Source drain side junction capacitance grading coefficient"),
 IOPA("cj",    BSIM1_MOD_CJ,        IF_REAL,
       "Source drain bottom junction capacitance per unit area"),
 IOPA("cjsw",  BSIM1_MOD_CJSW,      IF_REAL,
       "Source drain side junction capacitance per unit area"),
 IOP( "wdf",   BSIM1_MOD_DEFWIDTH,  IF_REAL,
       "Default width of source drain diffusion in um"),
 IOP( "dell",  BSIM1_MOD_DELLENGTH, IF_REAL,
       "Length reduction of source drain diffusion"),
 IP(  "nmos",  BSIM1_MOD_NMOS,      IF_FLAG,
       "Flag to indicate NMOS"),
 IP(  "pmos",  BSIM1_MOD_PMOS,      IF_FLAG,
       "Flag to indicate PMOS"),
};

static char *B1names[] = {
    "Drain",
    "Gate",
    "Source",
    "Bulk"
};

static char *B1modNames[] = {
    "nmos",
    "pmos",
    NULL
};

static IFkeys B1keys[] = {
    { 'm', NUMELEMS(B1names), B1names, 0, 0 },
};


static int B1kSize = NUMELEMS(B1keys);
static int B1pTSize = NUMELEMS(B1pTable);
static int B1mPTSize = NUMELEMS(B1mPTable);
static int B1iSize = sizeof(B1instance);
static int B1mSize = sizeof(B1model);


SPICEdev B1info = {
    {   "B1",
        "Berkeley Short Channel IGFET Model",

        &B1kSize,
        B1keys,
        8,      /* level 4 */
        B1modNames,
        GENmosParse,

        &B1pTSize,
        B1pTable,

        &B1mPTSize,
        B1mPTable,
    },

    B1param,
    B1mParam,
    B1load,
    B1setup,
    B1setup,
    B1temp,
    B1trunc,
    NULL,
    B1acLoad,
    NULL,
    GENdestroy,
    GENmDelete,
    GENdelete, 
    B1getic,
    B1ask,
    B1mAsk,
    B1pzLoad,
    B1convTest,
    B1disto,
    NULL, /* NOISE */

    &B1iSize,
    &B1mSize
};
