extern "C" {
  #include "b3soidef.h"
  #define DEV_b3soi
  #include "b3soiitf.h"
}

#define info	 B3SOIinfo
#define INSTANCE B3SOIinstance
#define MODEL	 B3SOImodel

#define SPICE_LETTER "M"
#define DEVICE_TYPE "bsimsoi3|bsimsoi3p11"
#define MIN_NET_NODES 4
#define MAX_NET_NODES 7
#define INTERNAL_NODES 23
#define MODEL_TYPE "nmos9|pmos9|nmos56|pmos56"
<<<<<<< HEAD:gnucap/models-bsim/BSIMSOI3p11/wrapper.h
#define UNCONNECTED_NODES FLOAT
=======
#define UNCONNECTED_NODES uFLOAT
>>>>>>> gnucap-main/master:gnucap/models-bsim/BSIMSOI3p11/wrapper.h

static std::string port_names[] = {"d", "g", "s", "bg", "p", "body", "temp"};
static std::string state_names[] = {};
