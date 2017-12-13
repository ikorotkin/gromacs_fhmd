#ifndef FHMD_PARAMS_H_
#define FHMD_PARAMS_H_

#define FHMD_VERSION        3.00            /* FHMD model version */

#define FHMD_MAX_LENGTH     1000.0          /* Maximum length scale [nm] -- for control purpose only */

#define FHMD_kB             0.00831451      /* Boltzmann constant [kJ/(mol*K)] */

#define FHMD_DEBUG                          /* Write debug information */
//#define FHMD_DEBUG_GRID                     /* Print FH grid coordinates to the screen */
//#define FHMD_DEBUG_FH                       /* Print FH debug information */
//#define FHMD_DEBUG_COM                      /* Print protein COM coordinates */
//#define FHMD_DEBUG_INTERPOL                 /* Print interpolated values for debugging */

enum FHMD_EOS {eos_argon, eos_spce};        /* Equation of state enumeration */

/* EOS coefficients for RIGID SPC/E water */
#define FHMD_EOS_SPCE_MU    409.496
#define FHMD_EOS_SPCE_KAPPA 933.41
#define FHMD_EOS_SPCE_A     0.010102137
#define FHMD_EOS_SPCE_B    -10.133069
#define FHMD_EOS_SPCE_C     2428.9203

/* Color codes */
#ifndef _MSC_VER
    #define RESET_COLOR     "\e[m"
    #define MAKE_RED        "\e[91m"
    #define MAKE_GREEN      "\e[92m"
    #define MAKE_YELLOW     "\e[93m"
    #define MAKE_BLUE       "\e[94m"
    #define MAKE_PURPLE     "\e[95m"
    #define MAKE_LIGHT_BLUE "\e[96m"
#else
    #define RESET_COLOR
    #define MAKE_RED
    #define MAKE_GREEN
    #define MAKE_YELLOW
    #define MAKE_BLUE
    #define MAKE_PURPLE
    #define MAKE_LIGHT_BLUE
#endif

#endif /* FHMD_PARAMS_H_ */
