#ifndef FHMD_PARAMS_H_
#define FHMD_PARAMS_H_

#define fhmd_version        3.00            /* FHMD model version */

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
