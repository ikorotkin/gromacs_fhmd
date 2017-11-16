#include "data_structures.h"


int fhmd_init(t_commrec *cr)
{

    if(MASTER(cr)) {

        /* Initial output */
        printf(MAKE_GREEN "\n  Aston University, Department of Mathematics, Dmitry Nerukh Research Group\n");
        printf("  Queen Mary University of London, School of Engineering and Material Science\n\n");
        printf(MAKE_BLUE "     Hybrid Molecular Dynamics - 2-Way Coupling Parallel VERSION %4.2f\n\n", fhmd_version);

        /* Read FHMD parameters */
        printf(MAKE_GREEN "FHMD: Reading parameters from coupling.prm..." RESET_COLOR " ");

        /* Print parameters to screen and output file */

        /* Allocate memory */


        fflush(stdout);

    }

    /* Broadcast parameters to all threads */

    return 1;   // Success

}
