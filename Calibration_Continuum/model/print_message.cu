#include "declarations.cuh"


void print_message(InitialData& id, double sigma_data2)
{

	#if   ( MODEL == KSCMD )
        printf("\n##################################################################### \n \n");
        printf("Initializing 3D-KS model with Chemotaxis and ECM degradation \n(0,2.5)x(0,2.5)x(0,0.9) mm \n");
        printf("Total time: 14 days with %d implicit time steps\n", id.itmax);
        printf("dx,y,z = %f\n", id.dx);
	#if   ( IC == EXP)
	printf("Initial signal density has exponential gradient from bottom to top.\n");
	#elif ( IC == LIN)
	printf("Initial signal density has linear gradient from bottom to top.\n");
	#elif ( IC == ID)
	printf("Initial signal density is identical to cell density.\n");
	#endif
        printf("Da = %lf, s = %lf, Chi = %lf, Chi_ecm = %lf\n", id.da, id.s, id.chi, id.chi_ecm);
        printf("Db = %lf, r = %lf\n", id.db, id.r);
        printf("Dc = %lf\n", id.dc);
        printf("sigma2 = %.15f\n", sigma_data2);
        printf("\n######################################################################## \n");
        printf(" \n");
        #elif   ( MODEL == KSMD )
        printf("\n##################################################################### \n \n");
        printf("Initializing 3D-KS model with ECM degradation \n(0,2.5)x(0,2.5)x(0,0.9) mm \n");
        printf("Total time: 14 days with %d implicit time steps\n", id.itmax);
        printf("dx,y,z = %f\n", id.dx);
        printf("Da = %lf, s = %lf, Chi_ecm = %lf\n", id.da, id.s, id.chi_ecm);
        printf("Dc = %lf\n", id.dc);
        printf("sigma2 = %.15f\n", sigma_data2);
        printf("\n######################################################################## \n");
        printf(" \n");
        #elif ( MODEL == KSC || MODEL == TEST)
        printf("\n########################################## \n \n");
        printf("Initializing 3D-KS model with Chemotaxis \n(0,2.5)x(0,2.5)x(0,0.9) mm \n");
        printf("Total time: 14 days with %d implicit time steps\n", id.itmax);
        printf("dx,y,z = %f\n", id.dx);
	#if   ( IC == EXP)
        printf("Initial signal density has exponential gradient from bottom to top.\n");
        #elif ( IC == LIN)
        printf("Initial signal density has linear gradient from bottom to top.\n");
        #elif ( IC == ID)
        printf("Initial signal density is identical to cell density.\n");
        #endif
        printf("Da = %lf, s = %lf, Chi = %lf\n", id.da, id.s, id.chi);
        printf("Db = %lf, r = %lf\n", id.db, id.r);
        printf("sigma2 = %.15f\n", sigma_data2);
        printf("\n########################################## \n");
        printf(" \n");
        #elif ( MODEL == FK)
        printf("\n########################################## \n \n");
        printf("Initializing 3D-FK \n(0,2.5)x(0,2.5)x(0,0.9) mm \n");
        printf("Total time: 14 days with %d implicit time steps\n", id.itmax);
        printf("dx,y,z = %f\n", id.dx);
        printf("Da = %lf, s = %lf\n", id.da, id.s);
        printf("sigma2 = %.15f\n", sigma_data2);
        printf("\n########################################## \n");
        printf(" \n");
        #endif


}
