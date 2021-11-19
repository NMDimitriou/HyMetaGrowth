#include "declarations.cuh"

/* Just ignore this for now
*/

FILE* save_nrmse(DataArray& host)
{
	FILE* fnrmse;
    #if   ( MODEL == KSCMD  )
    char ferr[4096]; sprintf(ferr, "../NRMSE_KSCMD/nrmse_%s.txt", host.dat); 
//	FILE* fnrmse; 
	fnrmse = fopen(ferr, "wa");
    #elif ( MODEL == KSMD  )
    char ferr[4096]; sprintf(ferr, "../NRMSE_KSMD/nrmse_%s.txt", host.dat); 
	//FILE* fnrmse; 
	fnrmse = fopen(ferr, "wa");
    #elif ( MODEL == KSC )
    char ferr[4096]; sprintf(ferr, "../NRMSE_KSC/nrmse_%s.txt", host.dat); 
	//FILE* fnrmse; 
	fnrmse = fopen(ferr, "wa");
    #elif ( MODEL == FK )
    char ferr[4096]; sprintf(ferr, "../NRMSE_FK/nrmse_%s.txt", host.dat); 
	//FILE* fnrmse; 
	fnrmse = fopen(ferr, "wa");
    #endif

	return fnrmse;
}
