#include "declarations.cuh"

/* Calculates response of simulation output to parameters
   (used in sensitivity analysis)
   Nikolaos Dimitriou, McGill, 2021
*/

void calcResponse(InitialData& id, DataArray& host)
{
#if ( STUDY == SA)
	int k;

	for (k=0; k<total; k++)
	{
		if(host.a_intm[k] > 1e-03)
		{
			id.tot_num ++;
			if(k < zsize) id.bot_num ++;
		}
		//id.tot_num += host.a_intm[k]; 
		//if(k < zsize) id.bot_num += host.a_intm[k];

	}
//	id.tot_num *= total;
//	id.bot_num *= zsize;
#endif
}
