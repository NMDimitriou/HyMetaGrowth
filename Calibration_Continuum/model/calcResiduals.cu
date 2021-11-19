#include "declarations.cuh"

/* Calculate residuals between simulation and experimental data
 * Nikoalos Dimitriou, McGill, 2021
*/

void calcResiduals(int c, double se, InitialData& id, DataArray& host, ExpData& expdat,FILE* fnrmse)
{
#if (STUDY == CALIBRATION)
	#ifdef NRMSE
	double seF=0.0;
	double rmse, sd, nrmse; //RMSE, SD or Max-Min, NRMSE
	#endif
	int k;

	switch(c)
   	{
       	case 0 :
        for (k=0; k<total; k++)
		{
            se = host.a_intm[k]-expdat.dtwo[k];
            id.SSE += se*se;
			#ifdef NRMSE
            seF += se*se;
			#endif
        }
        printf("Day 2: SSE = %.15f\n", id.SSE);
        #ifdef NRMSE
		//save_nrmse(host);
        rmse  = sqrt(seF/total);
        sd    = calc_ExpSd(expdat.dtwo);
        nrmse = rmse/sd;
        printf("RMSE  = %.15f\n",  rmse);
        printf("NRMSE = %.15f\n", nrmse);
        fprintf(fnrmse,"%.15f\n", nrmse); fflush(fnrmse);
        seF=0.0;
        #endif
        break;

        case 1 :
        for (k=0; k<total; k++)
		{
            se = (host.a_intm[k]-expdat.dfive[k]);
            id.SSE += se*se;
			#ifdef NRMSE
            seF += se*se;
			#endif
        }
        printf("Day 5: SSE = %.15f\n", id.SSE);
        #ifdef NRMSE
        rmse  = sqrt(seF/total);
        sd    = calc_ExpSd(expdat.dfive);
        nrmse = rmse/sd;
        printf("RMSE  = %.15f\n",  rmse);
        printf("NRMSE = %.15f\n", nrmse);
        fprintf(fnrmse,"%.15f\n", nrmse); fflush(fnrmse);
        seF=0.0;
        #endif
        break;

        case 2 :
        for (k=0; k<total; k++)
		{
         	se= (host.a_intm[k]-expdat.dsvn[k]);
            id.SSE += se*se;
			#ifdef NRMSE
            seF += se*se;
			#endif
        }
		printf("Day 7: SSE = %.15f\n", id.SSE);
        #ifdef NRMSE
        rmse  = sqrt(seF/total);
        sd    = calc_ExpSd(expdat.dsvn);
        nrmse = rmse/sd;
        printf("RMSE  = %.15f\n",  rmse);
        printf("NRMSE = %.15f\n", nrmse);
        fprintf(fnrmse,"%.15f\n", nrmse); fflush(fnrmse);
        seF=0.0;
        #endif
        break;

        case 3 :
        for (k=0; k<total; k++)
		{
           	se= (host.a_intm[k]-expdat.dnine[k]);
            id.SSE += se*se;
			#ifdef NRMSE
            seF += se*se;
			#endif
        }
        printf("Day 9: SSE = %.15f\n", id.SSE);
        #ifdef NRMSE
        rmse  = sqrt(seF/total);
        sd    = calc_ExpSd(expdat.dnine);
        nrmse = rmse/sd;
        printf("RMSE  = %.15f\n",  rmse);
        printf("NRMSE = %.15f\n", nrmse);
        fprintf(fnrmse,"%.15f\n", nrmse); fflush(fnrmse);
        seF=0.0;
        #endif
        break;

        case 4 :
        for (k=0; k<total; k++)
		{
           	se= (host.a_intm[k]-expdat.dtwlv[k]);
            id.SSE += se*se;
			#ifdef NRMSE
            seF += se*se;
			#endif
        }
        printf("Day 12: SSE = %.15f\n", id.SSE);
        #ifdef NRMSE
        rmse  = sqrt(seF/total);
        sd    = calc_ExpSd(expdat.dtwlv);
        nrmse = rmse/sd;
        printf("RMSE  = %.15f\n",  rmse);
        printf("NRMSE = %.15f\n", nrmse);
        fprintf(fnrmse,"%.15f\n", nrmse); fflush(fnrmse);
        seF=0.0;
        #endif		
		break;

        case 5 :
        for (k=0; k<total; k++)
		{
           	se= (host.a_intm[k]-expdat.dfrtn[k]);
            id.SSE += se*se;
			#ifdef NRMSE
            seF += se*se;
			#endif
        }
        printf("Day 14: SSE = %.15f\n", id.SSE);
        #ifdef NRMSE
        rmse  = sqrt(seF/total);
        sd    = calc_ExpSd(expdat.dfrtn);
        nrmse = rmse/sd;
        printf("RMSE  = %.15f\n",  rmse);
        printf("NRMSE = %.15f\n", nrmse);
        fprintf(fnrmse,"%.15f\n", nrmse); fflush(fnrmse);
        seF=0.0;
		fclose(fnrmse);
        #endif
        break;
	}
#endif
}
