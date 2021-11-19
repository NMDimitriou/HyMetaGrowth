#include "declarations.cuh"

/* Initialize simulation parameters and process
 * Nikolaos Dimitriou, McGill, 2021
*/

void simulate(DataArray& host, DataArray& device1, DataArray& device2, CAArray& casp, vector<Cell>& cells, vector<Cell>& cellsTmp, InitialData& id)
{
	//define parameters
	int c  = 0 ; //checks for export time
  	int ct = 0 ; //counter of function calls
	double step_dt = 0.; //counter for explicit time-step

  	id.dx        = 2.5f/(ysize -1);	
	id.cfl       = 0.15;
	
	#if ( MODEL == FK || MODEL == KSC || MODEL == KSMD || MODEL == KSCMD || MODEL == TEST)
	double dt_a  = id.cfl*id.dx*id.dx/id.da;
	id.dt_imp    = dt_a;
    #endif
    #if ( MODEL == KSC || MODEL == KSCMD || MODEL == TEST)
	double dt_b  = id.cfl*id.dx*id.dx/id.db;
	id.dt_imp    = min(dt_a,dt_b);
    #endif
    #if ( MODEL == KSMD || MODEL == KSCMD )
    double dt_c  = id.cfl*id.dx*id.dx/id.dc;
    id.dt_imp    = min(id.dt_imp,dt_c);
    #endif

	id.dt_imp    = min(id.dt_imp,0.1);
 	id.itmax     = 14.0/id.dt_imp;

	#if ( MODEL == FK || MODEL == KSC || MODEL == KSMD || MODEL == KSCMD || MODEL == TEST)
    double cfl_a = id.da*id.dt_imp/(id.dx*id.dx);
    #endif
    #if ( MODEL == KSC || MODEL == KSCMD || MODEL == TEST)
    double cfl_b = id.db*id.dt_imp/(id.dx*id.dx);
    double chi_b_dx   = id.chi/id.dx;
    #endif
    #if   ( MODEL == KSMD || MODEL == KSCMD )
    double cfl_c = id.dc*id.dt_imp/(id.dx*id.dx);
    double chi_ecm_dx = id.chi_ecm/id.dx;
    #endif

	casp.adhes      = -5; //n=1 at day6
 	casp.delaydiv   =  0.0;
	casp.phenchange =  0; //2;
	casp.divpot     =  100; //

	printf("\n########################################## \n \n");
	printf("Initializing 3D-KSC model \n(0,2.5)x(0,2.5)x(0,0.9) mm \n");
	printf("Total time: 14 days with %d implicit time steps\n", id.itmax);
	printf("dx,y,z = %f\n", id.dx);
	printf("Da = %lf, s = %lf, Chi = %lf\n", id.da, id.s, id.chi);
	printf("Db = %lf, r = %lf\n", id.db, id.r);
	printf("\n########################################## \n");
	printf("\n");
	
	//INITIALIZE
	//-host
	hostInitialize(host);	
	if(*host.dt >= id.dt_imp/4.0) *host.dt = id.dt_imp/4.0;
	
	//-CA
	//scaleGrowth(id);
	casp.agediv   = log(2.0)/(id.s*id.dt_imp);
	casp.spdeath  = 0.005*id.dt_imp; //
 	initialize_ca(host, casp, cells);
 	printf("Age div : %d\n", casp.agediv);
	char fname[4096];
	char fcount[4096];
	sprintf(fname, "ca_res_%s.csv", host.dat);	
	sprintf(fcount, "num_cells_%s.txt", host.dat);
	FILE *fca;
	FILE *fna;
	fca = fopen(fname, "wa");
	fna = fopen(fcount, "wa");
 	if (fca==0) { printf("  error creating saving file\n"); }
	if (fna==0) { printf("  error creating saving file\n"); }

	//-device
	gpuInitialize(id, host, device1, device2);

	int exportInterval[14] = {int(1.0f/id.dt_imp) ,int(2.0f/id.dt_imp) ,int(3.0f/id.dt_imp) ,int(4.0f/id.dt_imp),
                              int(5.0f/id.dt_imp) ,int(6.0f/id.dt_imp) ,int(7.0f/id.dt_imp) ,int(8.0f/id.dt_imp),
                              int(9.0f/id.dt_imp) ,int(10.0f/id.dt_imp),int(11.0f/id.dt_imp),int(12.0f/id.dt_imp),
                              int(13.0f/id.dt_imp),int(14.0f/id.dt_imp)};

	// Export Day 0.
	ExportCA(fca,id,cells,0);
	ExportCellCount(fca,id,cells,0);
	ExportCheckpoint("a_double","b_double","ecm_double",host,device1,0,1);

	//events for time measurement
	float run_time=0.0f;
  	cudaEvent_t start, stop;
  	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	//COMPUTE
  	printf("Starting simulation...\n");
  	printf("Implicit time-step = %lf\n",id.dt_imp);
  	printf("Explicit time-step = %lf\n",*host.dt);
  	cudaEventRecord(start, 0);

	for(int i=1; i<=id.itmax; i++)
	{
    	//Strang splitting
		//EXPLICIT STEP
        #if ( MODEL == FK || MODEL == KSCMD || MODEL == KSMD || MODEL == KSC || MODEL == TEST)
		while(step_dt - 0.5*id.dt_imp < -1.0e-12)
		{
			ct ++;
			#if   ( MODEL == KSCMD )
            gpuStep_explicit(device1, device2, id, chi_b_dx, chi_ecm_dx, ct);
            #elif ( MODEL == KSMD  )
            gpuStep_explicit(device1, device2, id, 0.0     , chi_ecm_dx, ct);
            #elif ( MODEL == KSC   || MODEL == TEST)
            gpuStep_explicit(device1, device2, id, chi_b_dx, 0.0       , ct);
            #elif ( MODEL == FK    )
            gpuStep_explicit(device1, device2, id, 0.0     , 0.0       , ct);
            #endif
			ExportDT(host,device1);

			step_dt += *host.new_dt;
                  
		}
		#endif
		step_dt = 0.0;

		//IMPLICIT STEP
        #if   ( MODEL == KSCMD ) // not ready yet
        gpuStep_implicit(device1, device2, cfl_a, cfl_b, cfl_c);
        #elif ( MODEL == KSMD  ) // not ready yet
        gpuStep_implicit(device1, device2, cfl_a, 0.0  , cfl_c);
        #elif ( MODEL == KSC  || MODEL == TEST )
        gpuStep_implicit(device1, device2, cfl_a, cfl_b, 0.0  );
		ExportIter(host, device1);
        //test_CD(device1, device2, id);
        #elif ( MODEL == FK    ) // not ready yet
        gpuStep_implicit(device1, device2, cfl_a, 0.0  , 0.0  );
        #endif

		caStep(cells, cellsTmp, casp, host, id);
		ct ++;

		//EXPLICIT STEP
        #if ( MODEL == FK || MODEL == KSCMD || MODEL == KSMD || MODEL == KSC || MODEL == TEST)
        while(step_dt - 0.5*id.dt_imp < -1.0e-12)
        {
            ct ++;
            #if   ( MODEL == KSCMD )
            gpuStep_explicit(device1, device2, id, chi_b_dx, chi_ecm_dx, ct);
            #elif ( MODEL == KSMD  )
            gpuStep_explicit(device1, device2, id, 0.0     , chi_ecm_dx, ct);
            #elif ( MODEL == KSC   || MODEL == TEST)
            gpuStep_explicit(device1, device2, id, chi_b_dx, 0.0       , ct);
            #elif ( MODEL == FK    )
            gpuStep_explicit(device1, device2, id, 0.0     , 0.0       , ct);
            #endif
            ExportDT(host,device1);
            step_dt += *host.new_dt;    
        }
        #endif
        step_dt = 0.0;
		

		if((i-exportInterval[c]) > - 1.0e-6)
		{
        	cudaEventRecord(stop, 0);		
            cudaEventSynchronize(stop);
            float elapsed;
            cudaEventElapsedTime(&elapsed, start, stop); //gives in milliseconds
            elapsed /= 1000.0f; //convert to seconds
            float mpoints = (float)(exportInterval[c]/id.dt_imp) * host.xmax * host.ymax * host.zmax / 1.0e6f;//reportInterval , exportInterval[c]                  
            printf("t_sim = %.2f: %f Mpoints/s \n", i*id.dt_imp, mpoints / elapsed);
            printf("iteration = %d\n",ct);
            run_time += elapsed;

         	ExportCheckpoint("a_double","b_double","ecm_double",host,device1,ceil(i*id.dt_imp),1);
            ExportCA(fca,id,cells,i);
            casp.adhes ++;
			if(casp.adhes >=5) {casp.adhes=5;}
           	casp.delaydiv --;
            casp.phenchange --;
			if(c >= 10){ casp.spdeath += 0.06*id.dt_imp; }
            c++;

            cudaEventRecord(start, 0);
		}
		ExportCellCount(fna,id,cells,i);
	}
	printf("elapsed time = %lf seconds\n",run_time);

   	//FINISH
	gpuClose(device1, device2);
	hostClose(host);
    cudaEventDestroy(start);
    cudaEventDestroy(stop);

	fclose(fna);
	fclose(fca);

}//end

