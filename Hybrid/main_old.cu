#include "declarations.cuh"


//data arrays and simulation parameters
DataArray host, device1, device2;
InitialData id;
CAArray casp;
//Cell  currCell, newCell;
vector<Cell> cells;
vector<Cell> cellsTmp;

//variables for memory allocation and addressing
int zsize, ysize, total;

//mixmax_engine mxmx;
//uniform_real_distribution<float> drandom (0.0f,1.0f) ;

//Function to check for error after a CUDA call.
void CudaCheck()
{
        cudaError_t er = cudaGetLastError();
        if (er!=0)
        {       
                printf("CUDA Error: %s\n", cudaGetErrorString(er));
#ifdef _WIN32
                system("PAUSE");
#endif          
                exit(-1);
        };
};

//function to halt program with an error message
void Error(const char* text)
{
        printf("%s\n", text);
#ifdef _WIN32
                system("PAUSE");
#endif
        exit(-1);
}

//main program
int main(int argc, char *argv[])
{

        host.xmax=480;
        host.ymax=480;
        host.zmax=176;

        printf("Argument supplied is: %s\n", argv[1]);

	//compute 3D addressing variables
        total = host.xmax * host.ymax * host.zmax;
        zsize = host.xmax * host.ymax;
        ysize = host.xmax;
	//mxmx.seed(seedby_urandom());

        //define parameters
	double scale = 1.44245;
        id.da        = 0.00417;
        id.db        = 0.00456;
        id.chi       = 0.0166;
	id.s         = 0.201;
        id.r         = 0.154;
        id.dx        = 2.5/(ysize -1);
	id.cfl       = 0.15;
	double dt_a  = id.cfl*id.dx*id.dx/id.da;
	double dt_b  = id.cfl*id.dx*id.dx/id.db;
	id.dt_imp    = min(dt_a,dt_b);
	id.dt_imp    = min(id.dt_imp,0.1);
	//id.dt_imp    = max(id.dt_imp,0.001);
	id.itmax     = 14.0/id.dt_imp;
	double cfl_a = id.da*id.dt_imp/(id.dx*id.dx);
	double cfl_b = id.db*id.dt_imp/(id.dx*id.dx);
//      id.ka      = 0.0;
	casp.adhes      = -6;
	casp.delaydiv   = 2;
	casp.phenchange = 5;
	casp.divpot   = 10000*24;

	printf("\n########################################## \n \n");
        printf("Initializing 3D-KS model \n(0,2.5)x(0,2.5)x(0,0.9) mm \n");
        printf("Total time: 14 days with %d implicit time steps\n", id.itmax);
        printf("dx,y,z = %f\n", id.dx);
        printf("Da = %lf, s = %lf, Chi = %lf\n", id.da, id.s, id.chi);
	printf("Db = %lf, r = %lf\n", id.db, id.r);
        printf("\n########################################## \n");

        //INITIALIZE
        hostInitialize(host);
	if(*host.dt >= id.dt_imp/4.0) *host.dt = id.dt_imp/4.0;

	casp.agediv   = 1.0/(scale*id.s*id.dt_imp);
        casp.spdeath  =-1.0*(*host.dt);
        printf("Age div : %d\n", casp.agediv);

	initialize_ca(host, casp, cells);
	FILE *fca;
        fca = fopen("ca_res_AE.csv", "wa");
        if (fca==0) { printf("  error creating saving file\n"); return -1; }

	ExportCA(fca,id,cells,0);

        gpuInitialize(id, host, device1, device2);

        //Export time points for movie with 24 frames per day
//	int reportInterval = 1000 ;
/*	int   total_tp = 336;
        float exprt_dt = 1.0f/24.0f;
        int exportInterval[total_tp];
        for (int j=1; j<=total_tp; j++)
        {
                exportInterval[j] = int(j*exprt_dt/id.dt_imp);
                if (j<6) printf("exportInterval[%d] = %d \n", j, exportInterval[j]);
        }
*/
//	int exportInterval[6] = {int(2.0f/id.dt_imp),int(5.0f/id.dt_imp),int(7.0f/id.dt_imp),
//				 int(9.0f/id.dt_imp),int(12.0f/id.dt_imp),int(14.0f/id.dt_imp)} ; 

	int exportInterval[14] = {int(1.0f/id.dt_imp),int(2.0f/id.dt_imp),int(3.0f/id.dt_imp),int(4.0f/id.dt_imp),
				  int(5.0f/id.dt_imp),int(6.0f/id.dt_imp),int(7.0f/id.dt_imp),int(8.0f/id.dt_imp),
				  int(9.0f/id.dt_imp),int(10.0f/id.dt_imp),int(11.0f/id.dt_imp),int(12.0f/id.dt_imp),
				  int(13.0f/id.dt_imp),int(14.0f/id.dt_imp)};

//	int exportInterval[6] = {2,5,7,9,12,14} ;

        int    c = 0  ;
        double step_dt= *host.dt;
	int ct = 0;
        // Export Day 0.
        ExportCheckpoint("a_double","b_double",host,device1,0,1);
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

/*		while(step_dt <= id.dt_imp)
        	{
                        gpuStep_explicit(device1, device2);
                        ExportDT(host,device1);
                        step_dt += *host.new_dt;
			//printf("min dt = %1.2e \n", *host.new_dt);
                }
		step_dt = *host.dt;
*/
		
		//Strang splitting
		while(step_dt <= id.dt_imp/2.0)
		{
			gpuStep_explicit(device1, device2);
                	ExportDT(host,device1);
			step_dt += *host.new_dt;
			//ExportIter(host, device1);
                	//caStep(cells, cellsTmp, casp, host, id);
			ct ++;
		}
	
		gpuStep_implicit(device1, device2, cfl_a, cfl_b);

		ExportIter(host, device1);
                caStep(cells, cellsTmp, casp, host, id);

		ct ++;
	
		step_dt = *host.dt;
		while(step_dt <= id.dt_imp/2.0)
                {
                        gpuStep_explicit(device1, device2);
                        ExportDT(host,device1);
                        step_dt += *host.new_dt;
			//ExportIter(host, device1);
                	//caStep(cells, cellsTmp, casp, host, id);
			ct ++;
                }
		step_dt = *host.dt;
	
	//	printf("min dt = %1.2e \n", *host.new_dt);

                if((i-exportInterval[c]) > - 1.0e-6)
                {
                        cudaEventRecord(stop, 0);
                        cudaEventSynchronize(stop);
                        float elapsed;
                        cudaEventElapsedTime(&elapsed, start, stop); //gives in milliseconds
                        elapsed /= 1000.0f; //convert to seconds
                        float mpoints = (float)(exportInterval[c]/step_dt) * host.xmax * host.ymax * host.zmax / 1.0e6f;//reportInterval , exportInterval[c]                  
                        printf("t_sim = %.2f: %f Mpoints/s \n", i*id.dt_imp, mpoints / elapsed);
			printf("iteration = %d\n",ct);
			run_time += elapsed;

			ExportCheckpoint("a_double","b_double",host,device1,ceil(i*id.dt_imp),1);
			ExportCA(fca,id,cells,i);
			casp.adhes ++;
			casp.delaydiv --;
			casp.phenchange --;
                        c++;

			cudaEventRecord(start, 0);
                }
        }
	printf("elapsed time = %lf seconds\n",run_time);
        //FINISH
        gpuClose(device1, device2);
        hostClose(host);
        cudaEventDestroy(start);
        cudaEventDestroy(stop);
	fclose(fca);

        return 0;
}

