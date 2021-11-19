#include "declarations.cuh"

//data arrays and simulation parameters
DataArray host, device1, device2;

ExpData expdat;

InitialData id, curr;

//variables for memory allocation and addressing
int zsize, ysize, total;

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
int main(int argc, char* argv[])
{

	if (argc < 3) { // We expect 3 arguments: the parameter file, the name of the dataset
			// and the total number of gpus that will be used for parameter estimation
			// or sensitivity analysis
		std::cerr << "Usage: ./main [arg1] [arg2] [arg3]" << std::endl;
        std::cerr << "	     " << "[arg1] parameter file (params.txt)" << std::endl;
		std::cerr << "	     " << "[arg2] name of dataset (e.g. AE)" << std::endl;
		std::cerr << "	     " << "[arg3] number of GPUs (>1 works only for calibration or sensitivity test)" << std::endl;
		std::cerr << "Before running the simulation make sure that you have loaded the ICs or experimental data." << std::endl;
		std::cerr << "To do this, go to parent directory (cd parent) and run" << std::endl;
		std::cerr << "	      ./main_parent start [name of dataset, e.g. AE]" << std::endl;
        return 1;
    }

    host.xmax=480;
    host.ymax=480;
    host.zmax=176;
	host.dat =argv[2];
	//compute 3D addressing variables
    total = host.xmax * host.ymax * host.zmax;
    zsize = host.xmax * host.ymax;
    ysize = host.xmax;

	//Read the data
	printf("Reading model parameters\n");
	readData(id,argv[1]);

	//Set device
	printf("Setting GPU\n");
	int num_gpus = atoi(argv[3]);
	setDevice(id,num_gpus);

    //Initialize
	#if (STUDY == CALIBRATION)
	printf("Loading experimental data\n");
    hostExpInit(expdat);
	#endif

	//Simulate
	printf("Simulation starts...\n");
    curr.SSE= simulate(host,device1,device2,expdat,id);

	//Save
	#if (STUDY == CALIBRATION)
    FILE* fin;
    fin = fopen("loglike.txt", "w");
    fprintf(fin, "%.15f\n", curr.SSE);
    fflush( fin);
    fclose(fin);
	#elif (STUDY == SA)
    FILE* fin;
    fin = fopen("response.txt", "w");
    fprintf(fin, "%.15f\n", id.tot_num);
	fprintf(fin, "%.15f"  , id.bot_num);
    fflush( fin);
    fclose(fin);
	#endif


    return 0;
}

