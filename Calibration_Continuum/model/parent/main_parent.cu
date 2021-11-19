#include "declarations_parent.cuh"

//data arrays and simulation parameters
DataArray host;//, device1, device2;
#if ( STUDY == CALIBRATION )
ExpData expdat;
#endif
//InitialData id, curr;

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

	if (argc < 3) { // We expect 3 arguments: start or stop argument (loads or unloads the data)
                        // the path of the data
                        // the name of the dataset
                std::cerr << "Usage: ./main_parent [arg1] [arg2] [arg3]" << std::endl;
                std::cerr << "       " << "[arg1] start or stop (load or unloads the data)" << std::endl;
                std::cerr << "       " << "[arg2] path of the data" << std::endl;
                std::cerr << "       " << "[arg3] name of dataset" << std::endl;
		std::cerr << "Reminder: the file names of the data should follow the pattern:" << std::endl; 
		std::cerr << " 		corr_dens_[name of dataset]_[time-point as D0,D2 etc.].bin" << std::endl;
        	return 1;
        }

        host.xmax=480;
        host.ymax=480;
        host.zmax=176;
	host.path= argv[2];
	host.dat = argv[3];
	//compute 3D addressing variables
        total = host.xmax * host.ymax * host.zmax;
        zsize = host.xmax * host.ymax;
        ysize = host.xmax;

	char hostname[HOST_NAME_MAX + 1];
  	gethostname(hostname, HOST_NAME_MAX + 1);
  	printf("hostname: %s\n", hostname);
	printf("Dataset: %s\n", host.dat);


	if(!strcmp(argv[1],"start"))
	{
		hostInitializeParent(host,hostname);
		#if ( STUDY == CALIBRATION )
		hostExpInitParent(host,expdat,hostname);
		#endif
	}
	else if(!strcmp(argv[1],"stop"))//(ff==1)
	{
		printf("Erasing memory. \n");
		hostClose(host);
		#if ( STUDY == CALIBRATION )
		hostCloseExp(expdat);
		#endif
	}
	else 
	{
		printf("First argument must be start or stop. Exiting\n"); exit(-1);
	}

		

        return 0;
}

