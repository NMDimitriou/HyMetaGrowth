#include "declarations.cuh"

//data arrays and simulation parameters
DataArray host, device1, device2;
InitialData id, curr;
CAArray casp;
vector<Cell> cells;
vector<Cell> cellsTmp;

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
	host.xmax=480;
  	host.ymax=480;
  	host.zmax=176;
	//compute 3D addressing variables
  	total = host.xmax * host.ymax * host.zmax;
  	zsize = host.xmax * host.ymax;
  	ysize = host.xmax;
	host.dat = argv[2];
	
	//Read the data
	readData(id,argv[1]);

	//Simulate
  	simulate(host,device1,device2,casp,cells,cellsTmp,id);

  	return 0;
}



