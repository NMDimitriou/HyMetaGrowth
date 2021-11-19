#include "declarations.cuh"

/* Sets the GPU to be used for the simulation
 * Nikolaos Dimitriou, McGill, 2021
*/

void setDevice(InitialData& id, int tot_gpus)
{
	int devID, nDevices;
    cudaGetDeviceCount(&nDevices);

#if ( STUDY == SA )
	printf(" \n");
    printf("Number of GPUs in this node = %d\nrun = %d\n",nDevices,id.dev);
    printf(" \n");

	devID  = (id.dev) % nDevices;
	cudaSetDevice(devID);
    int myDevId = -1;
    cudaGetDevice( &myDevId ) ;
    printf(" \n");
    printf("Set Device = %d\n", myDevId);
    printf(" \n");
#else

	printf(" \n");
    printf("Number of GPUs in this node = %d\nrank = %d\n",nDevices,id.dev);
	printf(" \n");

	char hostname[HOST_NAME_MAX + 1];
    gethostname(hostname, HOST_NAME_MAX + 1);

	char filename[4096];
   	const char* name = "../devID";
    sprintf(filename, "%s_%d_%s.txt", name, id.pid, hostname);

    for (int i = 0; i < nDevices; i++) {
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, i);
        printf("Device Number: %d\n", i);
        printf("  Device name: %s\n", prop.name);
        printf("  Memory Clock Rate (KHz): %d\n",
                prop.memoryClockRate);
        printf("  Memory Bus Width (bits): %d\n",
                prop.memoryBusWidth);
        printf("  Peak Memory Bandwidth (GB/s): %f\n\n",
                2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6);
    }

	if(id.dev < tot_gpus)
	{
		devID = (id.dev) % nDevices;

		cudaSetDevice(devID);
        int myDevId = -1;
        cudaGetDevice( &myDevId ) ;
        printf(" \n");
        printf("Set Device = %d\n", myDevId);
        printf(" \n");

		FILE* fdev;
        fdev = fopen(filename, "w");
        fprintf(fdev, "%d", devID);
		fclose(fdev);


	}
	else
	{
		FILE* fdev;
        fdev = fopen(filename, "r");
        fscanf(fdev, "%d", &devID);
        fclose(fdev);

		cudaSetDevice(devID);
        int myDevId = -1;
        cudaGetDevice( &myDevId ) ;
        printf(" \n");
        printf("Set Device = %d\n", myDevId);
        printf(" \n");
	}

/*
	int fingpu[4];
        int dev_free = -1;
        //Read used memory of gpus
        FILE* fgpu;
        fgpu = fopen("gpu_id.txt","r");
        for(int i=0;i<nDevices;i++)
        {
        	fscanf(fgpu, "%d",&fingpu[i]);
                printf("GPU %d has used memory %d MiB\n",i, fingpu[i]);

                if(fingpu[i] < 500)
                {
                	dev_free = i ; //devices with free memory
                }

	}
        fclose(fgpu);
        printf(" \n");

	if(dev_free >=0) 
	{
		devID = dev_free;
	}

	cudaSetDevice(devID);
        int myDevId = -1;
        cudaGetDevice( &myDevId ) ;
        printf(" \n");
        printf("Set Device = %d\n", myDevId);
        printf(" \n");
*/
                // Check if selected device has free memory
/*                if(fingpu[devID]>500 && dev_free >= 0)
                {
                        printf("Potential conflict with GPUs, selecting an empty device. \n");
                        devID = dev_free;
                }

                if(fingpu[devID]>500 && dev_free < 0)
                {
                        printf("There is no available device. Exiting.");
                        exit(-1);
                }
*/
/*
		cudaSetDevice(devID);
                int myDevId = -1;
                cudaGetDevice( &myDevId ) ;
                printf(" \n");
                printf("Set Device = %d\n", myDevId);
                printf(" \n");

	}
*/

#endif

}
