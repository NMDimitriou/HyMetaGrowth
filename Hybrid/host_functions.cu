#include "declarations.cuh"

//allocate memory on host
void hostInitialize(DataArray& host)
{
	int size      = total * sizeof(double);
	double dt_min = 1e-08;

#if ( MODEL == FK || MODEL == KSC || MODEL == KSMD || MODEL == KSCMD || MODEL == TEST)
	host.a	     = (double*) malloc(size);
	host.a_rows  = (double*) malloc(size);
	host.a_cols  = (double*) malloc(size);
	host.a_intm  = (double*) malloc(size);
	if (host.a           == NULL) Error("Host.a memory allocation failed\n");
    if (host.a_rows      == NULL) Error("Host.a_rows memory allocation failed\n");
    if (host.a_cols      == NULL) Error("Host.a_cols memory allocation failed\n");
	if (host.a_intm      == NULL) Error("Host.a_intm memory allocation failed\n");
#endif

#if ( MODEL == KSC || MODEL == KSCMD || MODEL == TEST)
	host.b       = (double*) malloc(size);
    host.b_rows  = (double*) malloc(size);
    host.b_cols  = (double*) malloc(size);
	host.del_b   = (double*) calloc(total,sizeof(double));
	host.b_intm  = (double*) malloc(size);
	if (host.b           == NULL) Error("Host.b memory allocation failed\n");
    if (host.b_rows      == NULL) Error("Host.b_rows memory allocation failed\n");
    if (host.b_cols      == NULL) Error("Host.b_cols memory allocation failed\n");
    if (host.del_b       == NULL) Error("Host.del_b memory allocation failed\n");
	if (host.b_intm      == NULL) Error("Host.b_intm memory allocation failed\n");
#endif

#if ( MODEL == KSMD || MODEL == KSCMD)
	host.ecm       = (double*) malloc(size);
    host.ecm_rows  = (double*) malloc(size);
    host.ecm_cols  = (double*) malloc(size);
    host.del_ecm   = (double*) calloc(total,sizeof(double));
	host.ecm_intm  = (double*) malloc(size);
	if (host.ecm           == NULL) Error("Host.ecm memory allocation failed\n");
	if (host.ecm_rows      == NULL) Error("Host.ecm_rows memory allocation failed\n");
    if (host.ecm_cols      == NULL) Error("Host.ecm_cols memory allocation failed\n");
    if (host.del_ecm       == NULL) Error("Host.del_ecm memory allocation failed\n");
	if (host.ecm_intm    == NULL) Error("Host.ecm_intm memory allocation failed\n");
#endif
	host.d_max_block = (double*) malloc(256*sizeof(double));
	if (host.d_max_block == NULL) Error("Host.d_max_block memory allocation failed\n");

	cudaMallocHost((void**)&host.dt    , sizeof(double));
    cudaMallocHost((void**)&host.new_dt, sizeof(double));

#if ( MODEL == FK || MODEL == KSC || MODEL == KSMD || MODEL == KSCMD)
	char filename[4096];
	sprintf(filename, "IC/corr_dens_%s_D0.bin", host.dat);
	
	FILE *f = fopen(filename,"rb");
	if (f!=NULL)
	{
	     	printf("Loading Initial Conditions from %s ...\n", filename);
	     	fread(host.a,sizeof(double),total,f);
	}else
	{
	     printf("Initial Conditions not found. \r\n");
	     exit(-1);
	}
	fclose(f);
	memcpy(host.a_rows,host.a,size);
    memcpy(host.a_cols,host.a,size);
    double* max_a = max_element(host.a, host.a+int(zsize*host.zmax));
#elif ( MODEL == TEST)
	for (int k=0; k<host.zmax; k++) {
    for (int j=0; j<host.ymax; j++) {
    for (int i=0; i<host.xmax; i++) {
	    int ind = i + j*host.xmax + k*host.xmax*host.ymax;
		if(i>2.0*host.xmax/5.0 && i<3.0*host.xmax/5.0 && j>2.0*host.ymax/5.0 && j<3.0*host.ymax/5.0 && k>3.0*host.zmax/5.0 && k<4.0*host.zmax/5.0)        
		{        
			host.a[ind]=0.1;
        }else{
         		host.a[ind]=0.0;
        }
	}
    }  
	}
	printf("This is a test.\n");
	memcpy(host.a_rows,host.a,size);
        memcpy(host.a_cols,host.a,size);
        double* max_a = max_element(host.a, host.a+int(zsize*host.zmax));
#endif

#if ( MODEL == KSMD || MODEL == KSCMD)
	char fnameECM[4096];
	sprintf(fnameECM, "../../IC/ecm_%s_D0.bin", host.dat);
    FILE *ff = fopen(fnameECM,"rb");
    if (ff!=NULL)
    {
    	printf("Loading Initial Conditions ... \r\n");
       	fread(host.ecm,sizeof(double),total,ff);
    }else{
    	printf("Initial Conditions for ECM not found. \r\n");
       	exit(-1);
	}
	fclose(ff);
    memcpy(host.ecm_rows,host.ecm,size);
	memcpy(host.ecm_cols,host.ecm,size);
#endif
#if ( MODEL == KSC || MODEL == KSCMD)
	double count  = 0.0;
    int c         = 0;
    for (int j=0; j<host.zmax;j++){
    for (int i=0; i<zsize; i++){
        	//host.b[c] = (1.0-0.005681818*count);
            	//host.b[c] = exp(-count/20.0)*host.a[c]/(*max_a);
            	//host.b[c] = host.a[c];

        if(host.a[c]>1.0e-12){
            		#if   ( IC == EXP )
                	host.b[c] = exp(-count/50.0);
                	#elif ( IC == LIN )
                	host.b[c] = (1.0-0.005681818*count);
                	#elif ( IC == ID  )
                	host.b[c] = host.a[c];
                	#endif
		}else{
            		host.b[c] = 0.0;
       	}
        c ++;	
	}
	count ++;    
	}
	memcpy(host.b_rows,host.b,size);
    memcpy(host.b_cols,host.b,size);
#elif ( MODEL == TEST )
	double count  = 0.0;
	int c         = 0;
	for (int j=0; j<host.zmax;j++)
	{
    	for (int i=0; i<zsize; i++)
        {
        	host.b[c] = (1.0-0.005681818*count);
            c ++;
    	}
        count ++;   
	}
	memcpy(host.b_rows,host.b,size);
    	memcpy(host.b_cols,host.b,size);
#endif
	for(int i=0; i<256; i++)
    	{
        	host.d_max_block[i]=100.0;
    	}

    	*host.dt = dt_min; //a "small" initial dt	
}



//free memory on host
void hostClose(DataArray& host)
{
#if   ( MODEL == FK || MODEL == KSC || MODEL == KSMD || MODEL == KSCMD || MODEL == TEST)
	free(host.a    );
	free(host.a_rows);
    free(host.a_cols);
	free(host.a_intm);
#endif
#if   ( MODEL == KSC  || MODEL == KSCMD || MODEL == TEST)
	free(host.b     );
	free(host.b_rows);
	free(host.b_cols);
	free(host.del_b );
	free(host.b_intm);
#endif
#if   ( MODEL == KSMD || MODEL == KSCMD )
	free(host.ecm     );
    free(host.ecm_rows);
    free(host.ecm_cols);
    free(host.del_ecm );
	free(host.ecm_intm);
#endif        
	cudaFreeHost(host.dt);    
	cudaFreeHost(host.new_dt);
	free(host.d_max_block);
}
