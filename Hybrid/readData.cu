#include "declarations.cuh"

/* Reads the model parameters from params.txt file
 * Nikolaos Dimitriou, McGill 2021
*/ 

void readData(InitialData& id, char* argv)
{
#if ( MODEL == FK )
        float ff[2];
        FILE *f;
        f = fopen(argv, "r");
        for(int i=0;i<2;i++)
        {
                fscanf(f, "%f",&ff[i]);
        }
        fclose(f);
   
        id.da=ff[0]; id.s =ff[1]; 
#elif ( MODEL == KSC || MODEL == TEST)
	float ff[5];
        FILE *f;
        f = fopen(argv, "r");
        for(int i=0;i<5;i++)
        {
                fscanf(f, "%f",&ff[i]);
        }
        fclose(f);

        id.da=ff[0]; id.s =ff[1]; id.chi=ff[2]; id.db=ff[3]; id.r=ff[4]; 
#elif ( MODEL == KSMD )
	float ff[4];
        FILE *f;
        f = fopen(argv, "r");
        for(int i=0;i<5;i++)
        {
                fscanf(f, "%f",&ff[i]);
        }
        fclose(f);

        id.da=ff[0]; id.s =ff[1]; id.chi_ecm=ff[2];  id.dc=ff[3];  
#elif ( MODEL == KSCMD )
	float ff[7];
        FILE *f;
        f = fopen(argv, "r");
        for(int i=0;i<7;i++)
        {
                fscanf(f, "%f",&ff[i]);
        }
        fclose(f);

        id.da=ff[0]; id.s =ff[1]; id.chi=ff[2]; id.chi_ecm = ff[3]; id.db=ff[4]; id.r=ff[5]; id.dc=ff[6]; 
#endif
}
