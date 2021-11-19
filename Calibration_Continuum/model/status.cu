#include "declarations.cuh"

// Print status of simulation to be printed from Pi4U TMCMC
// -1 did not converge
// 0 Nan of Inf for SSE
// 1 Success
// Nikolaos Dimitriou, McGill, 2021

void status(int st)
{

	char filename[4096];
    const char* name = "status";
    sprintf(filename, "%s.txt",name);
	const char* message;

	switch(st)
  	{
        	case -1 :
		message="Did not converge.";
       		break;
                case 0 :
		message="SSE was Nan/Inf.";
		break;
		case 1 :
		message="Success.";
		break;
	}

        FILE* fst;
        fst = fopen(filename, "w");
        fprintf(fst, "%s", message);
        fclose(fst);
}

