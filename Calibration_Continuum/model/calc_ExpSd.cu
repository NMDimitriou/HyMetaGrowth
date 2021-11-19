#include "declarations.cuh"

/* Calculate max-min or SD of experimental data
 * Nikolaos Dimitriou, McGill, 2021
*/

double calc_ExpSd(double* expdata)
{
	double sum=0.0;
	double mean=0.0;
	double sd=0.0;
	double max=0.0,min=100.0;

	for(int i=0; i<total; i++)
	{
		sum += expdata[i];
		if(max < expdata[i]) { max = expdata[i]; }
		if(min > expdata[i]) { min = expdata[i]; }
	}

	mean = sum/total;

	for(int i=0; i<total; i++)
        {
                sd += (expdata[i]-mean)*(expdata[i]-mean);
        }

	//return sqrt(sd/(total-1));
	return (max-min) ;

}
