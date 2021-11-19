#ifndef DECLARATIONS_H
#define DECLARATIONS_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <math.h>
#include <algorithm>
#include <time.h>
#include <random>
#include <iterator>
#include <vector>

/*
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>
*/

#define NONE            0

#define FK              1
#define KSC             2
#define KSMD            3
#define KSCMD           4
#define TEST            5

#define CALIBRATION     1
#define SA              2
#define HYBRID		3

#define EXP             1
#define LIN             2
#define ID              3


using namespace std;

//structure to hold size of simulation and arrays on the host
typedef struct DataArray
{
	int xmax;
	int ymax;
	int zmax;

#if (MODEL == FK || MODEL == KSC || MODEL == KSMD || MODEL == KSCMD || MODEL == TEST)   
 	double* a; //cell density
    	double* a_intm; //intermediate values of a
    	//Parameters for tridiagonal solver
    	double* a_rows;
    	double* a_cols;
    	double* a_ax;
    	double* a_ay;
    	double* a_az;
    	double* a_bx;
    	double* a_by;
    	double* a_bz;
    	double* a_cx;
    	double* a_cy;
    	double* a_cz;
#endif
#if (MODEL == KSC || MODEL == KSCMD || MODEL == TEST)
    	double* b; //signal density
    	double* b_intm;
    	//Parameters for tridiagonal solver
    	double* b_rows;
    	double* b_cols;
    	double* b_ax;
    	double* b_ay;
    	double* b_az;
    	double* b_bx;
    	double* b_by;
    	double* b_bz;
    	double* b_cx;
    	double* b_cy;
    	double* b_cz;
    	double* del_b;  //signal density gradient
#endif
#if ( MODEL == KSMD || MODEL == KSCMD )
	double* ecm; //ecm density
   	double* ecm_intm;
    	double* ecm_rows;
    	double* ecm_cols;
    	double* ecm_ax;
    	double* ecm_ay;
    	double* ecm_az;
    	double* ecm_bx;
    	double* ecm_by;
    	double* ecm_bz;
    	double* ecm_cx;
    	double* ecm_cy;
    	double* ecm_cz;
    	double* del_ecm;//ECM density gradient
#endif
	double* d_max_block;
    	double* dt; //adaptive time step for the explicit part
    	double* new_dt; //stores the adaptive time-step in the host
    	char* dat;	
} DataArray;

//Initial data for the simulation
typedef struct InitialData
{
#if ( MODEL == FK || MODEL == KSC || MODEL == KSMD || MODEL == KSCMD || MODEL == TEST)
    	double da; //Diffusion coefficient for species A
    	double s ; //Growth coefficient for species A
	double scale; //division scaling factor from exponential with base e to base 2
#endif
#if ( MODEL == KSC || MODEL == KSCMD  || MODEL == TEST)
    	double db; //Diffusion of signal factors
    	double chi; //Chemotactic constant
    	double r ; //Increase of signal factors
#endif
#if ( MODEL == KSMD || MODEL == KSCMD )
    	double dc; //Degradation of ECM as diffusion process
    	double chi_ecm; //Chemotactic constant that drives cells towards max ecm density
#endif
//  double ka; //Cell death term
//  double dup; //uptake term 
    	double dt_imp; //Timestep for the implicit part
   	double cfl;
    	double dx; //Grid spacing
    	int itmax; //Max number of iteration steps
#if ( STUDY == SA )     
    	double tot_num;
    	double bot_num;
#endif
	
} InitialData;


//Structure to hold data for space of CA
typedef struct CAArray {
   	bool* lattice;
	bool* phen   ;
        vector<int> indcNeigh;
        int adhes     ;
        int agediv    ;
        int delaydiv  ;
	int divpot    ;
	int phenchange;
        double spdeath;
	int count_cell;
} CAArray;

//Structure to hold data of cells of CA
typedef struct Cell {
	int  place   ;
   	int  age     ;
   	int  ndiv    ;
   	bool is_alive;
} Cell;

//Kernel configuration parameter for boundary conditions. All dimensions must be multiple of this.
#define BLOCKSIZE 16
#define BDIMX 16
#define BDIMY 16
#define BDIMZ 16
#define radius 2
#define BXMAX 480
#define BYMAX 480
#define BZMAX 176

// (19-point central differences) 3D Laplace-computing kernel configuration parameters. Dimensions of a kernel is how many threads it contains
#define CELLW	8	//width of thread block
#define CELLH	8	//height of thred block
#define CELLD	8	//depth of thread block. At least 3, but value is arbitrary since CELLW*CELLH is already multiple of the warp size (32)
//NOTE: On devices with compute capability 1.0 and 1.1 (1st generation devices) CELLW should be 16, CELLH and CELLD: 5.
//      this causes reading more data in every step, but the read will be coalesced!
//      On 2nd generation devices all global memory operations are coalesced, block dimensions are chosen to
//      minimize read redundancy. 

//fundamental variables for addressing 3D coordinates in 1D arrays
extern int zsize, ysize, total;

//functions in main.cu
void CudaCheck();
void Error(const char* text);

//functions in host_functions.cu
void hostInitialize(DataArray& host);
void hostClose(DataArray& host);

//functions in ca_functions.cu
void initialize_ca(DataArray& host, CAArray& casp, vector<Cell>& cells);
int returnEmptyPlace(CAArray& casp, int indx);
int returnEmptyPlace2D(CAArray& casp, int indx);
int returnNeigh(CAArray& casp, int indx);
int returnNeighPhen(CAArray& casp, int indx);
int returnPosChance(CAArray& casp, DataArray& host, InitialData& id, int indx);
void caStep(vector<Cell>& cells, vector<Cell>& cellsTmp, CAArray& casp, DataArray& host, InitialData& id);
void ExportCellCount(FILE* fna, InitialData& id, vector<Cell>& cells, int i);
void ExportCA (FILE* fca, InitialData& id, vector<Cell>& cells, int i);

//functions in device_functions.cu
void gpuClose(DataArray& device1, DataArray& device2);
void gpuStep_explicit(DataArray& device1, DataArray& device2, InitialData& id, double chi_dx, double chi_ecm_dx, int iter);
void gpuStep_implicit(DataArray& device1, DataArray& device2, double cfl_a, double cfl_b, double cfl_c);
void gpuInitialize(InitialData& id, DataArray& host, DataArray& device1, DataArray& device2);
void ExportDT( DataArray& host, DataArray& device1 );
void ExportIter(DataArray& host, DataArray& device1);
void ExportCheckpoint(const char* name_a, const char* name_b, const char* name_c, DataArray& host, DataArray& device1, int l, int printout);

//functions in readData.cu
void readData(InitialData& id, char* argv);

//functions in scaleGrowth.cu
//int expb2_f (const gsl_vector * x, void *data, gsl_vector * f);
//int expb2_df (const gsl_vector * x, void *data, gsl_matrix * J);
//void callback(const size_t iter, void *params, const gsl_multifit_nlinear_workspace *w);
//void scaleGrowth(DataArray& host, InitialData& id);

//functions in simulate.cu
void simulate(DataArray& host, DataArray& device1, DataArray& device2, CAArray& casp, vector<Cell>& cells, vector<Cell>& cellsTmp, InitialData& id);

//function to save the data
template <typename T>
void SwapEnd(T& var);
void VTKsave(char fname[4096], int xmax, int ymax, int zmax, double* a);


#endif //declarations_h
