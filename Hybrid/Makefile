# Make file for the hybrid models

IDIR=include

CC=g++
NVCC=nvcc
RM=rm -f

LDIR=$(CUDA_LIB_PATH)
LIBS=-lcuda -lcudart -lm -lpthread

ODIR=obj

_DEPS=declarations.cuh
DEPS=$(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ=main.o readData.o host_functions.o device_functions.o ca_functions.o simulate.o VTKsave.o #scaleGrowth.o
OBJ=$(patsubst %,$(ODIR)/%,$(_OBJ))

# model: FK/KSC/KSMD/KSCMD/TEST, IC: EXP/LIN/ID, study: CALIBRATION/SA/NONE
SIMU_OPTION += -DSTUDY=NONE
SIMU_OPTION += -DMODEL=KSC
#SIMU_OPTION += -DSAVE_DATA
#SIMU_OPTION += -DNRMSE
SIMU_OPTION += -DIC=EXP

CFLAGS=-I $(IDIR)
CXXFLAGS=-std=c++11 $(SIMU_OPTION) #`gsl-config --cflags --libs` 
CUFLAGS=--ptxas-options=-v -arch=sm_70

MAIN=main

$(ODIR)/%.o: %.cu $(DEPS)
	$(NVCC) -c -o $@ $< $(CFLAGS) $(CXXFLAGS) $(CUFLAGS)

$(MAIN):$(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(CXXFLAGS) $(LIBS)


clean:
	$(RM) $(ODIR)/*.o *~ core* $(INCDIR)/*~

clear:
	$(RM) *.raw *.csv *.out


