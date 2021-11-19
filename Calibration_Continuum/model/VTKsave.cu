#include "declarations.cuh"

#ifdef SAVE_DATA

template <typename T>
void SwapEnd(T& var)
{
  char* varArray = reinterpret_cast<char*>(&var);
  for(long i = 0; i < static_cast<long>(sizeof(var)/2); i++)
    std::swap(varArray[sizeof(var) - 1 - i],varArray[i]);
}


void VTKsave(char fname[4096], int xmax, int ymax, int zmax, double* a)
{
    float tmp;
    ofstream vtkstream;
    vtkstream.open(fname, std::ios::out | std::ios::app | std::ios::binary);
    if (vtkstream) {
        //FILE *fs=fopen("test.vtk","w");
        //if(fs==NULL){printf("file vtk, failed to open\n");}
        vtkstream << "# vtk DataFile Version 2.0"<<"\n";
        vtkstream << "Volume data"<<"\n";
        vtkstream << "BINARY"<<"\n";
        vtkstream << "DATASET STRUCTURED_POINTS"<<"\n";
        vtkstream << "DIMENSIONS " <<  xmax << " " << ymax << " " << zmax <<"\n";
        vtkstream << "ASPECT_RATIO " << 1 << " " << 1 << " " << 1 <<"\n";
        vtkstream << "ORIGIN " << 0 << " " << 0 << " " << 0 <<"\n";
        vtkstream << "POINT_DATA " << total <<"\n";
        vtkstream << "SCALARS Density float " << 1 <<"\n";
        vtkstream << "LOOKUP_TABLE default" <<"\n";
        for(unsigned int i=0; i<total; i++)
        {
            tmp=float(a[i]);
            SwapEnd(tmp);
            vtkstream.write((char*)&tmp, sizeof(float));
        }
    }
    vtkstream.close();

}
#endif
