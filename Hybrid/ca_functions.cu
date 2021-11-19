#include "declarations.cuh"

//initialize CA
void initialize_ca(DataArray& host, CAArray& casp, vector<Cell>& cells) {

	int Nz = host.zmax;
	int size = total*sizeof(bool);
	casp.lattice = (bool*) malloc(size);
	casp.phen    = (bool*) malloc(size);

	srand(time(NULL)); // initialize random number generator

	if (casp.lattice == NULL) Error("Host memory allocation failed in CA\n");
	if (casp.phen    == NULL) Error("Host memory allocation failed in CA\n");

	casp.lattice[total] = {false};
	casp.phen[total]    = {false};

    for (int i=0; i<2*zsize; i++) {casp.lattice[i]=true; casp.phen[i]=true;};//filling base
    for (int i=zsize*(Nz-2); i<total; i++) {casp.lattice[i]=true;};//filling top
    for (int i=0; i<2*ysize; i++)
    {
        for (int j=i; j<=zsize*(Nz-2)+i; j += zsize) {casp.lattice[j]=true;}; //filling left boundary perpendicular to y-axis
    }
    for (int i=0; i<=zsize-2*ysize; i+=ysize)
    {
        for (int j=i; j<=zsize*(Nz-2)+i; j += zsize) {casp.lattice[j]=true;}; //filling back boundary perpendicular to x-axis
    }
    for (int i=ysize-2; i< ysize   ; i+=ysize)
    {
        for (int j=i; j<=zsize*(Nz-2)+i; j += zsize) {casp.lattice[j]=true;}; //filling front boundary perpendicular to x-axis
    }
    for (int i=zsize-2*ysize; i<zsize; i++)
    {
        for(int j=i; j< total ; j += zsize) {casp.lattice[j]=true;}; //filling right boundary perpendicular to y-axis
    }


	casp.indcNeigh    = {-2,2,-2*ysize,2*ysize,-2*zsize,2*zsize}; //neighborhood

    random_device rand_dev;
    mt19937 generator(rand_dev());

	char filename[4096];
    sprintf(filename, "IC/coord_idx_%s_D0.txt", host.dat);
    int n; //cell spatial index
    ifstream read(filename);
    while(read>>n)
    {
        int   rage = casp.agediv * (double)rand()/(double)RAND_MAX;//drandom(mxmx);//(double)rand()/(double)RAND_MAX;
        casp.lattice[n] = true; //initial cell distribution from indices file
        Cell initialCell = {n, rage, 0, true}; //cell state
        cells.push_back(initialCell);
		casp.count_cell ++;
    }
}//end of initialize_ca





int returnEmptyPlace(CAArray& casp, int indx) 
{
    int neigh[6], nF = 0;
    for(int j=0;j<6;j++) {//searching through neighborhood
        if (!casp.lattice[indx+casp.indcNeigh[j]]) {//if free spot
            neigh[nF] = indx+casp.indcNeigh[j]; //save the index
            nF++; //increase the number of found free spots
        }
    }
    if(nF) //selecting free spot at random
        return neigh[rand() % nF];
    else //no free spot
        return 0;
    
}//end of returnEmptyPlace




int returnEmptyPlace2D(CAArray& casp, int indx) {
        
    int neigh[4], nF = 0;
    for(int j=0;j<4;j++) {//searching through neighborhood
        if (!casp.lattice[indx+casp.indcNeigh[j]]) {//if free spot
            neigh[nF] = indx+casp.indcNeigh[j]; //save the index
            nF++; //increase the number of found free spots
        }
    }
    if(nF) //selecting free spot at random
        return neigh[rand() % nF];
    else //no free spot
        return 0;

}//end of returnEmptyPlace





int returnNeigh(CAArray& casp, int indx) {

    int nF = 0;
    for(int j=0;j<6;j++) {//searching through neighborhood
        if (casp.lattice[indx+casp.indcNeigh[j]]) //if neighbour
            nF++; //increase the number of found neighbours        
    }
    return nF;
}//end of returnNeigh




int returnNeighPhen(CAArray& casp, int indx) {
    
    int nF = 0;
    for(int j=0;j<6;j++) {//searching through neighborhood
    
        if (casp.lattice[indx+casp.indcNeigh[j]] && casp.phen[indx+casp.indcNeigh[j]]) //if neighbour
            nF ++ ; //increase the number of found neighbours
 
    }
    
    return nF;
}//end of returnNeighPhen




int returnPosChance(CAArray& casp, DataArray& host,InitialData& id, int indx) {

    
    int neigh = 0;
	int scl   = 2; //scale parameter accounts for the size of the cell
	double dt = id.dt_imp;//*host.dt;//id.dt_imp;

	// Generate random number
    random_device rand_dev;
    mt19937 generator(rand_dev());
    uniform_real_distribution<float> drandom(0.0f,1.0f);

//	double DaDt_dxdx   = id.da*id.dt_imp     /(id.dx*id.dx);
//	double Da6Dt_dxdx  = id.da*id.dt_imp*6.0f/(id.dx*id.dx);
//	double ChiDt_dxdx  = id.chi*id.dt_imp    /(id.dx*id.dx);
//	double ChiDt_4dxdx = id.chi*id.dt_imp    /(4.0f*id.dx*id.dx);

	double DaDt_dxdx   = id.da*dt     /(id.dx*id.dx);
//  double Da6Dt_dxdx  = id.da*dt*6.0 /(id.dx*id.dx);
//  double ChiDt_dxdx  = id.chi*dt    /(id.dx*id.dx);
    double ChiDt_4dxdx = id.chi*dt    /(4.0*id.dx*id.dx);

//	double F  = (host.b[indx+1]+host.b[indx-1]+host.b[indx+ysize]+
//  host.b[indx-ysize]+host.b[indx+zsize]+host.b[indx-zsize]-6.0f*host.b[indx]);

	float F1 = (host.b_intm[indx+1*scl    ]-host.b_intm[indx-1*scl    ]);
	float F2 = (host.b_intm[indx+ysize*scl]-host.b_intm[indx-ysize*scl]);
	float F3 = (host.b_intm[indx+zsize*scl]-host.b_intm[indx-zsize*scl]);

/*	double F1f = (host.b[indx+2      ]-host.b[indx        ]);
	double F1b = (host.b[indx        ]-host.b[indx-2      ]);
    double F2f = (host.b[indx+2*ysize]-host.b[indx        ]);
	double F2b = (host.b[indx        ]-host.b[indx-2*ysize]);
    double F3f = (host.b[indx+2*zsize]-host.b[indx        ]);
	double F3b = (host.b[indx        ]-host.b[indx-2*zsize]);
*/
/*
	double Fxi   = (host.b_intm[indx        ]-host.b_intm[indx-1      ]);
	double Fxip1 = (host.b_intm[indx+1      ]-host.b_intm[indx        ]); 
	double Fxim1 = (host.b_intm[indx-1      ]-host.b_intm[indx-2      ]);
	double Fyi   = (host.b_intm[indx        ]-host.b_intm[indx-ysize  ]);
    double Fyip1 = (host.b_intm[indx+ysize  ]-host.b_intm[indx        ]);
    double Fyim1 = (host.b_intm[indx-ysize  ]-host.b_intm[indx-2*ysize]);
	double Fzi   = (host.b_intm[indx        ]-host.b_intm[indx-zsize  ]);
    double Fzip1 = (host.b_intm[indx+zsize  ]-host.b_intm[indx        ]);
    double Fzim1 = (host.b_intm[indx-zsize  ]-host.b_intm[indx-2*zsize]);
*/	

//	printf("F=%.15f, F1=%.15f, F2=%.15f, F3=%.15f \n", F,F1,F2,F3);

//	float P0 =  1.0+dt*id.s-Da6Dt_dxdx - ChiDt_dxdx*(max(0.0,-Fxi)+max(0.0,-Fyi)+max(0.0,-Fzi)+
//  max(0.0, Fxi)+max(0.0, Fyi)+max(0.0, Fzi));  //(1.0 - Da6Dt_dxdx - ChiDt_dxdx*F);
    
    float P1 = max(0.0f, DaDt_dxdx - ChiDt_4dxdx*F1);//Back
    float P2 = max(0.0f, DaDt_dxdx + ChiDt_4dxdx*F1);//Front
    float P3 = max(0.0f, DaDt_dxdx - ChiDt_4dxdx*F2);//Left 
    float P4 = max(0.0f, DaDt_dxdx + ChiDt_4dxdx*F2);//Right
    float P5 = max(0.0f, DaDt_dxdx - ChiDt_4dxdx*F3);//Down 
    float P6 = max(0.0f, DaDt_dxdx + ChiDt_4dxdx*F3);//Up
	float P0 = 1.0f-(P1+P2+P3+P4+P5);


/*
	float P0 =  1.0+dt*id.s-Da6Dt_dxdx - ChiDt_dxdx*(max(0.0,-Fxi)+max(0.0,-Fyi)+max(0.0,-Fzi)+
	 max(0.0, Fxi)+max(0.0, Fyi)+max(0.0, Fzi));  //(1.0 - Da6Dt_dxdx - ChiDt_dxdx*F);
	float P1 = (DaDt_dxdx + ChiDt_dxdx*max(0.0,-Fxip1)); 
	float P2 = (DaDt_dxdx + ChiDt_dxdx*max(0.0, Fxim1)); 
	float P3 = (DaDt_dxdx + ChiDt_dxdx*max(0.0,-Fyip1)); 
	float P4 = (DaDt_dxdx + ChiDt_dxdx*max(0.0, Fyim1)); 
	float P5 = (DaDt_dxdx + ChiDt_dxdx*max(0.0,-Fzip1)); 
	float P6 = (DaDt_dxdx + ChiDt_dxdx*max(0.0, Fzim1)); 
*/
	if (P0 < 0.0f) Error("P0 is negative!\n");

//	printf("P0=%.5f, P1=%.5f, P2=%.5f, P3=%.5f, P4=%.5f, P5=%.5f, P6=%.5f\n",P0,P1,P2,P3,P4,P5,P6);

	float normR3D = P0+P1+P2+P3+P4+P5+P6;
	float normR2D = P0+P1+P2+P3+P4;
//	vector<float> Pi     = {P1,P2,P3,P4,P5,P6};
	vector<float> normRi3D = {P0/normR3D, 
                             (P0+P1)/normR3D, 
                             (P0+P1+P2)/normR3D, 
                             (P0+P1+P2+P3)/normR3D, 
                             (P0+P1+P2+P3+P4)/normR3D, 
                             (P0+P1+P2+P3+P4+P5)/normR3D, 
                             1.0f};
	vector<float> normRi2D = {P0/normR2D, 
                             (P0+P1)/normR2D, 
                             (P0+P1+P2)/normR2D, 
                             (P0+P1+P2+P3)/normR2D, 
                             1.0f};

//	printf("R1=%.5f, R2=%.5f, R3=%.5f, R4=%.5f, R5=%.5f, R6=%.5f\n",normRi[1],normRi[2],normRi[3],normRi[4],normRi[5],normRi[6]);
        
    // Generate a random number
	float randmov   = drandom(generator);//drandom(mxmx);//drandom(generator);
/*
	if(randmov <= normRi[0]){// Perform a random movement if free spot exists
		neigh = returnEmptyPlace(casp,indx);
		//printf("Random movement\n");
	}
	else{// Find the direction of max gradient and move if free spot
		int maxGradInd = max_element(Pi.begin(),Pi.end()) - Pi.begin();
		if(!casp.lattice[indx+casp.indcNeigh[maxGradInd]]){
			neigh = indx+casp.indcNeigh[maxGradInd];
			//printf("Deterministic movement\n");
		}
		else{// Find the direction of the next max gradient and move if free spot
			for(int i=1;i<6;i++){
				Pi[maxGradInd]=0;
                maxGradInd = max_element(Pi.begin(),Pi.end()) - Pi.begin();
                if(!casp.lattice[indx+casp.indcNeigh[maxGradInd]]){
                   	neigh = indx+casp.indcNeigh[maxGradInd];
					break;
                }
			}
		}
	}
*/

    // Search for available free spots with respect to moving probabilities
	if(indx <= 6*zsize && casp.phenchange <= 0){
		if(!casp.lattice[indx+casp.indcNeigh[0]] && randmov>=normRi2D[0] && randmov<normRi2D[1])
                neigh = indx+casp.indcNeigh[0]; //printf("Back \n");
       	 	
        else if(!casp.lattice[indx+casp.indcNeigh[1]] && randmov>=normRi2D[1] && randmov<normRi2D[2])
                neigh = indx+casp.indcNeigh[1]; //printf("Front \n");
      
        else if(!casp.lattice[indx+casp.indcNeigh[2]] && randmov>=normRi2D[2] && randmov<normRi2D[3])
               	neigh = indx+casp.indcNeigh[2]; //printf("Left \n");

       	else if(!casp.lattice[indx+casp.indcNeigh[3]] && randmov>=normRi2D[3] && randmov<normRi2D[4])
               	neigh = indx+casp.indcNeigh[3]; //printf("Right \n");
	}
	else{
        if(!casp.lattice[indx+casp.indcNeigh[0]] && randmov>=normRi3D[0] && randmov<normRi3D[1])
               	neigh = indx+casp.indcNeigh[0]; //printf("Back \n");
        
        else if(!casp.lattice[indx+casp.indcNeigh[1]] && randmov>=normRi3D[1] && randmov<normRi3D[2])
               	neigh = indx+casp.indcNeigh[1]; //printf("Front \n");
        
        else if(!casp.lattice[indx+casp.indcNeigh[2]] && randmov>=normRi3D[2] && randmov<normRi3D[3])
               	neigh = indx+casp.indcNeigh[2]; //printf("Left \n");
        	
        else if(!casp.lattice[indx+casp.indcNeigh[3]] && randmov>=normRi3D[3] && randmov<normRi3D[4])                	neigh = indx+casp.indcNeigh[3]; //printf("Right \n");
        
        else if(!casp.lattice[indx+casp.indcNeigh[4]] && randmov>=normRi3D[4] && randmov<=normRi3D[5])                	neigh = indx+casp.indcNeigh[4]; //printf("Down \n");
        
        else if(!casp.lattice[indx+casp.indcNeigh[5]] && randmov>=normRi3D[5] && randmov<=normRi3D[6])                	neigh = indx+casp.indcNeigh[5]; //printf("Up \n");
	}
	return neigh;

}//end of returnPosChance






void caStep(vector<Cell>& cells, vector<Cell>& cellsTmp, CAArray& casp, DataArray& host, InitialData& id)
{

    int Ae, nPhen, nnPhen, newSite, newPos;
	double d_prob; //spontaneous death probability
	Cell  currCell, newCell;
    random_shuffle(cells.begin(), cells.end()); //shuffling cells

  	while (!cells.empty()) {

		currCell = cells.back()  ; //pick the cell
        cells.pop_back()         ;

		if(casp.delaydiv <= 0)	currCell.age ++  ; // increase age of cells
		

		d_prob = (double)rand()/(double)RAND_MAX;

	   	if(d_prob <= casp.spdeath)
	   	{
			currCell.is_alive = false   ;
			casp.lattice[currCell.place]=false;
//          ellsTmp.push_back(currCell);

	   	}
		else
		{
           	Ae    = returnNeigh(casp,currCell.place); //return #neighbors Ae>=Ai for migration
			nPhen = returnNeighPhen(casp,currCell.place); //return phenotype of neighbours
 	
			//divide
           	if(currCell.age >= casp.agediv && currCell.is_alive && currCell.ndiv <= casp.divpot){
			//printf("Dividing...\n");
				if(currCell.place <= 6*zsize)//check if is close to the bottom
					newSite = returnEmptyPlace2D(casp,currCell.place); //check for new site
				
				else
					newSite = returnEmptyPlace(casp,currCell.place); //check for new site
				
				if (newSite) {//if there is a new spot

	      			newCell           = currCell ;
               	    newCell.place     = newSite  ;
             		newCell.age       = 0        ;
//             		newCell.is_alive  = true     ;
               		currCell.age      = 0        ;
//             		currCell.is_alive = true     ;
					currCell.ndiv ++ 	     ;
					newCell.ndiv ++		     ;
					casp.lattice[newSite]  = true;

					nnPhen = returnNeighPhen(casp,newCell.place); //return phenotype of neighbours
					if((newSite <= 6*zsize || nnPhen > 0) && casp.phenchange <= 0){
						casp.phen[newCell.place]  = true;
						casp.phen[currCell.place] = true;
					}
               		    cellsTmp.push_back(currCell);
               			cellsTmp.push_back(newCell );

           			}else{//no free spot

               			currCell.is_alive = false   ;
               			cellsTmp.push_back(currCell);
          	     	}	
           		}else if(Ae>=casp.adhes && currCell.is_alive && casp.phenchange > 0){
				    //printf("Migrating...\n");
                    newPos = returnPosChance(casp,host,id,currCell.place);
                    if(newPos){
                        newCell       = currCell;
                        newCell.place = newPos  ;
                        casp.lattice[currCell.place]=false;
                        casp.lattice[newPos        ]=true ;
					    cellsTmp.push_back(newCell);
				    }else  //no migration - do nothing
                        cellsTmp.push_back(currCell);
            
			   }else if(Ae>=casp.adhes && currCell.is_alive && casp.phenchange <=0 && !casp.phen[currCell.place] && nPhen == 0){ //migrate
				//printf("Migrating...\n");
				    newPos = returnPosChance(casp,host,id,currCell.place);
      				if(newPos){

 					    newCell       = currCell;
       				    newCell.place = newPos  ;
					    casp.lattice[currCell.place]=false;
                        casp.lattice[newPos        ]=true ;
					    nnPhen = returnNeighPhen(casp,newCell.place); //return phenotype of neighbours
					    if(newPos <= 6*zsize || nnPhen > 0)
						    casp.phen[newCell.place]  = true;
					    
               			cellsTmp.push_back(newCell);

             	    } else  //no migration - do nothing
                        cellsTmp.push_back(currCell);
             			
	   		    }else   
                    cellsTmp.push_back(currCell);  
		    }

       	}
       	cells.swap(cellsTmp);
    
}//end of ca_step



void ExportCellCount(FILE* fna, InitialData& id, vector<Cell>& cells, int i){

	int size = cells.size();
   	fprintf(fna, "%lf, %d\n", i*id.dt_imp, size);
  	fflush(fna);
}



void ExportCA (FILE* fca, InitialData& id, vector<Cell>& cells, int i){

	
	int size = cells.size();
	int j;	
	cout << "Elapsed time: " << (int)ceil(i*id.dt_imp) << ", #cells: " << size << endl;
	for(j=0; j<size; j++)
	{ 
		fprintf(fca, "%d, %d, %d, %d \n", (int)ceil(i*id.dt_imp), cells[j].place, cells[j].age, cells[j].is_alive);
       	fflush(fca);
   	}
	
}//end




