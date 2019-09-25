//Fract8_Con: Visualize the connectivity as a function of phi and rho;
//Fract8: 1) Apply triangular lattice, 2)use rectangular cell, h=1, w=sqrt(3)/2; 3) 35 by 201 biofilm; 
//Fract7* .3: repeat 100 times, only print the data for S-Rg 
//Fract7* .2: repeat 1000 times, don't generate the data for S-Rg
//Fract7: output the largest cluster size in each frame. Square lattice;
//Fract6* _2: Always visualize the S-Rg despite the central cell is not in the largest cluster.
//Fract6: _1: Vertical Lineage Inheritance;
//Fract5: Fix the bug: Considering the cell shape in visualizing the gyration radius.
//Fract4://Fract4: Visulize the Gyration radius and S to calculate the correlation function.
//Fract2: Use the MT random number generator. Biofilm filter: give up these whose phi is deviated more than 0.02%
//Fract1: update the algorithm to calculate the gyration radius.
//Fract_S_Rg_Phi: fractal dimension is visulized form the cluster size as a function of gyration radius. 
//Fractal dimension D: M(L)~L^D. 
//_n6: define anchor size as the part, which is inside the convex hull of the off-cell cluster which has the most cells inside the convex hull.
//_n5: consider all the off-cell clusters. Find out the one which is most likely to be the anchor (The cluster has the most off-cells inside the anchor). Use the cluster size of this off-cell cluster as the anchor size n5.
//n4: exclude all the clusters which are less than 200.
//_n3: calculate both n1 and n2
//Spi7 _n1: anchor size the the average of all te cluster sizes. Each cluster has the same probability to be picked up.
//Spi7: Visulize the structural anchor size. 
//Spi6*b5: period and convex hull, time duration: 50->TimLen.
//Spi6*b5: Convex Hull: typo is verified. consider wave from LocSta->TimLen. if the cells with most peaks are less than 10, then return failure. Anchor size is 0. 
//Spi6*b4: period and convex hull: find out all the cells with the most peaks. (there are too many cells if use the cells with peraks of MCP-1)
//spi6*b3: 1. updated method to categorize transmission: use the ratio of firing cells before the wave reaches the edge to tell break up. 
//         2. period and convex hull: find out all the cells with the most peaks (MCP) and peaks of MCP-1.  
//Spi6*b2: use the area of the convex hull as the anchor size.
//Spi6*b1: fix the bug in ConvexHull: j>2 in finding the convex hull.
//Spi61_b2: anchor is defined as the sum of cluster size containing any off cells inside the Convex Hull
//Spi6 #1. New method to calculate the period: 20 cells, average the 10 cells whose periods are less deviated from the mean; #2. anchor size is defined as the off cells inside the convex hull.
//Spi5. 1. repeat 500 times for each parameter combination; When calculating the period, consider only the wave between 100 to 300.
//Spi4.7. Test the relationsip between period and phi. Calculate the period of each simulaton.
//Spi4.2. Use the parameters of the distribution as the input variables.
//Spi4.1. Select 64 sampling cells uniformly in the biofilm.
//Spi4. Use mean and variance as the inputs indtead of parameters of the distribution. Fix a bug in deciding the wave pattern: 
        //if phi is too small to find out 1 capable cell, the code will just use non-capable cell.
//Spi3.1. New algorithm to find brek up: N_firing/N_sample< 0.8;
//Spi3: antomatically catigrize the wave style: 0 died; 2, spiral wave+ bereak up; 3, break up, 5. normal; 8. secondary+Spiral; 9. secondary; 
//binary noise
//Spi2: radial lineage inheritance for any distributions
//_Lognor: use \nu and \sigma as the input arguments of lognormal distribution
//_Gau: parameter, for example tau is log normal distributed
//Spi ->Spi1: Eliminate parameter v_0;
//SpiWave_SS_WT.c is from Fit17Wave_SS_WT_spiral.c. Use square unit, rectangular lattice, reflective boundary in top and both sides, obsorbing in bottom. 
//Fit17: updated method to visulize wavelength.
//_SS: same tructure. In one repeating, all combinations of (a,b) will use the same structure.
//Fit16: replace b to tau. tau=1/b;
//FitClu15: dt_write=0.5; only consider these periods from firing capable cells; T_th=0.1;
//FitClu11->FitClu14: boundary condition changed. FitClu11: 4 absorbing boundaries (dxdt_Tri); FitClu14: obsorbing in the bottom, 3 reflective boundaries in the other three sides (dxdt_Tri_Top).
//FitClu11->FitClu14: Biofilm is by 200, which is located inside a 100 by 200 biofilm. row 0 is triggered, row 1 to row 35 is the observed area. The left rows is used to reduce the influence of bottom obsorbing boundary.
//11Loop_a_1: explore the parameter space of F-N model to explain the trkA, wild type and ktrA.
//FitClu10->11: Directionality or transmission speed as a function of phi; Trigger in the middle: row and column 48 49 50 51 52  ## Delete neighbor in SearchTree_G4
//FitClu9->10: introduce parameter c: 
//FitClu8->9: multi-generation inheritance to distribute the parameter.
//FitClu7->FitClu8: Triangle lattice takes the place of square lattice.
//FitClu6->FitClu7Introduce noise which follows Gaussian distribution or uniform distribution in log space.
//Introduce the parameter gamma: {a,b,gamma}: wt=b*(v-gamma*w);
//Attention: M, the number of columns; N, the number of rows. Size of matrix: N by M; 
//Delete the limit that v>=0;
//Change the model: When v decreased less than a/4, make w=0;
//The bug1: CluGen: Initiation of cluster,neighbor solved; 
//The original code, for all the c files in folder "Fitzhugh"
//aNoise contains the distribution of parameter a;

#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "mersenne_twister.hpp"


#define RAND_MAX 2147483647
#define PI 3.1415926
#define EPSILON 0.0000000001
#define EPSILON2 0.00001


//Use the equation dEdt=-delta(-D*deltaE);
//Add noise to explain the double spiking wich  occurses in the experiment of Gurol. Applied to 2D grid. Parallel the Runge Kunta part within while loop.
// Parallel Algorithm, run in cluster: No iteractive jobs.
//Trigger 9 cells in the ceter of cell community, use Runge Kutta Algorithm.
// Reproducing the work in Fig 3d and 3e in Pindle et al's Nature, 2015,"Ion channels enable electrical communication in bacterial communities" 
//May 11, 2016, Purdue; 
//2d N*M cells,Protocol b: Equilibrate for 100 min, then pulse then pulse S=100 uM ?or E=200*1000 uM?
FILE * name(int N,int M,int *coefx, int *coefy,char *filename);
void Island(double *aNoise,double percentage,int N,int M,double lisland,double F);
int InOut(int target,int *source);
int CluGen(double *aNoise,double phi,double psi, int N,int M,double aLow,double DHigh);
int CluGen_G5_1(double *nvth,int height,int width,double phi,double rho,double DLow,double DHigh,
              double alpha, double beta, double a1, double a2, double a3);
int CluRadLin(double *nvth, int height, int width, char mode, double rho, double phi, double DLow, 
             double median, double lambda);
int mothergrandmother(char mode, double *nvth, int HeiMat, int WidMat, int height, double DLow, double DHigh,
           double *possibility, double ponoffoff, double ponoffon, double pononoff, double pononon);
int SearchTree_G4(double *nvth,int *state,int i,int j,int height,int width,int HeiMat,int WidMat,
                  double DHigh,int *cluster);
int SearchClu_G4_Squ(double *nvth,int *state,int i,int j,int height,int width,int HeiMat,int WidMat,
                  double DHigh,int *cluster);
int SearchClu_G4_Squ_Frame(double *nvth,int *state,int i,int j, int RowSta, int RowEnd, int ColSta, int ColEnd, 
       int HeiMat,int WidMat, double DHigh,int *cluster);
int SearchClu_G4_Tri_Frame(double *nvth,int *state,int i,int j,int RowSta, int RowEnd, int ColSta, int ColEnd, 
       int HeiMat,int WidMat, double DHigh,int *cluster);
int SearchCon_G4_Tri_Frame(double *nvth,int *state,int i,int j,int RowSta, int RowEnd, int ColSta, int ColEnd, 
       int HeiMat,int WidMat, double DHigh,int *cluster);
int LargestCluster(double *nvth, int HeiMat, int WidMat, int height,int width, int row, double DLow, double DHigh, int square);

void dxdt_Squ(double *EVnST,double *deri,double *aNoise,double *tauNoise,double *cNoise, double gamma, int coefx,int coefy,int N,int M,int Leng, double deltax, double deltay);
void dxdt_Squ_o4(double *EVnST,double *deri,double *aNoise,double *tauNoise,double *cNoise, double gamma, int coefx,int coefy,int N,int M,int Leng, double deltax, double deltay);
void dxdt_Tri(double *EVnST,double *deri,double *aNoise,double *tauNoise, double *cNoise, double gamma, int coefx,int coefy,int N,int M,int Leng, double deltax, double deltay);   //deri[] is used to store the data of dVdt,dndt......
void dxdt_Tri_Top(double *EVnST,double *deri,double *aNoise,double *tauNoise,double *cNoise, double gamma, int coefx,int coefy,int N,int M,int Leng, double deltax, double deltay);
void print1(int N,int M,int Length, double t, double tInterval,double deltat,double *EVnST, FILE *fpt);
void print1_1(int N,int M,int Length, double t, double tInterval,double deltat,double *EVnST, FILE *fpt);
int print5(int N, int M, double t, double *EVnST, FILE *fpt);
void Mplus(int n, int m, int Length, double deltat, double *Matin, double *Matout, double *Deri);     // addition of N*M matrix;
void Transport(int N,int M,int triy,int trix,int Ln,int Lm,int n,double *source,double *Tprint);

int MeaLam(double *wave, int height, int width, int TimStep, double threshold); //propagation length
int wavefront(double *wave, int height, int width, int ColSta, int ColEnd, double threshold);
int waveback(double *wave, int height, int width, int ColSta, int ColEnd, double threshold);
double wavelength(double *wave, int height, int width, int LenSec, double threshold);
int NFir(double *wave, int height, int width, double threshold);

int FindPeaks(double deltat, double TimLen, double dt_write, double *trace, int LocSta, int LocEnd, const int TimStep, double v_th, double *PeakTim, int LenPeakTime);
int WheDied(double * wave, int height, int width, int TimStep, double v_th);
int WheDied_Tim(double *wave, int height, int width, int TimStep, double v_th);   //died, -1; non-died, return the time sequence.
int WheNeverDied(double *wave, int height, int width, int TimStep, double v_th);

struct ConHulPoint{
	double xaxis;
	double yaxis;
	double cos;
	};
int ConvexHull(double *SetX, int SizeSetX, int length, double *ConHull, int LenConHull);	
int ExchangeD(double * a, double *b);
int WheInConHull(double *ConHull, int NumPoints, int LenConHull, double x, double y);
double CroPro2D(double x1, double y1, double x2, double y2);
double PerCelMaxPea(double *PerCMP, int NumCMP);

double TMean(double *Tij, int height, int width, double threshold);
int MinArrXZ_I(int *array, int length);
int MaxArrXZ_I(int *array, int length);
double SumMat_D(double *source, int height, int width, int HeiMat, int WidMat, double threshold, char mode);
void IniMat_D(int n,int m, double constant,double *matrix);
void IniMat_I(int n,int m, int constant,int *matrix);
int PriMatD(double *source,int height,int width,int HeiMat,int WidMat,FILE *fpt);
int PriMatI(int *source,int height,int width,int HeiMat,int WidMat,FILE *fpt);
int WheIn(int target,int *source, int length);

double minimum(double a,double b);
int equal(double a, double b, double error);
int WheMultiple(double numerator, double denominator, double error);

double AssignPhi(double phi,double low,double high);
double AssignPhi_MT(double phi,double low,double high);
double DisUni(double mean, double HalWid);
double gaussrand_K(void);
double GyrationRadius(int *array, int Height, int Width, int length); //input the array of locations (location= i*Width+j), return the gyration radius
double GyrationRadius_Tri(int *array, int Height, int Width, double dy, double dx, int length);

MTRand RanGen( (unsigned)time(NULL)); //Generate an instance from object MTRand, initiate it with the current time tt.


int main(void){ 

int N=35,M=230;   //N, the number of rows, M, the number of columns. L: the width of the frame.
//double dy=1, dx=sqrt(3.0)/2; //The size of single cell, used to visualize the gyration radius.
//$$ 1. change the value of parameter a and b,c;                                             
double tau=1, tau_off=0, phi, rho=0, phi_deviation=1.0/M; //when rho=1, phi_deviation need to be 1.0/M otherwise process will be frozen.

//int height, width, heights[]={N},widths[]={M};
double phis[126], rhos[101];
for (int i=0; i<126; i++){ //0, 4/1000, 0.5
	phis[i]=i*4.0 /1000.0;}
	
for (int i=0; i<101; i++){ //0, 1/100, 1;
	rhos[i]=i* 1.0 /100.0;}


int  Len_phis=sizeof(phis)/sizeof(phis[0]), Len_rhos=sizeof(rhos)/sizeof(rhos[0]), NumRep=1000;
int i,j, iphi, irho, ip, iheight; //i, j is only used for any short, self-closed loop. 
double tauNoise[N*M];   //row 0 stores the v, row 1 stores w.

FILE *fptFrac, *fptMap, *fptCon;     //fpt: ThT; fpt2: parameter;
char filename[100], filename2[100];
char  num[100]; 

time_t tt=time(NULL);


//srand((unsigned) time(&tt));  //Initiating random number generator only once. 

//MTRand RNG( (unsigned)tt);
//double x1=RNG.randDblExc(); //uniform random number: (0,1)
//double r1= RNG.rand();      //Uniform random number: [0, 1]

for(irho=0; irho< Len_rhos; irho++){
	rho=rhos[irho];
	strcpy(filename, "Fra8_Con_Phi");
	strcat(filename,"Rho");   sprintf(num, "%d",(int)(rho*100+EPSILON) );strcat(filename,num);
	strcat(filename, "_TH35W230_Fig3B.dat");	
    fptCon=fopen(filename, "w");     
    double connection[Len_phis];    
    IniMat_D(1, Len_phis, 0, &connection[0] );  //why do we initiate? The initiating value may be printed into file if something anormal happens.
	
for(iphi=0; iphi<Len_phis; iphi++){ 	
	phi=phis[iphi];  
	
	
    strcpy(filename, "Fra8");     
    strcat(filename,"Phi");   sprintf(num, "%d",(int)(phi*10000+EPSILON) );strcat(filename,num);
    strcat(filename,"Rho");   sprintf(num, "%d",(int)(rho*100+EPSILON));strcat(filename,num); 
    strcpy(filename2, filename);
    strcat(filename,"TH35W230_map.dat");	strcat(filename2, "TH35W230_SRg.dat");
	//fptMap=fopen(filename,"w");
	//fptFrac=fopen(filename2, "w");


for(ip=0;ip<NumRep; ip++){         
  //Initiate the biofilm
  //CluGen_G5_1(double *nvth,int height,int width,double phi,double rho,double DLow,double DHigh,
              //double alpha, double beta, double a1, double a2, double a3)
  double phi_actual;
  do{
	phi_actual=0;
	CluGen_G5_1(tauNoise, N, M, phi, rho, tau_off, tau, 1.0, 1.0, 0.0, 1.0, 0.0);
	//for(i=0; i<N*M; i++){	 
	      //if( RanGen.randExc()< phi ){
	          //tauNoise[i]=tau;
	       //}else
	          //tauNoise[i]=tau_off;
	 //}	
	for(i=0; i<N*M; i++){
		phi_actual+= (tauNoise[i]== tau);
		}  	  
   } while ( fabs(phi_actual/N/M- phi)> phi_deviation );

  int state[N*M], cluster[N*M], CluSiz, EffClu[2][N*M], End_cluster, End_EffClu;
  //*************************************************************
  //int SearchCon_G4_Tri_Frame(double *nvth,int *state,int i,int j,int RowSta, int RowEnd, int ColSta, int ColEnd, 
  //     int HeiMat,int WidMat, double DHigh,int *cluster)
  
  IniMat_I(1, N*M, 0, state);  //0 means unchecked. 1: checked, not target; 2: checked, target.
  int stop=0;
  for (j=0; (j< M) && (stop==0); j++) {
      if ( SearchCon_G4_Tri_Frame(tauNoise, state, 0, j, 0, N, 0, M, N, M, tau, &cluster[0]) == 1 ) {  
	    connection[iphi]+=1;
		stop=1;
	  } 
  }	

}   //end of loop ip  

  connection[iphi]/= NumRep; 
  printf("rho %lf tau=%lf phi %lf\n",rho, tau, phi);  
}  //end of loop iphi

PriMatD( &connection[0], 1, Len_phis, 1, Len_phis, fptCon); 		
fclose(fptCon);
//fclose(fptMass);
}  //end of loop irho

printf("running this tourine costs %ld seconds of time\n",time(NULL)-tt);
return 0;	
}

//Reflective boundary in the top and two sides, observing in the bottom.
void dxdt_Squ(double *EVnST,double *deri,double *aNoise,double *tauNoise,double *cNoise, double gamma, int coefx,int coefy,int N,int M,int Leng, double deltax, double deltay)
 {  
 int i=0;
 #pragma omp parallel
{
#pragma omp for
 for(i=0;i<N*M;i++){  
  deri[Leng+i]=1/tauNoise[i]*(EVnST[i]-gamma*EVnST[Leng+i]);                       //dwdt;   
   deri[i]=cNoise[i]*(  EVnST[i]*(1-EVnST[i])*(EVnST[i]-aNoise[i])-EVnST[Leng+i] );   //dvdt: This step also initiates the deri matrix.
 
    
  
  if(i/M==0)
	 deri[i]+=(EVnST[i+M]-EVnST[i])/pow(deltay,2);
  else if(i/M==N-1)	  
	 deri[i]+=(0-EVnST[i])/pow(deltay,2)+(EVnST[i-M]-EVnST[i])/pow(deltay,2);
  else
     deri[i]+=(EVnST[i+M]-EVnST[i])/pow(deltay,2)+(EVnST[i-M]-EVnST[i])/pow(deltay,2);
     
   if(i%M==0)
      deri[i]+=(EVnST[i+1]-EVnST[i])/pow(deltax,2);
   else if(i%M==M-1)
      deri[i]+=(EVnST[i-1]-EVnST[i])/pow(deltax,2); 
   else
      deri[i]+=(EVnST[i+1]-EVnST[i])/pow(deltax,2)+(EVnST[i-1]-EVnST[i])/pow(deltax,2);  	     
 }     //end of for loop
 } // end of parallel
}

void dxdt_Squ_o4(double *EVnST,double *deri,double *aNoise,double *tauNoise,double *cNoise, double gamma, int coefx,int coefy,int N,int M,int Leng, double deltax, double deltay)
 {  
 int i=0;
 #pragma omp parallel
{
#pragma omp for
 for(i=0;i<N*M;i++){  
  deri[Leng+i]=1/tauNoise[i]*(EVnST[i]-gamma*EVnST[Leng+i]);                       //dwdt;   
   deri[i]=cNoise[i]*(  EVnST[i]*(1-EVnST[i])*(EVnST[i]-aNoise[i])-EVnST[Leng+i] );   //dvdt: This step also initiates the deri matrix.
 
    
  
  if(i/M==0)
	 deri[i]+=(0-EVnST[i])/pow(deltay,2)+(EVnST[i+M]-EVnST[i])/pow(deltay,2);
  else if(i/M==N-1)	  
	 deri[i]+=(0-EVnST[i])/pow(deltay,2)+(EVnST[i-M]-EVnST[i])/pow(deltay,2);
  else
     deri[i]+=(EVnST[i+M]-EVnST[i])/pow(deltay,2)+(EVnST[i-M]-EVnST[i])/pow(deltay,2);
     
   if(i%M==0)
      deri[i]+=(0-EVnST[i])/pow(deltax,2)+(EVnST[i+1]-EVnST[i])/pow(deltax,2);
   else if(i%M==M-1)
      deri[i]+=(0-EVnST[i])/pow(deltax,2)+(EVnST[i-1]-EVnST[i])/pow(deltax,2); 
   else
      deri[i]+=(EVnST[i+1]-EVnST[i])/pow(deltax,2)+(EVnST[i-1]-EVnST[i])/pow(deltax,2);  	     
 }     //end of for loop
 } // end of parallel
}

void dxdt_Tri(double *EVnST,double *deri,double *aNoise,double *tauNoise,double *cNoise, double gamma, int coefx,int coefy,int N,int M,int Leng, double deltax, double deltay)
 {  
 int i=0;
 #pragma omp parallel
{
#pragma omp for
for(i=0;i<N*M;i++){  
   deri[Leng+i]=1/tauNoise[i]*(EVnST[i]-gamma*EVnST[Leng+i]);                      //dwdt;
   
   deri[i]=cNoise[i]*(EVnST[i]*(1-EVnST[i])*(EVnST[i]-aNoise[i])-EVnST[Leng+i]);   //dvdt: This step also initiates the deri matrix.
   
   if(i/M==0)    //v_yy
	 deri[i]+=(EVnST[i+M]-2*EVnST[i])/pow(deltay,2);
   else if(i/M==N-1)	  
	 deri[i]+=(0-EVnST[i])/pow(deltay,2)+(EVnST[i-M]-EVnST[i])/pow(deltay,2);
   else
     deri[i]+=(EVnST[i+M]-EVnST[i])/pow(deltay,2)+(EVnST[i-M]-EVnST[i])/pow(deltay,2); 	 
    
  if((i%M)%2==0){ 
	       
   if(i/M==N-1&&i%M==0)
      deri[i]+=(EVnST[i+1]-EVnST[i])/pow(deltax,2)/2+(0-EVnST[i])/pow(deltax,2)/2+(0-EVnST[i])/pow(deltax,2)/2+(0-EVnST[i])/pow(deltax,2)/2;
   else if(i/M==N-1)
      deri[i]+=(EVnST[i+1]-EVnST[i])/pow(deltax,2)/2+(EVnST[i-1]-EVnST[i])/pow(deltax,2)/2+(0-EVnST[i])/pow(deltax,2)/2+(0-EVnST[i])/pow(deltax,2)/2;
   else if(i%M==0)
      deri[i]+=(EVnST[i+1]-EVnST[i])/pow(deltax,2)/2+(0-EVnST[i])/pow(deltax,2)/2+(EVnST[i+M+1]-EVnST[i])/pow(deltax,2)/2+(0-EVnST[i])/pow(deltax,2)/2;
   else
      deri[i]+=(EVnST[i+1]-EVnST[i])/pow(deltax,2)/2+(EVnST[i-1]-EVnST[i])/pow(deltax,2)/2+(EVnST[i+M+1]-EVnST[i])/pow(deltax,2)/2+(EVnST[i+M-1]-EVnST[i])/pow(deltax,2)/2;
            }
            
            
   else if((i%M)%2==1){
	       
   if(i/M==0&&i%M==M-1)
      deri[i]+=(0-EVnST[i])/pow(deltax,2)/2+(EVnST[i-1]-EVnST[i])/pow(deltax,2)/2+(0-EVnST[i])/pow(deltax,2)/2+(0-EVnST[i])/pow(deltax,2)/2;
   else if(i/M==0)
      deri[i]+=(EVnST[i+1]-EVnST[i])/pow(deltax,2)/2+(EVnST[i-1]-EVnST[i])/pow(deltax,2)/2+(0-EVnST[i])/pow(deltax,2)/2+(0-EVnST[i])/pow(deltax,2)/2;
   else if(i%M==M-1)
      deri[i]+=(0-EVnST[i])/pow(deltax,2)/2+(EVnST[i-1]-EVnST[i])/pow(deltax,2)/2+(0-EVnST[i])/pow(deltax,2)/2+(EVnST[i-M-1]-EVnST[i])/pow(deltax,2)/2;
   else
      deri[i]+=(EVnST[i+1]-EVnST[i])/pow(deltax,2)/2+(EVnST[i-1]-EVnST[i])/pow(deltax,2)/2+(EVnST[i-M+1]-EVnST[i])/pow(deltax,2)/2+(EVnST[i-M-1]-EVnST[i])/pow(deltax,2)/2;
            }     
}     //end of for loop
 } // end of parallel
}

void dxdt_Tri_Top(double *EVnST,double *deri,double *aNoise,double *tauNoise,double *cNoise, double gamma, int coefx,int coefy,int N,int M,int Leng, double deltax, double deltay)
 {  
 int i=0;
 #pragma omp parallel
{
#pragma omp for
for(i=0;i<N*M;i++){     
   deri[Leng+i]=1/tauNoise[i]*(EVnST[i]-gamma*EVnST[Leng+i]);                      //dwdt;
   
   deri[i]=cNoise[i]*(EVnST[i]*(1-EVnST[i])*(EVnST[i]-aNoise[i])-EVnST[Leng+i]);   //dvdt: This step also initiates the deri matrix.
   
   if(i/M==0)
	 deri[i]+=(EVnST[i+M]-EVnST[i])/pow(deltay,2);   //Reflective on the top
   else if(i/M==N-1)	  
	 deri[i]+=(0-EVnST[i])/pow(deltay,2)+(EVnST[i-M]-EVnST[i])/pow(deltay,2);  //Obsorbing on the bottom
   else
     deri[i]+=(EVnST[i+M]-EVnST[i])/pow(deltay,2)+(EVnST[i-M]-EVnST[i])/pow(deltay,2); 	 
    
  if((i%M)%2==0){ 
	       
   if(i/M==N-1&&i%M==0)   //Last row, first column
      deri[i]+=(EVnST[i+1]-EVnST[i])/pow(deltax,2)/2+(0-EVnST[i])/pow(deltax,2)/2+(0-EVnST[i])/pow(deltax,2)/2;
   else if(i/M==N-1)  //Last row
      deri[i]+=(EVnST[i+1]-EVnST[i])/pow(deltax,2)/2+(EVnST[i-1]-EVnST[i])/pow(deltax,2)/2+(0-EVnST[i])/pow(deltax,2)/2+(0-EVnST[i])/pow(deltax,2)/2;
   else if(i%M==0)
      deri[i]+=(EVnST[i+1]-EVnST[i])/pow(deltax,2)/2+(EVnST[i+M+1]-EVnST[i])/pow(deltax,2)/2;
   else
      deri[i]+=(EVnST[i+1]-EVnST[i])/pow(deltax,2)/2+(EVnST[i-1]-EVnST[i])/pow(deltax,2)/2+(EVnST[i+M+1]-EVnST[i])/pow(deltax,2)/2+(EVnST[i+M-1]-EVnST[i])/pow(deltax,2)/2;
            }
            
            
   else if((i%M)%2==1){  //Odd columns
	       
   if(i/M==0&&i%M==M-1)   //First row, last column
      deri[i]+=(EVnST[i-1]-EVnST[i])/pow(deltax,2)/2;
   else if(i/M==0)
      deri[i]+=(EVnST[i+1]-EVnST[i])/pow(deltax,2)/2+(EVnST[i-1]-EVnST[i])/pow(deltax,2)/2;
   else if(i%M==M-1)
      deri[i]+=(EVnST[i-1]-EVnST[i])/pow(deltax,2)/2+(EVnST[i-M-1]-EVnST[i])/pow(deltax,2)/2;
   else
      deri[i]+=(EVnST[i+1]-EVnST[i])/pow(deltax,2)/2+(EVnST[i-1]-EVnST[i])/pow(deltax,2)/2+(EVnST[i-M+1]-EVnST[i])/pow(deltax,2)/2+(EVnST[i-M-1]-EVnST[i])/pow(deltax,2)/2;
            }     
}     //end of for loop
 } // end of parallel
}

double minimum(double a,double b)
{if (a<=b)
	return a;
else 
	return b;	
}


void print1(int N,int M,int Length, double t, double tInterval,double deltat,double *EVnST, FILE *fpt)
{         //print the data of every cell into file.
   
    int j=0; 	      
    	      
    if(WheMultiple(t,tInterval,deltat/10)|| equal(t,0.5,deltat/10)){   
        	for(j=0;j<N*M;j++){
        		if(j%M==M-1)
        		   fprintf(fpt,"%f\n",*(EVnST+j));	
        		else
		           fprintf(fpt,"%f  ",*(EVnST+j));
	        }	          
      	
        } 
         
}  

void print1_1(int N,int M,int Length, double t, double tInterval,double deltat,double *EVnST, FILE *fpt)
{         //print the data of every cell into file.
   
    int j=0; 	      
    	      
    if(WheMultiple(t,tInterval,deltat/10)|| equal(t,2*deltat,deltat/10)){   
        	for(j=0;j<N*M;j++){
        		if(j%M==M-1)
        		   fprintf(fpt,"%f\n",*(EVnST+j));	
        		else
		           fprintf(fpt,"%f  ",*(EVnST+j));
	        }	          
      	
        } 
         
} 

int print5(int N, int M, double t, double *EVnST, FILE *fpt)   //print row 1:5:201;
{  int i=0,j=0;
	if(N!=1&&M!=1&&fabs(fmod(t-0.00005,0.5)-0.5)<0.00008){   
        	for(i=0;i<N;i=i+5){
			for(j=0;j<M;j=j+1){
				if(j==M-1)
        		   fprintf(fpt,"%lf\n",*(EVnST+i*M+j));	
        		else
		           fprintf(fpt,"%lf  ",*(EVnST+i*M+j));
	          }}          
    }
return 1;
}

//Define the addition of Matrix.
void Mplus(int n, int m, int Length, double deltat, double *Matin, double *Matout, double *Deri)
{    // addition of N*M matrix;
   int i,j;
   for(i=0;i<n;i++)
       for(j=0;j<m;j++)
           *(Matout+i*Length+j)=*(Matin+i*Length+j)+*(Deri+i*Length+j)*deltat;
    }
//function Transport: write the data in aNoise or EVnST into Tprint
void Transport(int N,int M,int triy,int trix,int Ln,int Lm,int n,double *source,double *Tprint)
{
	int i=0,j=0,m=0;
	for(i=(N-1)/2-(Ln-1)/2;i<=(N-1)/2+(Ln-1)/2;i++)      //the triy  by trix grid in the left side.  Write the data to the Tprint(m,n); source==&EVnST[4][0]or &aNoise[0][0]
		for(j=(M-1)/4-(Lm-1)/2;j<=(M-1)/4+(Lm-1)/2;j++)
		  *(Tprint+400*(m++)+n)=*(source+i*M+j);	
      for(i=(N-1)/2-(Ln-1)/2;i<=(N-1)/2+(Ln-1)/2;i++)      //the triy  by trix grid in the right side.
		for(j=3*(M-1)/4-(Lm-1)/2;j<=3*(M-1)/4+(Lm-1)/2;j++)
		  *(Tprint+400*(m++)+n)=*(source+i*M+j);
		for(i=(N-1)/4-(Ln-1)/2;i<=(N-1)/4+(Ln-1)/2;i++)      //the triy  by trix grid in the top side.
	       for(j=(M-1)/2-(Lm-1)/2;j<=(M-1)/2+(Lm-1)/2;j++)
		  *(Tprint+400*(m++)+n)=*(source+i*M+j);
	  	for(i=3*(N-1)/4-(Ln-1)/2;i<=3*(N-1)/4+(Ln-1)/2;i++)      //the triy  by trix grid in the bottom side.
	       for(j=(M-1)/2-(Lm-1)/2;j<=(M-1)/2+(Lm-1)/2;j++)
		  *(Tprint+400*(m++)+n)=*(source+i*M+j);
}

//IniMat_D: make every element of the matrix to be a constant.
void IniMat_D(int n,int m, double constant,double *matrix)
{int i=0,j=0;
   for(i=0;i<n;i++)
	   for(j=0;j<m;j++)
		   *(matrix+i*m+j)=constant;
	
}

void IniMat_I(int n,int m, int constant,int *matrix)
{int i=0,j=0;
   for(i=0;i<n;i++)
	   for(j=0;j<m;j++)
		   *(matrix+i*m+j)=constant;
	
}

int PriMatD(double *source,int height,int width,int HeiMat,int WidMat,FILE *fpt){
	int i,j;
	
		for(i=0;i<height;i++){
		for(j=0;j<width;j++){
			  if(j==width-1)
	          fprintf(fpt,"%lf\n",*(source+i*WidMat+j));
	          else
	          fprintf(fpt,"%lf ",*(source+i*WidMat+j));}}
	 return 0;
	
	}

int PriMatI(int *source,int height,int width,int HeiMat,int WidMat,FILE *fpt){
	int i,j;
	
		for(i=0;i<height;i++){
		for(j=0;j<width;j++){
			  if(j==width-1)
	          fprintf(fpt,"%d\n",*(source+i*WidMat+j));
	          else
	          fprintf(fpt,"%d ",*(source+i*WidMat+j));}}
	          
	 return 0;
	
	}
//Define the function, name to build a new file.
FILE * name(int N,int M,int *coefx, int *coefy,char *filename)
{
static FILE * fpt;
 if(N==1&&M==1){     
      *coefy=0;
      *coefx=0;
	   fpt=fopen(filename,"w"); 
     }
    else if(N==1&&M!=1) {
      *coefx=1;
	  *coefy=0;
	  fpt=fopen(filename,"w");
      }
      
      else if(N!=1&&M!=1){
      *coefx=1;
      *coefy=1;
	   fpt=fopen(filename,"w");}
   
   	 else if(N!=1&&M==1){
     *coefx=0;
	  *coefy=1;
	  fpt=fopen(filename,"w");
	  }
	  else
  		printf("Wrong");
  		return fpt; 
}

void Island(double *aNoise,double percentage,int N,int M,double lisland,double F)
   {
	//int gk=30; //min^-1
	//double gl=0.2; //min^-1
	//int vk0=-380; //mV
	//int vl0=-156; //mv
	//int sth=40;   //uM
	 int vth=-150; //mV
	//int alpha0=2;  //min^-1
	//double beta=1.3;  //min^-1
	//int mglobal=1;
	//int Fglobal=5600;  //uM/mV
	//double sigma=0.2;  //mV
	//double deltak=0.001; //mV/uM;
	//double deltal=0.008; //mV/uM;
	//double gammas=0.1;  //min^-1
	//int gammae=10;   //min^-1
	//int gammat=4;  //min^-1
	//int alphas=1;  //uM/min/mV
	//int alphat=1;   //uM/min/mV
	//int Dglobal=82800;   //1380*60um^2/min
	int mark[201][201],irdom=0,jrdom=0,badposition=0,Nrdom=0,i=0,j=0,Ntotal=0;
	double Temvth,cos,sin,width,wij,Rij0,kwidth=0.2;
	time_t tij;                  
	srand((unsigned) time(&tij));
	width=lisland/3;
	for(i=0;i<N;i++)
      for(j=0;j<M;j++)
        mark[i][j]=0;		  
	
	
	while(Nrdom<percentage*N*M&&Ntotal<=1000000){   // avoid infinite loop
		Ntotal=Ntotal+1;
		 
		if(Ntotal%100000==0) 
		printf("%d   %d   %lf\n",Ntotal,Nrdom,lisland);   
		
		
        irdom=rand()%(N-(int)lisland/3*4)+(int)lisland/3*2;   //Make sure there is engouth space to form a island with center in this point.
        jrdom=rand()%(M-(int)lisland/3*4)+(int)lisland/3*2;   //If lisland is even, then there are some small probblems: the distribution of irdom and jrdom will not be symmetric.
		Rij0=sqrt(pow(irdom-(N-1)/2,2)+pow(jrdom-(M-1)/2,2));
		
		if(Rij0>=lisland/2){    //Dont't want the island cross the center;		
		cos=(irdom-(N-1)/2)/Rij0;
		sin=(jrdom-(M-1)/2)/Rij0;
		badposition=0;
		for(i=irdom-(int)lisland/3*2;i<=irdom+(int)lisland/3*2;i++){
			for(j=jrdom-(int)lisland/3*2;j<=jrdom+(int)lisland/3*2;j++){
			wij=kwidth*((i-irdom)*cos+(j-jrdom)*sin)+width;   // write wij as a function of radial distance between two points.
			if(mark[i][j]==1&&fabs((i-irdom)*cos+(j-jrdom)*sin)<=(lisland)/2.0+0.00001&&fabs((i-irdom)*sin-(j-jrdom)*cos)<=(wij)/2+0.00001){
			badposition=1;
            break;	
			}   
			}
		    if(badposition==1)
		    break;
			}
		
		if(badposition==0){
			Temvth=-exp((double)rand()/2147483647.0*2*log(F)-log(F)+log(-vth)); 
			for(i=irdom-(int)lisland/3*2;i<=irdom+(int)lisland/3*2;i++){
			for(j=jrdom-(int)lisland/3*2;j<=jrdom+(int)lisland/3*2;j++){
			wij=kwidth*((i-irdom)*cos+(j-jrdom)*sin)+width;   // write wij as a function of R
			if(fabs((i-irdom)*cos+(j-jrdom)*sin)<=(lisland)/2+0.00001&&fabs((i-irdom)*sin-(j-jrdom)*cos)<=(wij)/2+0.00001){
			mark[i][j]=1;
            Nrdom++;
            *(aNoise+i*M+j)=Temvth;	
			}	
			}}
		}  //end of giving value to island.
	   }
	}  //end of while loop.
	
	
}

int CluGen(double *aNoise,double phi,double psi, int N,int M,double aLow,double DHigh)
{int cluster[N*M],neighbor[N*M],NClu,NNei,NNeiE,NHigh=0,i;
	for(i=0;i<N*M;i++)  cluster[i]=-1;       //Initiate it before the first use.
	for(i=0;i<N*M;i++)  neighbor[i]=-1;
	
  while(NHigh<phi*N*M){         //compare two doubles: Be careful for possible calculation errors of computer.
	 i=0;while(cluster[i]!=-1&&i<N*M)	  //Initiate the cluster. -1 is the end of effective part. -2 means that position is concealed.
		  cluster[i++]=-1;
     i=0; while(neighbor[i]!=-1&&i<N*M)
		neighbor[i++]=-1;
		 NClu=0;
		 NNei=0;       //The length of array neighbor;
		 NNeiE=0;       //The number of neighbor cells.
	 
	  do{
		cluster[NClu]=(rand()%N)*M+rand()%M;  
	  } while(*(aNoise+cluster[NClu])>aLow);   
	  
		*(aNoise+cluster[NClu++])=DHigh;
	     NHigh+=1;		
		
		while((double)rand()/(RAND_MAX+0.1)<psi&&NHigh<phi*N*M){     
		  for(i=NClu-1;i<NClu;i++){
			  if(*(aNoise+cluster[i]-1)<=aLow&&cluster[i]%M!=0&&InOut(cluster[i]-1,&neighbor[0]))
			  {neighbor[NNei++]=cluster[i]-1;
			  NNeiE++;}
			  if(*(aNoise+cluster[i]+1)<=aLow&&cluster[i]%M!=M-1&&InOut(cluster[i]+1,&neighbor[0]))
			  {neighbor[NNei++]=cluster[i]+1;
			  NNeiE++;}
			  if(*(aNoise+cluster[i]-M)<=aLow&&cluster[i]/M!=0&&InOut(cluster[i]-M,&neighbor[0]))
			  {neighbor[NNei++]=cluster[i]-M;
		       NNeiE++;}
			  if(*(aNoise+cluster[i]+M)<=aLow&&cluster[i]/M!=N-1&&InOut(cluster[i]+M,&neighbor[0]))
			  {neighbor[NNei++]=cluster[i]+M;
		        NNeiE++;}			  
			  }   
		  
		   if(NHigh>=N*M-3)  printf("%d %d %d %d\n",NHigh,NClu,NNei,NNeiE);//$$BreakPoint     
		  
		       if(NNeiE==0)
		        break;
			  
 			     i=0;
		         while(1){
				   if(neighbor[i]!=-2&&(double)rand()/(RAND_MAX+0.1)<1.0/(double)NNeiE){
				    *(aNoise+neighbor[i])=DHigh;
					cluster[NClu++]=neighbor[i];
				    NHigh++;
					neighbor[i]=-2;      //canceal this selected point from neighbor.
					NNeiE=NNeiE-1;
					break;}
					 if(neighbor[++i]==-1)
				     i=0;}
				     
		 
			    
		}  //end of while((double)rand()/2147483647.0<=psi&&space==1);	end of generating one cluster.  
  }  //end of while(NHigh<phi*N*M);  completing all clusters.
  return 1;
}

int CluGen_G5_1(double *nvth,int height,int width,double phi,double rho,double DLow,double DHigh,
              double alpha, double beta, double a1, double a2, double a3){

 int i,j,m, HeiMat=height+50, WidMat=width;
double pon, poff, ponon,ponoff, poffon, poffoff, pononon, pononoff, ponoffon, ponoffoff, possibility[width],Prenvth[HeiMat*WidMat];

    pon=phi;
    poff=1-pon;
	ponon=phi/(1-rho+rho*phi);   //p(d:on|m:on): rho=0, no spatial correlation. rho=1, perfect lineage inheritance
	poffon=1-ponon;	
	ponoff=ponon*(1-rho);        //p(d:on|m:off)
    poffoff=1-ponoff;
	
    pononoff=ponon/(alpha*ponon+ponoff*poff/pon);
    pononon=alpha*pononoff;      //alpha, beta: multiple generation inheritance: the influence of gradmother. if alpha=1, no direct inheritance from grandmother;
    ponoffoff=ponoff/(beta*poffon*pon/poff+poffoff);
	ponoffon=beta*ponoffoff;   

  for(i=0;i<HeiMat;i++){
	if(i==0){
		for(j=0;j<width;j++)
		Prenvth[j]=AssignPhi_MT(phi,DLow,DHigh);
	 }
				
	else if(i==1){			
		for(j=0;j<width;j++){
			m=i*WidMat+j;
			if(Prenvth[m-WidMat]==DHigh)   
			 Prenvth[m]=AssignPhi_MT(ponon,DLow,DHigh);				
			else
			 Prenvth[m]=AssignPhi_MT(ponoff,DLow,DHigh);    
		      }
	}
	
	else{
		mothergrandmother('O', Prenvth, HeiMat, WidMat, i,DLow, DHigh, possibility, ponoffoff, ponoffon, pononoff, pononon);
		
		for(j=1;j<width;j=j+2){    //a1,a2,a3: describe the mobility of biofilm: the daughter has probability a2 to be a daughter of the top mother, a1 of the left mother and a2 of the right mother.
			m=i*WidMat+j;
			if(j==0)
			Prenvth[m]=AssignPhi_MT((a1+a2)*possibility[j]+a3*possibility[j+1],DLow, DHigh);
			else if(j==width-1)
			Prenvth[m]=AssignPhi_MT((a2+a3)*possibility[j]+a1*possibility[j-1],DLow, DHigh);
			else
			Prenvth[m]=AssignPhi_MT(a1*possibility[j-1]+a2*possibility[j]+a3*possibility[j+1],DLow, DHigh);
			}
			
		mothergrandmother('E', Prenvth, HeiMat, WidMat, i,DLow, DHigh, possibility, ponoffoff, ponoffon, pononoff, pononon);
		for(j=0;j<width;j=j+2){
			m=i*WidMat+j;
			if(j==0)
			Prenvth[m]=AssignPhi_MT((a1+a2)*possibility[j]+a3*possibility[j+1],DLow, DHigh);
			else if(j==width-1)
			Prenvth[m]=AssignPhi_MT((a2+a3)*possibility[j]+a1*possibility[j-1],DLow, DHigh);
			else
			Prenvth[m]=AssignPhi_MT(a1*possibility[j-1]+a2*possibility[j]+a3*possibility[j+1],DLow, DHigh);
			}			
		}
	}
	
	
	for(i=0;i<height;i++){
		for(j=0;j<width;j++){
			m=i*width+j;
	        nvth[m]=Prenvth[(i+HeiMat-height)*WidMat+j]; } }
		
	return 0;
 }
 
int CluRadLin(double *nvth, int height, int width, char mode, double rho, double phi, double DLow, 
             double median, double lambda){
				 
int cluster[height*width], state[height*width], NClu=0, NClu0, NEffClu=0, NEffClu0, NMother=0, daughter,mother, StoLoc, ABLR,m, i,j;
	IniMat_I(height, width, -999, &cluster[0]);
	IniMat_I(height, width, -999, &state[0]);
	IniMat_D(height, width, -999, nvth);				 

	
 if(mode=='s'){
	 for (m=0; m<height*width; m++)
	    nvth[m]=median;          
	          }
 else{
	 
	//Initiation: seeding with 6*6 grid of cells, same as the trigger, drawn independently from p(m).
	for (i=height/2-3; i<=height/2+2; i++){
			for(j=width/2-3; j<=width/2+2; j++){
	mother=i*width+j;
	if (mode=='b')
	  nvth[mother]=AssignPhi_MT(phi, DLow, median);
	else if(mode=='N')
	  nvth[mother]=exp( gaussrand_K()*lambda+ median );
	else if(mode=='U')
	  nvth[mother]=exp( DisUni( median, lambda ) );
	 cluster[NClu]=mother;
	 state[NClu++]=0;
	 NEffClu++;
           }
    }
    NClu0=NClu;
	NEffClu0=NEffClu;

	do{  //produce the 2nd, 3rd... generations until the height by width box is full.
		
		while(NMother<NEffClu0)  //produce the next generation from the mothers.  
		{ 
			do{StoLoc=rand()%NClu0;  //Randomly select one mother which has unoccupied neighbors.
			  } while(state[StoLoc]!=0);
			mother=cluster[StoLoc];
			NMother++;
			//Decide whether this mother have unoccupied neighbors to give birth to one daughter
			if( (mother/width!=0&& nvth[mother-width]==-999) || (mother/width!=(height-1) && nvth[mother+width]==-999) || 
			    (mother%width!=0&& nvth[mother-1]==-999) || (mother%width!=(width-1) && nvth[mother+1]==-999) ){
					 do{ABLR=rand()%4;
					    daughter=(mother-width)*(ABLR==0)+(mother+width)*(ABLR==1)+(mother-1)*(ABLR==2)+(mother+1)*(ABLR==3);
					   }  while( ((mother/width==0)&&(ABLR==0))|| ((mother/width==(height-1))&&(ABLR==1)) || ((mother%width==0)&&(ABLR==2)) 
					             || ((mother%width==(width-1))&&(ABLR==3)) || (nvth[daughter]!=-999) );
					 //Assign a number to the daughter: has probability of rho to be the same with the mother; 1-rho to be uncorrelated.
				if(rand()/(RAND_MAX+0.1)<rho)
				   nvth[daughter]=nvth[mother];
				else{ 	
					 if (mode=='b')
						  nvth[daughter]=AssignPhi_MT(phi, DLow, median);
						else if(mode=='N')
						  nvth[daughter]=exp( gaussrand_K()*lambda+ median );
						else if(mode=='U')
						  nvth[daughter]=exp( DisUni( median,lambda ) ); 
				    }
				     //Put daughter into cluster, increase the NClu and NEffClu	    
			    cluster[NClu]=daughter;
			    state[NClu++]=0;
			    NEffClu++;
			    state[StoLoc]=1;  //mark that this mother has generated daughter 
					}
			else  {        //Mark the mother as unable to give birth to daughters.
				state[StoLoc]=-1;
		        NEffClu--;       }   
		}
		//printf("NClu %d, NClu0 %d, NEffClu %d, NEffClu0 %d, NMother %d\n",NClu, NClu0, NEffClu, NEffClu0, NMother);	
		//for (m=0;m<NClu;m++)		  
		//  printf("Storage location %d, mother %d, value of mother %lf, state of mother %d\n",m,cluster[m], nvth[cluster[m]], state[m]);
	      //if(NClu==height*width){		
	    	//exit(0);}
		//Recover the variables to be ready for the production of the next generation.
		NMother=0;
		for (m=0;m<NClu0;m++)
		  state[m]=0*(state[m]==1)-1*(state[m]==-1);
		NEffClu0=NEffClu;
		NClu0=NClu;
			
			
		
	  } while (NEffClu>0);    
 } //end of mode 'b','U','N'
	
	return 0;	 
}
 
 
 int mothergrandmother(char mode, double *nvth, int HeiMat, int WidMat, int height, double DLow, double DHigh,
           double *possibility, double ponoffoff, double ponoffon, double pononoff, double pononon){
int j,m;
if (mode=='O'){
	for(j=0;j<WidMat;j++){
	 m=height*WidMat+j;
	 if(nvth[m-WidMat]==DLow&& nvth[m-2*WidMat]==DLow)
	   possibility[j]=ponoffoff;
	 else if(nvth[m-WidMat]==DLow&& nvth[m-2*WidMat]==DHigh)
	   possibility[j]=ponoffon;
	 else if(nvth[m-WidMat]==DHigh&& nvth[m-2*WidMat]==DLow)
	   possibility[j]=pononoff;
	 else 
	   possibility[j]=pononon;
	 }
	 }
	 
else{
  for(j=0;j<WidMat;j=j+2){
	 m=height*WidMat+j;
	 if(nvth[m-WidMat]==DLow&& nvth[m-2*WidMat]==DLow)
	   possibility[j]=ponoffoff;
	 else if(nvth[m-WidMat]==DLow&& nvth[m-2*WidMat]==DHigh)
	   possibility[j]=ponoffon;
	 else if(nvth[m-WidMat]==DHigh&& nvth[m-2*WidMat]==DLow)
	   possibility[j]=pononoff;
	 else 
	   possibility[j]=pononon;
	 }
	 
  for(j=1;j<WidMat;j=j+2){
	 m=height*WidMat+j;
	 if(nvth[m]==DLow&& nvth[m-WidMat]==DLow)
	   possibility[j]=ponoffoff;
	 else if(nvth[m]==DLow&& nvth[m-WidMat]==DHigh)
	   possibility[j]=ponoffon;
	 else if(nvth[m]==DHigh&& nvth[m-WidMat]==DLow)
	   possibility[j]=pononoff;
	 else 
	   possibility[j]=pononon;
	 }
	 }
	return 0;
}

double AssignPhi(double phi,double low,double high){
	         if((double)rand()/(RAND_MAX+0.1)<phi)
			    return high;
			    else
			    return low;
}

double AssignPhi_MT(double phi,double low,double high){
	   if( RanGen.randExc()< phi ){
	          return high;
	   }else
	      return low;
}


int LargestCluster(double *nvth, int HeiMat, int WidMat, int height,int width, int row, double DLow, double DHigh, int square){
	int j,state[HeiMat][WidMat],cluster[HeiMat][WidMat], clustersize[WidMat];
	int MaxColumn=0;
	
	IniMat_I(HeiMat,WidMat, 0, &state[0][0]);
	IniMat_I(HeiMat,WidMat, 0, &cluster[0][0]);
	 for(j=0;j<=width-square;j++){
	  clustersize[j]=SearchTree_G4(nvth, &state[0][0],row,j,height, width, HeiMat, WidMat,
                   DHigh, &cluster[0][0]);
      if (j==0)       
       MaxColumn=0;
       else if(clustersize[j]>clustersize[MaxColumn])
       MaxColumn=j;       
             }
     
    return MaxColumn;              
      
	   }

//We need to initiate the matrix state before using this function. 0: unchecked. 1: checked to be non-targe. 2. checked to be the target. 	   
int SearchClu_G4_Squ(double *nvth,int *state,int i,int j,int height,int width,int HeiMat,int WidMat,
                  double DHigh,int *cluster)  //il is the dimension of effective area of nvth and state;
{  int column=i*WidMat+j,i1,j1,m,p;  //i1: the No. of row; j1, the No. of column; m, the No. in the srorage, m= i1*WidMat+j1;
	if(state[column]==0&& nvth[column]==DHigh){
		  
	     int NClu0=0,NClu=0,NNei=0;  //The length of array neighbor;    //The number of effictive.neighbor cells.       	
		  cluster[NClu++]=column;
		  NNei+=1;
			state[column]=2;     //for(p=0;p<NClu;p++) printf("%d ",cluster[p]); printf("\n");
		     
		   while(1){
		   for(p=NClu0;p<NClu;p++){    //look for the neighbors of these new cells in cluster. 
			   m=cluster[p];
			   i1=m/WidMat;
			   j1=m%WidMat;
			   
		  if(nvth[m-1]==DHigh&& state[m-1]==0&& j1!=0)
		  {cluster[NNei++]=m-1;	
			 state[m-1]=2; }
			 
		  if(nvth[m+1]==DHigh&& state[m+1]==0&& j1!=width-1)
		  {cluster[NNei++]=m+1;	
			 state[m+1]=2; }
			 
		 if(nvth[m-WidMat]==DHigh&& state[m-WidMat]==0&& i1!=0)
		  {cluster[NNei++]=m-WidMat;	
			 state[m-WidMat]=2; }	
			 
		 if(nvth[m+WidMat]==DHigh&& state[m+WidMat]==0&& i1!=height-1)
		  {cluster[NNei++]=m+WidMat;	
			 state[m+WidMat]=2; }		
		                            } 
		                                          		       
		          if(NNei==NClu)			   		     
				  return NClu;
								
 					  
              NClu0=NClu;
              NClu=NNei;
		    }
	}
	
   else if(*(nvth+column)!=DHigh)
        *(state+column)=1;       
    
   return 0;	   
}

int SearchClu_G4_Squ_Frame(double *nvth,int *state,int i,int j, int RowSta, int RowEnd, int ColSta, int ColEnd, 
       int HeiMat,int WidMat, double DHigh,int *cluster)  //il is the dimension of effective area of nvth and state;
{  int column=i*WidMat+j,i1,j1,m,p;  //i1: the No. of row; j1, the No. of column; m, the No. in the srorage, m= i1*WidMat+j1;
	if(state[column]==0&& nvth[column]==DHigh){
		  
	     int NClu0=0,NClu=0,NNei=0;  //The length of array neighbor;    //The number of effictive.neighbor cells.       	
		  cluster[NClu++]=column;
		  NNei+=1;
			state[column]=2;     //for(p=0;p<NClu;p++) printf("%d ",cluster[p]); printf("\n");
		     
		   while(1){
		   for(p=NClu0;p<NClu;p++){    //look for the neighbors of these new cells in cluster. 
			   m=cluster[p];
			   i1=m/WidMat;
			   j1=m%WidMat;
			   
		  if(j1!=ColSta&& nvth[m-1]==DHigh&& state[m-1]==0)  //Not in the first column
		  {cluster[NNei++]=m-1;	
			 state[m-1]=2; }
			 
		  if(j1!=ColEnd-1&& nvth[m+1]==DHigh&& state[m+1]==0) //Not in the last column
		  {cluster[NNei++]=m+1;	
			 state[m+1]=2; }
			 
		 if(i1!=RowSta&& nvth[m-WidMat]==DHigh&& state[m-WidMat]==0)
		  {cluster[NNei++]=m-WidMat;	
			 state[m-WidMat]=2; }	
			 
		 if(i1!=RowEnd-1&& nvth[m+WidMat]==DHigh&& state[m+WidMat]==0)
		  {cluster[NNei++]=m+WidMat;	
			 state[m+WidMat]=2; }		
		                            } 
		                                          		       
		          if(NNei==NClu)			   		     
				  return NClu;
								
 					  
              NClu0=NClu;
              NClu=NNei;
		    }
	}
	
   else if(*(nvth+column)!=DHigh)
        *(state+column)=1;       
    
   return 0;	   
}

int SearchClu_G4_Tri_Frame(double *nvth,int *state,int i,int j,int RowSta, int RowEnd, int ColSta, int ColEnd, 
       int HeiMat,int WidMat, double DHigh,int *cluster){
	int column=i*WidMat+j,i1,j1,m,p;  //i1: the No. of row; j1, the No. of column; m, the No. in the srorage, m= i1*WidMat+j1;
	if(state[column]==0&& nvth[column]==DHigh){
		  
	     int NClu0=0,NClu=0,NNei=0;  //The length of array neighbor;    //The number of effictive.neighbor cells.       	
		  cluster[NClu++]=column;
		  NNei+=1;
			state[column]=2;     //for(p=0;p<NClu;p++) printf("%d ",cluster[p]); printf("\n");
		     
		   while(1){
		   for(p=NClu0;p<NClu;p++){    //look for the neighbors of these new cells in cluster. 
			   m=cluster[p];
			   i1=m/WidMat;
			   j1=m%WidMat;
			   
		  if(j1!=ColSta&& nvth[m-1]==DHigh&& state[m-1]==0)
		  {cluster[NNei++]=m-1;	
			 state[m-1]=2; }
			 
		  if(j1!=ColEnd-1&& nvth[m+1]==DHigh&& state[m+1]==0)
		  {cluster[NNei++]=m+1;	
			 state[m+1]=2; }
			 
		 if(i1!=RowSta&& nvth[m-WidMat]==DHigh&& state[m-WidMat]==0)
		  {cluster[NNei++]=m-WidMat;	
			 state[m-WidMat]=2; }	
			 
		 if(i1!=RowEnd-1&& nvth[m+WidMat]==DHigh&& state[m+WidMat]==0)
		  {cluster[NNei++]=m+WidMat;	
			 state[m+WidMat]=2; }	
				   
		 if(j1%2==0){
			 if(i1!=RowEnd-1&& j1!=ColSta&& nvth[m+WidMat-1]==DHigh&& state[m+WidMat-1]==0)
		     {cluster[NNei++]=m+WidMat-1;	
			 state[m+WidMat-1]=2; }
			 
			 if(i1!=RowEnd-1&& j1!=ColEnd-1&& nvth[m+WidMat+1]==DHigh&& state[m+WidMat+1]==0)
		     {cluster[NNei++]=m+WidMat+1;	
			 state[m+WidMat+1]=2; }	   }
			 
			else{
			if(i1!=RowSta&& j1!=ColSta&& nvth[m-WidMat-1]==DHigh&& state[m-WidMat-1]==0)
		     {cluster[NNei++]=m-WidMat-1;	
			 state[m-WidMat-1]=2; }
			 
			 if(i1!=RowSta&& j1!=ColEnd-1&& nvth[m-WidMat+1]==DHigh&& state[m-WidMat+1]==0)
		     {cluster[NNei++]=m-WidMat+1;	
			 state[m-WidMat+1]=2; }			}	
		                            } 
		                                          		       
		          if(NNei==NClu)			   		     
				      return NClu;
								
 					  
              NClu0=NClu;
              NClu=NNei;
		    }
	}
	
   else if(*(nvth+column)!=DHigh)
        *(state+column)=1;       
    
   return 0;
}

int SearchCon_G4_Tri_Frame(double *nvth,int *state,int i,int j,int RowSta, int RowEnd, int ColSta, int ColEnd, 
       int HeiMat,int WidMat, double DHigh,int *cluster){
	int column=i*WidMat+j,i1,j1,m,p;  //i1: the No. of row; j1, the No. of column; m, the No. in the srorage, m= i1*WidMat+j1;
	if(state[column]==0&& nvth[column]==DHigh){
		  
	     int NClu0=0,NClu=0,NNei=0;  //The length of array neighbor;    //The number of effictive.neighbor cells.       	
		  cluster[NClu++]=column;
		  NNei+=1;
			state[column]=2;     //for(p=0;p<NClu;p++) printf("%d ",cluster[p]); printf("\n");
		     
		  while(1){
		   for(p=NClu0;p<NClu;p++){    //look for the neighbors of these new cells in cluster. 
			   m=cluster[p];
			   i1=m/WidMat;
			   j1=m%WidMat;
			   
		  if(j1!=ColSta&& nvth[m-1]==DHigh&& state[m-1]==0)
		  {cluster[NNei++]=m-1;	
			 state[m-1]=2; }
			 
		  if(j1!=ColEnd-1&& nvth[m+1]==DHigh&& state[m+1]==0)
		  {cluster[NNei++]=m+1;	
			 state[m+1]=2; }
			 
		 if(i1!=RowSta&& nvth[m-WidMat]==DHigh&& state[m-WidMat]==0)
		  {cluster[NNei++]=m-WidMat;	
			 state[m-WidMat]=2; }	
			 
		 if(i1!=RowEnd-1&& nvth[m+WidMat]==DHigh&& state[m+WidMat]==0)
		  {cluster[NNei++]=m+WidMat;	
			 state[m+WidMat]=2; }	
				   
		 if(j1%2==0){
			 if(i1!=RowEnd-1&& j1!=ColSta&& nvth[m+WidMat-1]==DHigh&& state[m+WidMat-1]==0)
		     {cluster[NNei++]=m+WidMat-1;	
			 state[m+WidMat-1]=2; }
			 
			 if(i1!=RowEnd-1&& j1!=ColEnd-1&& nvth[m+WidMat+1]==DHigh&& state[m+WidMat+1]==0)
		     {cluster[NNei++]=m+WidMat+1;	
			 state[m+WidMat+1]=2; }	   }
			 
			else{
			if(i1!=RowSta&& j1!=ColSta&& nvth[m-WidMat-1]==DHigh&& state[m-WidMat-1]==0)
		     {cluster[NNei++]=m-WidMat-1;	
			 state[m-WidMat-1]=2; }
			 
			 if(i1!=RowSta&& j1!=ColEnd-1&& nvth[m-WidMat+1]==DHigh&& state[m-WidMat+1]==0)
		     {cluster[NNei++]=m-WidMat+1;	
			 state[m-WidMat+1]=2; }			}	
		                            }
		     //all cells are collected in cluster.                      
		        if(NNei==NClu){
				   for(p=0;p<NClu;p++){
				     if(cluster[p]/WidMat==HeiMat-1)  //If one element is in the last row;
				     return 1;   //Connected
				   }				     
				   return 0; //unconnected
				}						
 					  
              NClu0=NClu;
              NClu=NNei;
		    }
	}
	
   else if(*(nvth+column)!=DHigh)
        *(state+column)=1;       
    
   return 0; //unconnected
}


int SearchTree_G4(double *nvth,int *state,int i,int j,int height,int width,int HeiMat,int WidMat,
                  double DHigh,int *cluster)  //il is the dimension of effective area of nvth and state;
{  int column=i*WidMat+j,i1,j1,m,p;  //i1: the No. of row; j1, the No. of column; m, the No. in the srorage, m= i1*WidMat+j1;
	if(state[column]==0&& nvth[column]==DHigh){
		  
	     int NClu0=0,NClu=0,NNei=0;  //The length of array neighbor;    //The number of effictive.neighbor cells.       	
		  cluster[NClu++]=column;
		  NNei+=1;
			state[column]=2;     //for(p=0;p<NClu;p++) printf("%d ",cluster[p]); printf("\n");
		     
		   while(1){
		   for(p=NClu0;p<NClu;p++){    //look for the neighbors of these new cells in cluster. 
			   m=cluster[p];
			   i1=m/WidMat;
			   j1=m%WidMat;
			   
		  if(nvth[m-1]==DHigh&& state[m-1]==0&& j1!=0)
		  {cluster[NNei++]=m-1;	
			 state[m-1]=2; }
			 
		  if(nvth[m+1]==DHigh&& state[m+1]==0&& j1!=width-1)
		  {cluster[NNei++]=m+1;	
			 state[m+1]=2; }
			 
		 if(nvth[m-WidMat]==DHigh&& state[m-WidMat]==0&& i1!=0)
		  {cluster[NNei++]=m-WidMat;	
			 state[m-WidMat]=2; }	
			 
		 if(nvth[m+WidMat]==DHigh&& state[m+WidMat]==0&& i1!=height-1)
		  {cluster[NNei++]=m+WidMat;	
			 state[m+WidMat]=2; }	
				   
		 if(j1%2==0){
			 if(nvth[m+WidMat-1]==DHigh&& state[m+WidMat-1]==0&& i1!=height-1&& j1!=0)
		     {cluster[NNei++]=m+WidMat-1;	
			 state[m+WidMat-1]=2; }
			 
			 if(nvth[m+WidMat+1]==DHigh&& state[m+WidMat+1]==0&& i1!=height-1&& j1!=width-1)
		     {cluster[NNei++]=m+WidMat+1;	
			 state[m+WidMat+1]=2; }	   }
			 
			else{
			if(nvth[m-WidMat-1]==DHigh&& state[m-WidMat-1]==0&& i1!=0&& j1!=0)
		     {cluster[NNei++]=m-WidMat-1;	
			 state[m-WidMat-1]=2; }
			 
			 if(nvth[m-WidMat+1]==DHigh&& state[m-WidMat+1]==0&& i1!=0&& j1!=width-1)
		     {cluster[NNei++]=m-WidMat+1;	
			 state[m-WidMat+1]=2; }			}	
		                            } 
		                                          		       
		          if(NNei==NClu)			   		     
				  return NClu;
								
 					  
              NClu0=NClu;
              NClu=NNei;
		    }
	}
	
   else if(*(nvth+column)!=DHigh)
        *(state+column)=1;       
    
   return 0;	   
}


int InOut(int target,int *source){   //If target is not in the array source, return 1; otherwise return 0.
	int i=0,sign=1;
	while(source[i]!=-1){
		if(target==source[i++]){
			sign=0;
			break;
		}			
	}
	return sign;		
}

int WheIn(int target,int *source, int length){   //If target is not in the array source, return 1; otherwise return 0.
	int i=0,sign=0;
	for(i=0; i<length; i++){
		if(target==source[i]){
			sign=1;
			break;
		}			
	}
	return sign;		
}

double DisUni(double mean, double HalWid){
	double RanNum;
	RanNum=(double)rand()/(RAND_MAX);
	RanNum=RanNum*2*HalWid-HalWid+mean;
	return RanNum;	
	}
	
int equal(double a, double b, double error){
	if (fabs(a-b)<=error)
	 return 1;
	else 
	 return 0;}  
	 
int WheMultiple(double numerator, double denominator, double error){
	if (fabs(numerator-denominator*(int)((numerator+denominator/2)/denominator))<=error)
	  return 1;
	else
	  return 0;}
	  
int MeaLam(double *wave, int height, int width, int TimStep, double threshold){
  int i,j, it;
  for (i=height-1;i>=0; i--){
	    for(j=0;j<width;j++){
			for(it=0;it<TimStep;it++){
				if( wave[(it*height+i)*width+j]>=threshold )
				  return i+1;}}
	  }
	  
	  return 0;}
	  
int wavefront(double *wave, int height, int width, int ColSta, int ColEnd, double threshold){
  int i,j;
  for (i=height-1;i>=1; i--){
	    for(j=ColSta;j<=ColEnd;j++){
		if( wave[i*width+j]>=threshold )
				  return i;}}
  
	  return 0;}	
	  
int waveback(double *wave, int height, int width, int ColSta, int ColEnd, double threshold){
  int i,j;
  for (i=1;i<height; i++){
	    for(j=ColSta;j<=ColEnd;j++){
		if( wave[i*width+j]>=threshold )
				  return i;}}
  
	  return 99;}

double wavelength(double *wave, int height, int width, int LenSec, double threshold){
int i, WF, WB;
double WL=0;
	for(i=0;i<width;i=i+LenSec){
		WF=wavefront(wave, height, width, i, i+LenSec-1, threshold);
		WB=waveback(wave, height, width, i, i+LenSec-1, threshold);
		if(WF!=0 && WB!=99)
		WL=WL+ WF-WB+1;
		else
		WL=WL+0;               }
		
	WL=WL/(width/LenSec);
	return WL;
		}
	  
int NFir(double *wave, int height, int width, double threshold){
	int i,j, N_firing=0;
	
	for (i=0;i<height; i++)
	    for(j=0;j<width;j++)
		 N_firing+= ( wave[i*width+j]>=threshold );
		 
		 return N_firing;
			}	
			
double TMean(double *Tij, int height, int width, double threshold){
	int i,j, N_firing=0;
	double T_sum=0;
	
	for (i=0;i<height; i++)
	    for(j=0;j<width;j++){
		 N_firing+= ( Tij[i*width+j]>=threshold );
		 T_sum+= Tij[i*width+j]*( Tij[i*width+j]>=threshold );
		   }
		 
		 if(N_firing>=1)
		 return T_sum/N_firing;
		 else 
		 return 0;}  
		 
double gaussrand_K()
{
	static double V1, V2, S;
	static int phase = 0;
	double X;

	if(phase == 0) {
		do {
			double U1 = (double)rand() / RAND_MAX;
			double U2 = (double)rand() / RAND_MAX;

			V1 = 2 * U1 - 1;
			V2 = 2 * U2 - 1;
			S = V1 * V1 + V2 * V2;
			} while(S >= 1 || S == 0);

		X = V1 * sqrt(-2 * log(S) / S);
	} else
		X = V2 * sqrt(-2 * log(S) / S);

	phase = 1 - phase;

	return X;
}


double GyrationRadius(int *array, int Height, int Width, int length){
double GR=0, i0=0, j0=0;  //Accumulative variable, initiating to be 0
int k;
int location_ij[2][length];  //Use a matrix to store i,j to avoid repeating calculating.
for (k=0; k<length; k++){
	location_ij[0][k]=array[k]/Width;  //location of cell k, r_k
	location_ij[1][k]=array[k]%Width;
	i0+=location_ij[0][k];
	j0+=location_ij[1][k];
}
i0/=length;
j0/=length;

for (k=0; k<length; k++)
	GR=GR+ (location_ij[0][k]- i0)*(location_ij[0][k]- i0)+(location_ij[1][k]- j0)*(location_ij[1][k]- j0);
GR= sqrt(GR/length); //Rg=1/ 2N^2 sum_i,j (r_i-r_j)^2 

return GR;
}	 

double GyrationRadius_Tri(int *array, int Height, int Width, double dy, double dx, int length){
double GR=0,  y0=0,x0=0;  //Accumulative variable, initiating to be 0
int k,i=0, j=0;
double location_ij[2][length];  //Use a matrix to store y,x to avoid repeating calculating.
for (k=0; k<length; k++){
	i=array[k]/Width;  //#row
	j=array[k]%Width;  //#column
	
	location_ij[0][k]= (i*dy+dy/2)*(j%2 ==0)+(i*dy)*(j%2 ==1);  //y(i)
	location_ij[1][k]= j*dx;                                    //x(j)
	
	y0 +=location_ij[0][k];
	x0 +=location_ij[1][k];
}
y0/=length;
x0/=length;

for (k=0; k<length; k++)
	GR=GR+ (location_ij[0][k]- y0)*(location_ij[0][k]- y0)+(location_ij[1][k]- x0)*(location_ij[1][k]- x0);
GR= sqrt(GR/length); //Rg=1/ 2N^2 sum_i,j (r_i-r_j)^2 

return GR;
} 
			  
    
int FindPeaks(double deltat, double TimLen, double dt_write, double *trace, int LocSta, int LocEnd, const int TimStep, double v_th, double *PeakTim, int LenPeakTime){
	int i=0, NumPeaks=0, state;	
	for (i=LocSta; i<LocEnd; i++){
		if (NumPeaks>= LenPeakTime-1)
		  break;
		  
	  state=(i+1)*(i<=1)-(TimStep-i)*(i>=TimStep-2);  //First, second and the last secod last moments.
	  switch(state) {
      case 1 :
        if(trace[i]>= trace[i+1] && trace[i]>= trace[i+2]&& trace[i]>= v_th)
		  PeakTim[NumPeaks++]=2*deltat; 
         break;
      case 2 :
         if(trace[i]>= trace[i-1] && trace[i]>= trace[i+1] && trace[i]>= trace[i+2]&& trace[i]>= v_th)			
		  PeakTim[NumPeaks++]=i*dt_write; 
		  break;
      case ( -1 ) :
         if(trace[i]>= trace[i-1] && trace[i]>= trace[i-2]&& trace[i]>= v_th)
		  PeakTim[NumPeaks++]=i*dt_write; 
         break;
      case ( -2 ) :
         if(trace[i]>= trace[i-1] && trace[i]>= trace[i-2] && trace[i]>= trace[i+1]&& trace[i]>= v_th)			
		  PeakTim[NumPeaks++]=i*dt_write;
         break;
      
      default :
         if(trace[i]>= trace[i-1] && trace[i]>= trace[i-2] && trace[i]>= trace[i+1]&& trace[i]>= trace[i+2]&& trace[i]>= v_th)			
		  PeakTim[NumPeaks++]=i*dt_write;
                  }		
		}
		
	PeakTim[LenPeakTime-1]=NumPeaks;
	
	return 0;	
	}
	
double PerCelMaxPea(double *PerCMP, int NumCMP){
	int i=0, j;
	double mean=0, deviation[NumCMP];
	
	for (i=0; i<NumCMP; i++)
		mean+= PerCMP[i];
	mean/= NumCMP;
	
	for (i=0; i<NumCMP; i++)
	    deviation[i]= fabs( PerCMP[i]- mean); 
	    
	for(i=1; i< NumCMP; i++){
		j=i;
		while(j>=1 && deviation[j]< deviation[j-1]){
			ExchangeD(&deviation[j], &deviation[j-1]);
			ExchangeD(&PerCMP[j], &PerCMP[j-1]);
			j--;
			}
	}
	
	mean=0;
	for(i=0; i< NumCMP/2; i++)
	  mean+= PerCMP[i];
	mean/= (NumCMP/2);
	return mean;  
	}
	
int WheDied(double *wave, int height, int width, int TimStep, double v_th){
	int i,j, it;
   for(i=0; i<=height-1; i=i+height-1){
	for (j=0; j<width; j++){
		for (it=0; it<TimStep; it++){
			if (wave[ (i+it*height) *width +j]>=v_th) 
			  return 0;} 
  } 	}
			  
    for(i=0; i<height; i++){
	for (j=0; j<width; j=j+width-1){
		for (it=0; it<TimStep; it++){
			if (wave[ (i+it*height) *width +j]>=v_th) 
			  return 0;} 
  } 	}
  
  return 1;	
	}
	
int WheDied_Tim(double *wave, int height, int width, int TimStep, double v_th){   //died, -1; non-died, return the time sequence. 
	int i,j, it;
  for (it=0; it<TimStep; it++){	
	  
   for(i=0; i<height; i=i+height-1){
	for (j=0; j<width; j++){		
			if (wave[ (i+it*height) *width +j]>=v_th) 
			  return it;} 
   }
  
  for(i=0; i<height; i++){
	for (j=0; j<width; j=j+width-1){
		if (wave[ (i+it*height) *width +j]>=v_th) 
			  return it;
  } 	
  } 	
  
  
  }
			  
    
  
  return -1;	
	}
	
int WheNeverDied(double *wave, int height, int width, int TimStep, double v_th){
 int i,j, it=TimStep-1;
	for(i=0; i<height; i++){
	for (j=0; j<width; j++){
		if ( wave[ (i+it*height) *width +j]>=v_th) 
			  return 1;} 
  } 
  return 0;
  	}
  	
int ConvexHull(double *SetX, int SizeSetX, int length, double *ConHull, int LenConHull){
	struct ConHulPoint SorPoi[SizeSetX];
	int i, j;
	for (i=0; i< SizeSetX; i++){
		SorPoi[i].xaxis= SetX[i];
		SorPoi[i].yaxis= SetX[i+length];
		SorPoi[i].cos= 0;
		}
	//Find out A0: the 1)Lowest 2)most left cell 	
	int A0=0;
	for (i=0; i<SizeSetX; i++){
		if( SorPoi[i].yaxis < SorPoi[A0].yaxis  ||  (SorPoi[i].yaxis == SorPoi[A0].yaxis && SorPoi[i].xaxis < SorPoi[A0].xaxis) )  
	      A0=i;
	 }
	//Put A0 at the beginning of the array. 
	ExchangeD(&SorPoi[0].xaxis, &SorPoi[A0].xaxis);
	ExchangeD(&SorPoi[0].yaxis, &SorPoi[A0].yaxis);
	A0=0;
	
	//Calculate the cosine of angle between A0Ai and x axis, vector (1,0)
	SorPoi[A0].cos=1;
	for (i=1; i<SizeSetX ; i++){
	   SorPoi[i].cos=(SorPoi[i].xaxis - SorPoi[A0].xaxis )
	           /sqrt( pow(SorPoi[i].xaxis -SorPoi[A0].xaxis, 2) + pow(SorPoi[i].yaxis - SorPoi[A0].yaxis, 2) ); }
	           
	//Sort the array SorPoi[i] by deceeding the SorPoi[i].cos
	for (i=1; i<SizeSetX; i++){
		j=i;
		//If angle(A0Aj,x) < angle(A0Aj-1, x) switch. If the angle is the same, but the distance(A0Aj)< distance(A0Aj-1), switch.  
		while (j>=1 &&(SorPoi[j].cos > SorPoi[j-1].cos || (SorPoi[j].cos == SorPoi[j-1].cos && pow(SorPoi[j].xaxis- SorPoi[A0].xaxis, 2)+pow(SorPoi[j].yaxis- SorPoi[A0].yaxis, 2) < pow(SorPoi[j-1].xaxis- SorPoi[A0].xaxis, 2)+pow(SorPoi[j-1].yaxis- SorPoi[A0].yaxis, 2) ) ) ){
		   ExchangeD(&SorPoi[j].xaxis, &SorPoi[j-1].xaxis);
	       ExchangeD(&SorPoi[j].yaxis, &SorPoi[j-1].yaxis);
	       ExchangeD(&SorPoi[j].cos, &SorPoi[j-1].cos);
		   j--;
		   }
		j++;		   
	 }
	 
	 //puts("sorted array");
	 //for (i=0; i<SizeSetX; i++)
	   //printf("%lf %lf\n", SorPoi[i].xaxis, SorPoi[i].yaxis);
	
	 
	 
	 
	//Find out the minmal cells necessary for the convex hull. Colinear cells are discarded. Use the structure of stack: last in, first out. Push and pop.
	ConHull[0]=SorPoi[0].xaxis; ConHull[0+LenConHull]=SorPoi[0].yaxis;
	ConHull[1]=SorPoi[1].xaxis; ConHull[1+LenConHull]=SorPoi[1].yaxis;
	j=2;
	for (i=2; i<SizeSetX; i++){
	  ConHull[j]=SorPoi[i].xaxis; ConHull[j+LenConHull]=SorPoi[i].yaxis;	
	  
	  while (j>=2&& ((ConHull[j-1]-ConHull[j-2])*(ConHull[j+LenConHull]- ConHull[j-1+LenConHull])-(ConHull[j-1+LenConHull]-ConHull[j-2+LenConHull])*(ConHull[j]- ConHull[j-1])<= 0) ){ //turn right or colinear
	    ConHull[j-1]=ConHull[j]; 
	    ConHull[j-1+LenConHull]=ConHull[j+LenConHull] ;
	    j--;
	    
	   // printf("discard one cell %lf %lf \n", ConHull[j], ConHull[j+LenConHull]);
	   }
	   j++; 
	}
	
return j;
}

int WheInConHull(double *ConHull, int NumPoints, int LenConHull, double x, double y){
	int i;
	//On the right or colinear with AN-1A0
	if ( CroPro2D(ConHull[0]-ConHull[NumPoints-1], ConHull[0+LenConHull]-ConHull[NumPoints-1+LenConHull], x- ConHull[0], y- ConHull[0+LenConHull]) <=0 )
		  return 0;
	for (i=1; i<NumPoints; i++){
		if ( (ConHull[i]-ConHull[i-1])*(y- ConHull[i+LenConHull])-(ConHull[i+LenConHull]-ConHull[i-1+LenConHull])*(x- ConHull[i])<= 0)   //On the right or colinear with Ai-1Ai
		  return 0;
	}
	
	return 1;
	}

int ExchangeD(double * a, double *b){
	double TemExc;
	TemExc= *a;
	*a= *b;
	*b=TemExc;
	return 0; }
	
double CroPro2D(double x1, double y1, double x2, double y2){
	return x1*y2-y1*x2;}

int MinArrXZ_I(int *array, int length){
	int k=0, mini;
	for(k=0; k<length; k++){
		//the above: minimal i
		if(array[k]< array[mini]) 
		  mini=k;
	 }
	 
	 return mini;	
	}

int MaxArrXZ_I(int *array, int length){
	int k=0, max=0;
	for(k=0; k<length; k++){
		//the above: minimal i
		if(array[k]> array[max]) 
		  max=k;
	 }
	 
	 return max;	
	}
	
double SumMat_D(double *source, int height, int width, int HeiMat, int WidMat, double threshold, char mode){
	int i, j;
	double sum=0;
	if (mode== 't'){  //With threshold
	  for(i=0; i<height; i++){
		  for(j=0; j<width; j++){
			  sum+= (source[i*WidMat+j]>= threshold ); 
		  }
	   }
	 }
	 
	else if (mode== 's'){  //Without threshold
	  for(i=0; i<height; i++){
		  for(j=0; j<width; j++){
			  sum+= source[i*WidMat+j]; 
		  }
	   }
	 }
	return sum;	
	}
