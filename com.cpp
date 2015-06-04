//#include <stdio.h> /*an input-output library*/
//#include <stdlib.h> /*other functions*/
#include <math.h> /*math library*/
#include <string> /*String Function Library*/
#include <sstream>  //include this to use string streams
#include<cstring>
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <new>
#define MAXT 500
#define MAXTAU 50000
#define MAXAA 50
#define NMOL 2
#define MAXACOR 10 // Max number of molecules to calculate dipole-dipole correlation function
#define MAXBONDLENGTH 0.2
#define CUB(x) ((x)*(x)*(x))
#define SQR(x) ((x)*(x))
#define pi  3.14159
using namespace std;


// Classes
class Gromacs {

	private:
	int i,j;

	public:
		float x1,y1,z1; /*Centroid*/
		float cumsum;
		float VEC1[MAXAA],VEC2[MAXAA],VEC3[MAXAA],Distance[MAXAA-1]; // Vectors from to COM to each Atom in Gromacs  Distance is Vector Magnitude
		int adjgraph[MAXAA][MAXAA]; //Adjacency Graph 
		double DipDipCor[MAXTAU];

	

	// Calculate Weighted Average
	void Weighted_Ave(int Natoms,float O1X[MAXAA]) {		
		for(int i=0;i<Natoms;i++){
	   		cumsum+=O1X[i];
		}
		cumsum/=Natoms;	
	}
		
	// Calculate Center of Mass
	void COM(int Natoms, float O1X[MAXAA],float O1Y[MAXAA],float O1Z[MAXAA]){
	x1=y1=z1=0.0;	
		for(i=0;i<Natoms;i++){
	   		x1+=O1X[i];
			y1+=O1Y[i];
			z1+=O1Z[i];
		}

		// Centers of Mass and distance
		x1/=Natoms;
		y1/=Natoms;
		z1/=Natoms;
		for(i=0;i<Natoms;i++){
			VEC1[i] = O1X[i]-x1;
			VEC2[i] = O1Y[i]-y1;
			VEC3[i] = O1Z[i]-z1;
		}

		

	}

	// Generate Adjacency Graph
	void AGraph(float O1X[MAXAA],float O1Y[MAXAA],float O1Z[MAXAA],int Natoms) {

		//Initialize
		for(i=0;i<Natoms;i++){
			for(j=0;j<Natoms;j++){
				adjgraph[i][j]=0;
			}
		}
		//
		for(i=0;i<Natoms;i++){

			for(j=0;j<Natoms;j++){
			Distance[j] = 0;
			// First Calculate Pair-wise Distances
		if (i!=j){ 
		Distance[j] = sqrt((O1X[i]-O1X[j])*(O1X[i]-O1X[j])+(O1Y[i]-O1Y[j])*(O1Y[i]-O1Y[j])-(O1Z[i]-O1Z[j])*(O1Z[i]-O1Z[j]));	
			}
			}		

			for(j=0;j<Natoms;j++){
			if(i!=j && Distance[j]<=MAXBONDLENGTH) adjgraph[i][j]=1;
			}


		}

}

//end Function
/*
//Calculate Dipole-Dipole Correlation Function
void DDcor(double** mux[NMOL][MAXTAU], double** muy[NMOL][MAXTAU], double** muz[NMOL][MAXTAU]){
double Norm1, Norm2;
int tau;
	for(i=0;i<MAXTAU;i++) {
		DipDipCor[i] = 0;
	}

for(i=0;i<NMOL-1;i++){

	for(j=0;j<MAXTAU;j++){
		
		for(tau=0;tau<MAXTAU-j;tau++){
			Norm1 = sqrt(mux[i][j]*mux[i][j]+muy[i][j]*muy[i][j]+muz[i][j]*muz[i][j]);
			Norm2 = sqrt(mux[i][j+tau]*mux[i][j+tau]+muy[i][j+tau]*muy[i][j+tau]+muz[i][j+tau]);				
			DipDipCor[tau] += mux[i][j]*mux[i][j+tau]+muy[i][j+tau]*muy[i][j+tau]+muz[i][j]*muz[i+tau];
			DipDipCor[tau]/(Norm1*Norm2);
		}		
	
	DipDipCor[tau]/=(MAXTAU-j);

	}

}
*/

};

int main() {
// Read in Gromacs Coordinates
// 	Into
//	an	
//	Array
// string line;
float O1X[MAXAA],O1Y[MAXAA],O1Z[MAXAA]; /*Atom Cartesian Coordinates*/
float x1,y1,z1; /*Centroid*/
char Atom_type[2],ResID[8];
int Natoms=2;
int AtomID;

// Read in Dipole Moments
double** mux = new double*[NMOL];
double** muy = new double*[NMOL];
double** muz = new double*[NMOL];
for(int i = 0; i < NMOL; ++i){
    mux[i] = new double[MAXTAU];
    muy[i] = new double[MAXTAU];
    muz[i] = new double[MAXTAU];
}

//
// Grab dimer distances
//
ifstream myfile;
myfile.open("minus_minus_dist.xvg");
int i = 0;
double time,dx,dy,dz;
double distance[MAXTAU];
string line;
while ( i < MAXTAU) {
	getline (myfile,line);
	myfile >> time >> distance[i] >> dx >> dy >> dz ;
	i++;
    	}
    		
	myfile.close();
    
//
// Calculate radial distribution function
    int Length,MaxRadius,ir,Radius;
    float Bin_Width,Nfact,rho;
    Bin_Width = 0.05;
    MaxRadius = 6;
    Length = MaxRadius/Bin_Width; /*Total Numer of Distance Samples in Analysis*/
    double GoR[Length];
    for(i=0;i<Length;i++){
        GoR[i]=0.0;
    }
    
    for(i=0;i<MAXTAU;++i){
    ir = distance[i]/Bin_Width;
    GoR[ir] += 1.0;
    }
    
    ofstream outfile5;
    outfile5.open("minus_minus_rdf.txt");
    Nfact = 4*pi*MAXTAU*CUB(Bin_Width);
    // Ilans Rdf Code
    for (i = 0; i < Length; i++){
        Radius=i*Bin_Width;
        GoR[i] /= (Nfact*(i*(i+1)+1.0/3.0));
        outfile5 << Radius << "\t" << GoR[i] << "\n";
    }
    outfile5.close();

//
//
//	Grab Dipole Moment Output From
//	g_dipoles
//
int timestep;
Gromacs mac;
double totalmu;
for(int j = 0;j < NMOL;j++){
	std::string filename = "monomer";
	std::string fileopen;
    	char c1[100]="";
 	std::ostringstream ostr; //output string stream
    	ostr << j+1; 
    	std::string number = ostr.str(); //the str() function of the stream
    	fileopen = filename+number+"_dipole.txt";
    	strcpy(c1, fileopen.c_str());
	ifstream myfile;
	myfile.open(c1);
	int i = 0;
	string line;
	    while ( i < MAXTAU) {
		
	 	getline (myfile,line);
		myfile >> timestep >> mux[j][i] >> muy[j][i] >> muz[j][i]  >> totalmu;
		if(j==2) cout << "Fuck you" << "\n";
		i++;
    		}
    		
	myfile.close(); 
}

// Correlation Variables and Arrays
double DipDipCor[MAXTAU];
double AutoCor[MAXTAU], AutoCorSum[MAXTAU];
double Norm1, Norm2,theta_int;
int tau;
for(int i=0;i<MAXTAU;i++){
AutoCorSum[i] = DipDipCor[i]=0.0;
}

    
// Histogram Variables and Arrays
double bincostheta = 0.1;
int Nbins;
Nbins = 1+(1/bincostheta)*2; /*Number of costheta bins,plus one for 0 bin*/
int itheta;
double Prob[Nbins];
	for(int j=0;j<Nbins;j++){
		Prob[j]=0.0;
	}

    
//
// Calculate CrossCorrelation and Bin costheta of Dot Product
//

int probnorm;
int cornum;
for(int i=0;i<NMOL-1;i++){
probnorm = 0;
	for(tau=0;tau<MAXTAU;tau++){
	DipDipCor[tau] = 0;
	Norm1=0;
	Norm2=0;
	Norm1 = sqrt(mux[i][tau]*mux[i][tau]+muy[i][tau]*muy[i][tau]+muz[i][tau]*muz[i][tau]);
	Norm2 = sqrt(mux[i+1][tau]*mux[i+1][tau]+muy[i+1][tau]*muy[i+1][tau]+muz[i+1][tau]*muz[i+1][tau]);
	theta_int = (mux[i][tau]*mux[i+1][tau]+muy[i][tau]*muy[i+1][tau]+muz[i][tau]*muz[i+1][tau])/(Norm1*Norm2);
	itheta = ((Nbins-1)/2)+theta_int/bincostheta;
	if(distance[tau] <=1.0) { 
		Prob[itheta]+=1.0;
		probnorm +=1;
		}

	cornum = 0;
	for(int j=0;j<MAXTAU-tau;j++){

		if(distance[j] < 1.0 ){
		cornum +=1;
		Norm1 = sqrt(mux[i][j]*mux[i][j]+muy[i][j]*muy[i][j]+muz[i][j]*muz[i][j]);
		Norm2 = sqrt(mux[i+1][j+tau]*mux[i+1][j+tau]+muy[i+1][j+tau]*muy[i+1][j+tau]+muz[i+1][j+tau]*muz[i+1][j+tau]);		
		DipDipCor[tau]+=(mux[i][j]*mux[i+1][j+tau]+muy[i][j]*muy[i+1][j+tau]+muz[i][j]*muz[i+1][j+tau])/(Norm1*Norm2);
		}

	}		
	
	DipDipCor[tau]/=(cornum);
	
	}

}
    



//
// AutoCorrelation Sum
//
for(int i=0;i<NMOL;i++){

	for(tau=0;tau<MAXTAU;tau++){
	
		for(int j=0;j<MAXTAU-tau;j++){
		Norm1=0.0;
		Norm2=0.0;
		Norm1 = sqrt(mux[i][j]*mux[i][j]+muy[i][j]*muy[i][j]+muz[i][j]*muz[i][j]);
		Norm2 = sqrt(mux[i][j+tau]*mux[i][j+tau]+muy[i][j+tau]*muy[i][j+tau]+muz[i][j+tau]*muz[i][j+tau]);				
		AutoCorSum[tau]+=(mux[i][j]*mux[i][j+tau]+muy[i][j+tau]*muy[i][j+tau]+muz[i][j]*muz[i][j+tau])/(Norm1*Norm2);
		}		
	
	AutoCorSum[tau]/=(MAXTAU-tau);
	
	}
	
	

}

//
// Single Autocorrelation
//
for(int i=0;i<NMOL;i++){


	for(tau=0;tau<MAXTAU;tau++){
	AutoCor[tau] = 0;
	Norm1=0;
	Norm2=0;
		for(int j=0;j<MAXTAU-tau;j++){
		Norm1 = sqrt(mux[i][j]*mux[i][j]+muy[i][j]*muy[i][j]+muz[i][j]*muz[i][j]);
		Norm2 = sqrt(mux[i][j+tau]*mux[i][j+tau]+muy[i][j+tau]*muy[i][j+tau]+muz[i][j+tau]*muz[i][j+tau]);				
		AutoCor[tau]+=(mux[i][j]*mux[i][j+tau]+muy[i][j+tau]*muy[i][j+tau]+muz[i][j]*muz[i][j+tau])/(Norm1*Norm2);
		}		
	
	AutoCor[tau]/=(MAXTAU-tau);

	}


	// Print AutoCorrelations to File
    	std::string filename = "monomer";
	std::string fileopen;
	char c1[100]="";
 	std::ostringstream ostr; //output string stream
    	ostr << i+1; 
    	std::string number = ostr.str(); //the str() function of the stream
    	fileopen = filename+number+"_autocor.txt";
    	strcpy(c1, fileopen.c_str());
	ofstream outfile1;
	outfile1.open(c1);
	for(int i=0;i<MAXTAU;i++){
	outfile1 << i << "\t" << AutoCor[i] << "\n";
	}
	outfile1.close();
}



// Write Data to File
//
//
ofstream outfile1;
ofstream outfile2;
ofstream outfile3;
outfile1.open("minus_minus_cc.txt");
outfile2.open("prob_costheta.txt");
outfile3.open("averautocor.txt");
for(int i=0;i<MAXTAU;i++){
outfile1 << i << "\t" << DipDipCor[i] << "\n";
outfile3 << i << "\t" << AutoCorSum[i]/NMOL << "\n";
}
for(int i=0;i<Nbins;i++){
outfile2 << i*bincostheta << "\t" << Prob[i]/probnorm << "\n"; 
}


outfile1.close();
outfile2.close();
outfile3.close();



//Release Dynamically Allocated Memory
cout << "Warp 1 Engage" << "\n";
for(int i = 0; i < NMOL; ++i) {
delete [] mux[i];
delete [] muy[i];
delete [] muz[i];
}
delete [] mux;
delete [] muy;
delete [] muz;


}




	











