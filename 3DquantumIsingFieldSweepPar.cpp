# include <cstdlib>
# include <fstream>
# include <iostream>
# include <cmath>
# include <ctime>
# include <string>
# include <omp.h>

using namespace std;

int main ( int argc, char *argv[] );

// ----------------------Glauber dynamics of Ising spins------------------
// to simulate magnetization reversal process in magnetic field sweeps in kamiokite:
// See	for details:			 https://arxiv.org/abs/2211.05028

// This cpp script is written by Evgenii Barts, University of Groningen
// It is based on the matlab version written by prof. Maxim Mostovoy 20-03-2022
// The adapted Matlab version can be found in Other_Codes folder

// Details can be found in PhD thesis "Unconventional magnetic states and defects" by Evgenii Barts

struct Observables{
	double Magnetization;
	double Polarization;
	double LB2Neel;
	double MagnetizationSTD;
	double PolarizationSTD;
	double LB2NeelSTD;
	double Acceptrate;
};

// Function prototypes
int *ising_3d_initialize ( int Nx, int Ny, int Nz, double thresh, int *seed );
bool *Zn_doping_3d_initialize ( int Nx, int Ny, int Nz, double xZn, int *seed );
void exchangefields(int Spin[], double HA[], int Nx, int Ny, int Nz, double J1, double Jperp, double mubAH,  double beta);
void r8mat_uniform_01 ( int Nx, int Ny, int Nz, int *seed, double r[] );
void r8mat_uniform_01a ( int NIT, int *seed, double r[] );
void r8mat_uniform_01b (int NIT, int *seed, double *r);
void timestamp ( );
double calMagnetization(int Nx, int Ny, int Nz, int *Spin, double HA[], bool ZnInd[],
		double muA, double muB, double muBohr, double eK34, double eK1);
double calPolarization(int Nx, int Ny, int Nz, int *Spin, double HA[], bool ZnInd[], double eK34, double eK1);
double calLB2Neel(int Nx, int Ny, int Nz, int *Spin);
		
void GlauberLoop(int Spin[], double HA[], bool ZnInd[], int Nx, int Ny, int Nz, int NIT, 
		double eK34, double eK1, double J12, double Jperp2, double ezee, int *seed, double rnd[], double *AR);
Observables GlauberLoopNMEASURES(int Spin[], double HA[], bool ZnInd[], int Nx, int Ny, int Nz, int NIT, int Nmeas0, int Nmeas,
			double mubBH, double mubAH, double muA, double muB, double muBohr, double J1, double Jperp,
			double K, double beta, int *seed);
		
double mean(double array[], int length);
double standardDeviation(double array[], int length);

void FScomputation(int iNsw);


//-----------------------------------------main----------------------
// Compile with:
// g++ -o FswPar -O3 -fopenmp 3DquantumIsingFieldSweepPar.cpp 
// g++ -o FswPar -O3 -fopenmp 3DquantumIsingFieldSweepPar.cpp -std=c++11 

int main ( int argc, char *argv[] )
{
	
	const int Nsw = 24;//4; // Number of independent sweeps
	
	// Start parallel loop
	#pragma omp parallel for
	for (int iNsw = 0; iNsw < Nsw; iNsw++){
		FScomputation(iNsw);
	}
	// End of parallel loop
	
	cout << "  End of execution.\n";
	cout << "\n";
	
	return 0;
}


// Perform a single field sweep experiment
void FScomputation(int iNsw){
		int i;
		int Nx, Ny, Nz; // Number of Ising sites
		int seed;
		double thresh = 0.5;
		
		int NIT = 1e7;//1e5;//2e5; // Number of iterations
		int Nmeas = 50; int Nmeas0 = 50;  // Number of measurements
		double J1 = 47.0;//49.0; // In-plane Heisenberg exhange coupling
		double H = 0.; // Applied mangetic field
		double muBohr = 0.67; //Bohr magneton in K/T
		double muABohr = 4.21; //muA (A sublattice) in Bohr magnetons
		double muBBohr = 4.51; //muB (B sublattice) in Bohr magnetons
		double muA = muABohr*muBohr;
		double muB = muBBohr*muBohr;
		double mubAH = muA*H; 
		double mubBH = muB*H;
		double Jperp = -0.2;//-1.2; // Interplane Heisenberg exhange coupling
		double K = 0.0;
		int NH = 64;//30; // Number of field points along the field sweep
		
		//initialize and set field points
		double *Hvalues;
		Hvalues = new double[NH];
		Hvalues[0] = -8.0;
		Hvalues[1] = -6.0;
		Hvalues[2] = -4.0;
		Hvalues[3] = -2.0;
		double Hfin, Hin, deltaH;
		Hin = 0.0001; Hfin = 8.0;
		int iStart = 4; int iEnd= NH ;
		deltaH = (Hfin - Hin)/(iEnd - iStart - 1);
		for (i=iStart; i<iEnd; i++){
			Hvalues[i] = Hin + deltaH*(i-iStart);
		}

		
		//cout << "  ETESTSTn.\n";
		
		double T, beta;
		T = 25.0; // Temperature
		beta = 1.0/T; // Inversed temperature
		
		// Observables: magnetization, electric polarization and B-sublattice Neel order parameter
		// and their standard deviations
		double *MagValues, *PolarizValues, *LB2Values;
		double *MagSTDValues, *PolarizSTDValues, *LB2STDValues;
		
		
		// Number of Ising sites
		Nx = 20; Ny = Nx ; Nz = Nx ; //iterations = 15; //seed = 123456789;
				
	
		int *Spin;
		bool *ZnInd;
		
		MagValues = new double[NH]; 
		PolarizValues = new double[NH];
		LB2Values = new double[NH];
		
		MagSTDValues = new double[NH]; 
		PolarizSTDValues = new double[NH];
		LB2STDValues = new double[NH];
		
		double *Acrate = new double[NH];
		
		//
		srand(unsigned(time(NULL))*(iNsw + 1));
		//cout<<(unsigned(time(NULL)))<<endl;
		//seed =  rand() + 123456789*iNsw;
		seed =  rand();
		//cout<< seed<<endl;
		//cout << "ISING_3D_ininital_state \n";
		
		
		//Uncomment to initialize random state
		//Spin = ising_3d_initialize(Nx, Ny, Nz, thresh, &seed );
		
		//Initialize all-down state
		Spin = new int[Nx*Ny*Nz];
		for ( i = 0; i < Nx*Ny*Nz; i++ ){ Spin[i] = -1; }
		
		// Initialize exchange fields
		double *HA; 
		HA = new double[Nx*Ny*Nz];
		
		// Introduce Zn doping 
		double xZn = 0.14; // Doping percentage
		ZnInd = Zn_doping_3d_initialize(Nx, Ny, Nz, xZn , &seed );
		
		Observables Coolingresults;	
		
		// Perform field sweep loop
		for (i=0; i<NH; i++){
			H = Hvalues[i]; // Current magnetic field
			mubAH = muA*H; 
			mubBH = muB*H;
			
			// Glauber loop
			Coolingresults = GlauberLoopNMEASURES(Spin, HA, ZnInd, Nx, Ny, Nz, NIT, Nmeas0, Nmeas, mubBH, mubAH,
							muA, muB, muBohr, J1, Jperp, K, beta, &seed);
			
			// Measure observables 
			PolarizValues[i] = Coolingresults.Polarization;
			MagValues[i] = Coolingresults.Magnetization;
			LB2Values[i] = Coolingresults.LB2Neel;
			//
			PolarizSTDValues[i] = Coolingresults.PolarizationSTD;
			MagSTDValues[i] = Coolingresults.MagnetizationSTD;
			LB2STDValues[i] = Coolingresults.LB2NeelSTD;
			//
			Acrate[i] = Coolingresults.Acceptrate;
			//cout << H << "	"<< LB2Values[i] << endl;
			//cout << Coolingresults.LB2Neel << endl;
		}
		
		// Save the data
		ofstream ofile;
		//string fileName = "FieldSweeps_" + to_string(iNsw) + ".txt";
		//NumberToString(T Number)
		//string fileName = "FieldSweeps_" + NumberToString(iNsw) + ".txt";
		//string fileName = "FieldSweeps_.txt";
		string fileName = "FSdata/FieldSweeps_" + to_string(iNsw) + ".txt";
		//ofile.open("CoolingResultsFS.txt");
		ofile.open(fileName);
		ofile << "FIELD" <<  "    "  << "MAGNETIZATION";
		ofile << "    "  << "MAGNETIZATIONstd" << "    " <<  "POLARIZATION";
		ofile << "    "  << "POLARIZATIONstd" << "    " <<  "LB2";
		ofile << "    "  << "LB2std" << "    " << "Accept rate"<< "    "<<  endl;
		for (i=0; i<NH; i++){
			ofile << Hvalues[i] <<  "    "  << MagValues[i]<< "    "  << MagSTDValues[i] << "    ";
			ofile << PolarizValues[i]<< "    "  << PolarizSTDValues[i] << "    ";
			ofile << LB2Values[i]<< "    "  << LB2STDValues[i] <<  "    " ;
			ofile << Acrate[i] << endl;
		}
		ofile.close();
		 
		delete [] Hvalues;
		delete [] Acrate;
		delete [] PolarizValues;
		delete [] MagValues;
		delete [] PolarizSTDValues;
		delete [] MagSTDValues;
		
		delete [] Spin;
		delete [] HA;
		delete [] ZnInd;
		//free ( Spin );
		//free ( HA );
		//free ( ZnInd );
}

// Glauber loop
void GlauberLoop(int Spin[], double HA[], bool ZnInd[], int Nx, int Ny, int Nz, int NIT, 
		double eK34, double eK1, double J12, double Jperp2, double ezee, int *seed, double rnd[], double *AR){	
	
	
	double accprob;
	int nx,ny,nz; 
	int myz, mxm, mxp, modnz, mym, myp, mzm, mzp;
	
	double hxm0, hxm1, hxp0, hxp1,hym0,hym1,hyp0,hyp1, hzm0, hzm1,hzp0,hzp1;
	double r;
	
	r8mat_uniform_01a (NIT, seed, rnd );
	
	double *xdice;
	xdice = new double[1];
	int kaccept = 0;
	for (int i = 0; i < NIT; i++)
	{
		r = 1.0;
		nx = int(Nx*rnd[i]); ny = int(Ny*rnd[i+NIT]); nz = int(Nz*rnd[i+2*NIT]);
		myz = (ny + nz) % 2;
		
		// along x
		mxm = nx - 1 + myz;
		if (mxm >= 0 && ZnInd[mxm + ny*Nx + nz*Nx*Ny]){
			hxm0 = HA[mxm + ny*Nx + nz*Nx*Ny];
			hxm1 = hxm0 + J12*Spin[nx+ny*Nx + nz*Nx*Ny];
			r = r*(cosh(hxm1)+cosh(0.5*hxm1)*eK34 + eK1)
				/(cosh(hxm0)+cosh(0.5*hxm0)*eK34 + eK1);
		}
		//
		mxp = nx  + myz;
		if (mxp < Nx  && ZnInd[mxp + ny*Nx + nz*Nx*Ny]){
			hxp0 = HA[mxp + ny*Nx + nz*Nx*Ny];
			hxp1 = hxp0 + J12*Spin[nx+ny*Nx + nz*Nx*Ny];
			r = r*(cosh(hxp1)+cosh(0.5*hxp1)*eK34 + eK1)
				 /(cosh(hxp0)+cosh(0.5*hxp0)*eK34 + eK1);
		}
		//along y
		modnz = nz % 2;
		if (modnz == 0){
			mym = ny - 1;
			if (mym >= 0 && ZnInd[nx + mym*Nx + nz*Nx*Ny]){
				hym0 = HA[nx + mym*Nx + nz*Nx*Ny];
				hym1 = hym0 + J12*Spin[nx+ny*Nx + nz*Nx*Ny];
				r = r*(cosh(hym1)+cosh(0.5*hym1)*eK34 + eK1)
					 /(cosh(hym0)+cosh(0.5*hym0)*eK34 + eK1);
			}
		} 
		else {
			myp = ny + 1;
			if (myp < Ny  && ZnInd[nx + myp*Nx + nz*Nx*Ny]){
				hyp0 = HA[nx + myp*Nx + nz*Nx*Ny];
				hyp1 = hyp0 + J12*Spin[nx+ny*Nx + nz*Nx*Ny];
				r = r*(cosh(hyp1)+cosh(0.5*hyp1)*eK34 + eK1)
					 /(cosh(hyp0)+cosh(0.5*hyp0)*eK34 + eK1);
			}
		}
		//along z
		mzm = nz - 1;
		if (mzm >= 0 && ZnInd[nx + ny*Nx + mzm*Nx*Ny]){
			hzm0 = HA[nx + ny*Nx + mzm*Nx*Ny];
			hzm1 = hzm0 + Jperp2*Spin[nx+ny*Nx + nz*Nx*Ny];
			r = r*(cosh(hzm1)+cosh(0.5*hzm1)*eK34 + eK1)
				 /(cosh(hzm0)+cosh(0.5*hzm0)*eK34 + eK1);
		}
		mzp = nz + 1;
		if (mzp < Nz  && ZnInd[nx + ny*Nx + mzp*Nx*Ny]){
			hzp0 = HA[nx + ny*Nx + mzp*Nx*Ny];
			hzp1 = hzp0 + Jperp2*Spin[nx+ny*Nx + nz*Nx*Ny];
			r = r*(cosh(hzp1)+cosh(0.5*hzp1)*eK34 + eK1)
				 /(cosh(hzp0)+cosh(0.5*hzp0)*eK34 + eK1);
		} 
		r = r*exp(-ezee*Spin[nx+ny*Nx + nz*Nx*Ny]);
		
		accprob = r / (1.0 + r);
		r8mat_uniform_01b(1, seed, xdice );

		if (xdice[0] <= accprob){
			Spin[nx+ny*Nx + nz*Nx*Ny] = -Spin[nx+ny*Nx + nz*Nx*Ny];
			kaccept++;
			
			if (mxm >= 0 && ZnInd[mxm+ny*Nx + nz*Nx*Ny]){
				HA[mxm+ny*Nx + nz*Nx*Ny]= hxm1;
			}
			if (mxp < Nx  && ZnInd[mxp+ny*Nx + nz*Nx*Ny]){
				HA[mxp+ny*Nx + nz*Nx*Ny] = hxp1;
			} 
			if (modnz == 0){
				if (mym >= 0 && ZnInd[nx+mym*Nx + nz*Nx*Ny]){
					HA[nx+mym*Nx + nz*Nx*Ny] = hym1;
				}
			}
			else{
				if(myp < Ny  && ZnInd[nx+myp*Nx + nz*Nx*Ny]){
					HA[nx+myp*Nx + nz*Nx*Ny] = hyp1;
				}
			}
			if (mzm >= 0 && ZnInd[nx+ny*Nx + mzm*Nx*Ny] ){
				HA[nx+ny*Nx + mzm*Nx*Ny] = hzm1;
			}
			if (mzp < Nz  && ZnInd[nx+ny*Nx + mzp*Nx*Ny] ){
				HA[nx+ny*Nx + mzp*Nx*Ny] = hzp1;
			}
		}
	}
	
	//cout << kaccept << endl;
	//AR[0] = float(kaccept)/NIT;
	*AR = float(kaccept)/NIT;
	//cout << *AR << endl;
	delete [] xdice;
	return;
}

// Perform several Glauber loops to accumulate statistically valuable averaged results
Observables GlauberLoopNMEASURES(int Spin[], double HA[], bool ZnInd[], int Nx, int Ny, int Nz, int NIT, int Nmeas0, int Nmeas,
			double mubBH, double mubAH, double muA, double muB, double muBohr, double J1, double Jperp, double K,
			double beta, int *seed){
	int i;
	double J12 = 2*beta*J1;
	double Jperp2 = 2*beta*Jperp;
	double ezee = 2*beta*mubBH;
	double eK1 = 0.5*exp(-beta*K);
	double eK34 = exp(-0.75*beta*K);

	Observables results;
	double *rnd;	rnd = new double[3*NIT];
	double* Mlist; double* LB2list; double* Pzlist; 
	Mlist = new double[Nmeas]; 
	LB2list = new double[Nmeas]; Pzlist = new double[Nmeas];
	double* ARlist= new double[Nmeas];
	double* AR = new double[Nmeas];
	
	exchangefields(Spin, HA, Nx, Ny, Nz, J1, Jperp, mubAH,  beta);
	for (i=0; i < Nmeas0; i++){
		GlauberLoop(Spin, HA, ZnInd, Nx, Ny, Nz, NIT,
					eK34,  eK1, J12, Jperp2, ezee, seed, rnd, AR);
	}
	
	for (i=0; i < Nmeas; i++){
		GlauberLoop(Spin, HA, ZnInd, Nx, Ny, Nz, NIT,
					eK34,  eK1, J12, Jperp2, ezee, seed, rnd, AR);
		
		Mlist[i]  =  calMagnetization(Nx, Ny, Nz, Spin, HA, ZnInd, muA, muB, muBohr, eK34, eK1);
		Pzlist[i] = calPolarization(Nx, Ny, Nz, Spin, HA, ZnInd, eK34, eK1);
		LB2list[i] = calLB2Neel(Nx, Ny, Nz, Spin);
		ARlist [i] = *AR; //AR[0]
	}
	results.Magnetization = mean(Mlist,Nmeas);
	results.Polarization  = mean(Pzlist,Nmeas);
	results.LB2Neel       = mean(LB2list,Nmeas);
	results.MagnetizationSTD = standardDeviation(Mlist,Nmeas); 
	results.PolarizationSTD  = standardDeviation(Pzlist,Nmeas) ;
	results.LB2NeelSTD 	     = standardDeviation(LB2list,Nmeas);
	results.Acceptrate    = mean(ARlist,Nmeas);
	delete [] Mlist;
	delete [] LB2list;
	delete [] Pzlist;
	
	delete [] AR;
	delete [] rnd;
	return results;
}

// Initialize random state
int *ising_3d_initialize ( int Nx, int Ny, int Nz, double thresh, int *seed )
{
	// (last element for k=0), max for i and j -> (Nx-1) +(Ny-1)*Nx= Ny*Nx -1
	// so start k=1 with Ny*Nx
	int *Spin;
	int i,j,k;
	double *r;
	r = new double[Nx*Ny*Nz];
	r8mat_uniform_01 ( Nx, Ny, Nz, seed, r );
	Spin = new int[Nx*Ny*Nz];
	for ( k = 0; k < Nz; k++ ) {for ( j = 0; j < Ny; j++ ) {for ( i = 0; i < Nx; i++ ){
		if ( r[i+j*Nx + k*Nx*Ny] <= thresh ){
			Spin[i+j*Nx + k*Nx*Ny] = -1;
		}
		else{
			Spin[i+j*Nx + k*Nx*Ny] = +1;
		}
	}}}
	delete [] r;
	return Spin;
}

// Initialize random Zn doping at the given concentration xZn
bool *Zn_doping_3d_initialize ( int Nx, int Ny, int Nz, double xZn, int *seed )
{
	// (last element for k=0), max for i and j -> (Nx-1) +(Ny-1)*Nx= Ny*Nx -1
	// so start k=1 with Ny*Nx
	bool *ZnInd;
	int i, iZn, NZn, NSITES;
	NSITES = Nx*Ny*Nz;
	ZnInd = new bool[NSITES];
	
	for (i=0; i<NSITES; i++){
		ZnInd[i] = true;
	}
	
	NZn = int(NSITES*xZn);
	double *xdice;
	xdice = new double[1];	
	iZn = 0;
	while(iZn<NZn){
		r8mat_uniform_01b (1, seed, xdice );
		i = int(xdice[0]*NSITES);
		if (ZnInd[i]){
			ZnInd[i] = false;
			iZn++;
		}
	}
	delete [] xdice;
	return ZnInd;
}

// Calculate exchange fields
void exchangefields(int Spin[], double HA[], int Nx, int Ny, int Nz, double J1, double Jperp, double mubAH, double beta){
	int i,j,k;
	int nx,ny,nz;
	int myz, mxm, mxp, modnz, mym, myp, mzm, mzp;
	double hxm0, hxm1, hxp0, hxp1,hym0,hym1,hyp0,hyp1, hzm0, hzm1,hzp0,hzp1;
	int NSITES;
	NSITES = Nx*Ny*Nz;
	
	for (i=0; i<NSITES; i++){
		HA[i] = mubAH;
	}
	
	for ( k = 0; k < Nz; k++ ) {for ( j = 0; j < Ny; j++ ) {for ( i = 0; i < Nx; i++ ){
		nx = i; ny = j; nz = k;
		myz = (ny + nz) % 2;
		
		// along x
		mxm = nx - 1 + myz;
		if (mxm >= 0){
			HA[mxm + ny*Nx + nz*Nx*Ny] += - J1*Spin[nx+ny*Nx + nz*Nx*Ny];
		}
		//
		mxp = nx  + myz;
		if (mxp < Nx ){
			HA[mxp + ny*Nx + nz*Nx*Ny] += - J1*Spin[nx+ny*Nx + nz*Nx*Ny];
		}
		//along y
		modnz = nz % 2;
		if (modnz == 0){
			mym = ny - 1;
			if (mym >= 0){
				HA[nx + mym*Nx + nz*Nx*Ny] += -J1*Spin[nx+ny*Nx + nz*Nx*Ny];
			}
		}
		else {
			myp = ny + 1;
			if (myp < Ny ){
				HA[nx + myp*Nx + nz*Nx*Ny] += -J1*Spin[nx+ny*Nx + nz*Nx*Ny];
			}
		}
		
		//along z
		mzm = nz - 1;
		if (mzm >= 0){
			HA[nx + ny*Nx + mzm*Nx*Ny] += -Jperp*Spin[nx+ny*Nx + nz*Nx*Ny];
		}
		mzp = nz + 1;
		if (mzp < Nz ){
			HA[nx + ny*Nx + mzp*Nx*Ny] += -Jperp*Spin[nx+ny*Nx + nz*Nx*Ny];
		}
    }}} 
	//for (i=0; i<20; i++){
	//	cout << HA[i] - mubAH<< endl;
	//}
	for (i=0; i<NSITES; i++){
		HA[i] = beta*HA[i];
	}
	//for (i=0; i<NSITES; i++){
	//	if(fabs(HA[i]) < 1e-5){
	//		HA[i] = 1e-6;
		//}
	//}
	//break;
	return;
}

// Calculate POLARIZATION
double calPolarization(int Nx, int Ny, int Nz, int *Spin, double HA[], bool ZnInd[],
						double eK34, double eK1){
	int i,j,k;
	int NSITES = Nx*Ny*Nz;
	double *B = new double[NSITES];
	
    double s = 0;
	for ( i = 0; i < NSITES; i++ ){
		if(ZnInd[i]){
			B[i] = (sinh(HA[i]) + 0.5*sinh(0.5*HA[i])*eK34)/(cosh(HA[i]) + cosh(0.5*HA[i])*eK34 + eK1);
		}
		else{
			B[i] = 0.0;
		}
	}
	
    for ( k = 0; k < Nz-1; k++ ) {for ( j = 0; j < Ny; j++ ) {for ( i = 0; i < Nx; i++ ){
		s = s + Spin[i+j*Nx + k*Nx*Ny]*B[i+j*Nx + (k+1)*Nx*Ny]
			+ B[i+j*Nx + k*Nx*Ny]*Spin[i+j*Nx + (k+1)*Nx*Ny];
    }}}
	delete []B;
    return s/float(2*Nx*Ny*(Nz-1));
}

// Calculate B-sublattice Neel vector squared
double calLB2Neel(int Nx, int Ny, int Nz, int *Spin){
	int i,j,k, nx, ny,nz;
    double s = 0;
    for ( k = 0; k < Nz; k++ ) {for ( j = 0; j < Ny; j++ ) {for ( i = 0; i < Nx; i++ ){
		nx = i; ny = j; nz = k;
		if(nz % 2 == 0){
			s = s + Spin[i+j*Nx + k*Nx*Ny];
		}
		else{
			s = s - Spin[i+j*Nx + k*Nx*Ny];
		}
    }}}
    return s/float(Nx*Ny*Nz);
}

// Calculate MAGNETIZATION
double calMagnetization(int Nx, int Ny, int Nz, int *Spin, double HA[], bool ZnInd[],
						double muA, double muB, double muBohr, double eK34, double eK1){
	int i,j,k, nx, ny,nz;
	double B;
    double s = 0;
	int NSITES = Nx*Ny*Nz;
	for ( i = 0; i < NSITES; i++ ){
		s = s + muB*Spin[i];
		if(ZnInd[i]){
			B = (sinh(HA[i]) + 0.5*sinh(0.5*HA[i])*eK34)/(cosh(HA[i]) + cosh(0.5*HA[i])*eK34 + eK1);
			s = s + muA*B ;
		}
	}
    return s/float(Nx*Ny*Nz)/muBohr;
}

double mean(double array[], int length){
    double s = 0;
    for (int i = 0; i < length; ++i)
        s = s + array[i];
    return s/float(length);
}

double standardDeviation(double array[], int length){
    double mu = mean(array,length);
    double s = 0;
    for (int i = 0; i < length; ++i)
        s = s + (array[i] - mu)*(array[i] - mu);
    s = s/float(length);
    return sqrt(s);
}

//****************************************************************************
void r8mat_uniform_01 ( int Nx, int Ny, int Nz, int *seed, double r[] ){
  int i, j, k;
  int i4_huge = 2147483647;
  int kseed;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R8MAT_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  for ( k = 0; k < Nz; k++ ) {for ( j = 0; j < Ny; j++ ) {for ( i = 0; i < Nx; i++ ){
      kseed = *seed / 127773;

      *seed = 16807 * ( *seed - kseed * 127773 ) - kseed * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }

      r[i+j*Nx + k*Nx*Ny] = ( double ) ( *seed ) * 4.656612875E-10;
    }}}
  return;
}

void r8mat_uniform_01a ( int NIT, int *seed,  double r[]){
  int i; 
  int i4_huge = 2147483647;
  int kseed;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R8MAT_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  for ( i = 0; i < 3*NIT; i++ ) {
      kseed = *seed / 127773;

      *seed = 16807 * ( *seed - kseed * 127773 ) - kseed * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }

      r[i] = ( double ) ( *seed ) * 4.656612875E-10;
    }

  return;
}

void r8mat_uniform_01b ( int NIT, int *seed,  double *r){
  int i; 
  int i4_huge = 2147483647;
  int kseed;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R8MAT_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  
      kseed = *seed / 127773;

      *seed = 16807 * ( *seed - kseed * 127773 ) - kseed * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }

      r[0] = ( double ) ( *seed ) * 4.656612875E-10;
    

  return;
}

void timestamp ( ){
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct std::tm *tm_ptr;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  //std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );
  std::strftime ( time_buffer, TIME_SIZE, "%I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
