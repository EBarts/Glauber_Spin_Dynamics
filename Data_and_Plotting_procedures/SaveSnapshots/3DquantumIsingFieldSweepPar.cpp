# include <cstdlib>
//# include <sstream>
# include <fstream>
# include <iostream>
# include <cmath>
# include <ctime>
# include <string>
//# include <omp.h>

using namespace std;

int main ( int argc, char *argv[] );

struct Observables{
	double Magnetization;
	double Polarization;
	double LB2Neel;
	double MagnetizationSTD;
	double PolarizationSTD;
	double LB2NeelSTD;
	double Acceptrate;
};

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
		double eK34, double eK1, double J12, double Jperp2, double ezee, int *seed, double rnd[],double *xdice, double *AR);
Observables GlauberLoopNMEASURES(int Spin[], double HA[], bool ZnInd[], int Nx, int Ny, int Nz, int NIT, int Nmeas0, int Nmeas,
			double mubBH, double mubAH, double muA, double muB, double muBohr, double J1, double Jperp,
			double K, double beta, int *seed, double *rnd);
		
double mean(double array[], int length);
double standardDeviation(double array[], int length);

void FScomputation(double beta);


//****************************************************************************80
// compile with g++ -o FswPar -O3 -fopenmp 3DquantumIsingFieldSweepPar.cpp 
//g++ -o Fsw -O3 3DquantumIsingFieldSweep.cpp -std=c++11 

int main ( int argc, char *argv[] )
{
	timestamp ( );
	double T, beta;
	T = 25.0; //30.0;//55;
	
	//cout << "T=" << T <<" input, 	";
	cout << "START SIMULATION" << endl;
	beta = 1.0/T;
	FScomputation(beta);
	cout << "  End of execution" << endl;
	timestamp ( );
	return 0;
}



void FScomputation(double beta){
		int i;
		int Nx, Ny, Nz;
		int seed;
		double thresh = 0.5;
		
		int NIT = 1e7;//1e6;//2e6;;//1e6;//5e5;//2e6;//1e5;//2e5;
		int Nmeas = 50; 
		int Nmeas0 = 50;
		double J1 = 47.0;//49.0;//49.0; 
		double H = 0.;
		double muBohr = 0.67; //in K/T
		double muABohr = 4.21; //muA in Bohr magnetons
		double muBBohr = 4.51; //muB in Bohr magnetons
		double muA = muABohr*muBohr;
		double muB = muBBohr*muBohr;
		double mubAH = muA*H; 
		double mubBH = muB*H;
		double Jperp = -0.2; //-1.2;//1.6;//-0.4;//-0.4;
		double K = 0.0;
		int NH = 32;//30;
		
		
		double *Hvalues;
		Hvalues = new double[NH];
		
		
		//cout << "  ETESTSTn.\n";
		//#pragma omp for
		//int HOtherValues = 4;
		
		//Hvalues[0] = -6.0;
		//Hvalues[1] = -4.0;
		//Hvalues[2] = -2.0;
		//Hvalues[NH-1] = 9.0;
		double Hfin, Hin, deltaH;
		//Hin = 0.0; Hfin = 7.0;
		//Hin = 0.0; Hfin = 10.0;
		
		//int iStart = 3; int iEnd= NH - 1;
		//deltaH = (Hfin - Hin)/(iEnd - iStart - 1);
		
		//for (i=0; i<NH; i++){
		//for (i=iStart; i<iEnd; i++){
			//Hvalues[i] = Hin + deltaH*(i-iStart);
		//}
		
		//Hin = 0.1; Hfin = 1.0;//8.0;
		//Hin = -2.0; Hfin = 8.0;//8.0;
		Hin = 0.0001;//0.001;//0.0001; 
		Hfin = 10.0;//10.0;//10.0;//8.0;
		deltaH = (Hfin - Hin)/(NH  - 1);
		for (i=0; i<NH; i++){
			Hvalues[i] = Hin + deltaH*i;
		}
		
		//cout << "  ETESTSTn.\n";
		
		
		
		double xZn = 0.14;//0.14;
		//Nx = 16; Ny = Nx; Nz = Nx; //iterations = 15; //seed = 123456789;
		//Nx = 8; Ny = 8; Nz = 8;//16;
		//Nx = 128; Ny = 128; Nz = 2; // 2^7 x 2^7 x 2 (power 15)
		
		//Nx = 8; Ny = 8; Nz = 512;//Nx = 2^3; Ny = 2^3; Nz = 2^9;
		//Nx = 256; Ny = 256; Nz = 256;
		Nx = 32; Ny = Nx; Nz = Nx;
		// Ising cooling
		double *MagValues, *PolarizValues, *LB2Values;
		double *MagSTDValues, *PolarizSTDValues, *LB2STDValues;
		//
	
		int *Spin;
		bool *ZnInd;
		//cout << Spin << endl;
		
		MagValues = new double[NH]; 
		PolarizValues = new double[NH];
		LB2Values = new double[NH];
		
		MagSTDValues = new double[NH]; 
		PolarizSTDValues = new double[NH];
		LB2STDValues = new double[NH];
		
		double *Acrate = new double[NH];
		//
		srand(unsigned(time(NULL)));
		//cout<<(unsigned(time(NULL)))<<endl;
		//seed =  rand() + 123456789*iNsw;
		seed =  rand();
		//cout<< seed<<endl;
		//cout << "ISING_3D_ininital_state \n";
		//random state
		//Spin = ising_3d_initialize(Nx, Ny, Nz, thresh, &seed );
		
		//All down state
		Spin = new int[Nx*Ny*Nz];
		for ( i = 0; i < Nx*Ny*Nz; i++ ){ Spin[i] = -1; }
		// AFM state
		//for ( int k = 0; k < Nz; k++ ) {for ( int j = 0; j < Ny; j++ ) {for ( i = 0; i < Nx; i++ ){
		//	if(k%2 == 0){
		//		Spin[i+j*Nx + k*Nx*Ny] = 1;
		//	}
		//	else{
		//		Spin[i+j*Nx + k*Nx*Ny] = -1;
		//	}
		//}}}
		
		//cout << iNsw <<"	"<< Spin[0] <<"		"<< Spin[5]<<"	" << Spin[27] <<endl;
		double *HA; 
		HA = new double[Nx*Ny*Nz];
		
		
		
		ZnInd = Zn_doping_3d_initialize(Nx, Ny, Nz, xZn , &seed );
		
		Observables Coolingresults;	
		double *rnd;	rnd = new double[3*NIT];
		for (i=0; i<NH; i++){
			H = Hvalues[i];// cout << T << endl;
			mubAH = muA*H; 
			mubBH = muB*H;
			
			Coolingresults = GlauberLoopNMEASURES(Spin, HA, ZnInd, Nx, Ny, Nz, NIT, Nmeas0, Nmeas, mubBH, mubAH,
							muA, muB, muBohr, J1, Jperp, K, beta, &seed, rnd);
			PolarizValues[i] = Coolingresults.Polarization;
			MagValues[i] = Coolingresults.Magnetization;
			LB2Values[i] = Coolingresults.LB2Neel;
			//
			PolarizSTDValues[i] = Coolingresults.PolarizationSTD;
			MagSTDValues[i] = Coolingresults.MagnetizationSTD;
			LB2STDValues[i] = Coolingresults.LB2NeelSTD;
			//
			Acrate[i] = Coolingresults.Acceptrate;
			cout << H  << endl;
			//cout << H << "	"<< LB2Values[i] << endl;
			//cout << Coolingresults.LB2Neel << endl;
		}
		ofstream ofile;
		//string fileName = "FieldSweeps_" + to_string(iNsw) + ".txt";
		//NumberToString(T Number)
		//string fileName = "FieldSweeps_" + NumberToString(iNsw) + ".txt";
		string fileName = "FieldSweepSingle.txt";
		//string fileName = "FSdata\\FieldSweeps_" + to_string(iNsw) + ".txt";
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
		delete [] rnd;
		
}

void GlauberLoop(int Spin[], double HA[], bool ZnInd[], int Nx, int Ny, int Nz, int NIT, 
		double eK34, double eK1, double J12, double Jperp2, double ezee, int *seed, double rnd[], double *xdice, double *AR){	
	
	
	double accprob;
	int nx,ny,nz; 
	int myz, mxm, mxp, modnz, mym, myp, mzm, mzp;
	
	double hxm0, hxm1, hxp0, hxp1,hym0,hym1,hyp0,hyp1, hzm0, hzm1,hzp0,hzp1;
	double r;
	
	r8mat_uniform_01a (NIT, seed, rnd );
	
	//cout << "    rnd[0] "  << rnd[0]<< endl;
	//cout << "    rnd[1] "  << rnd[1]<< endl;
	//cout << "    rnd[2] "  << rnd[2]<< endl;
	//cout << "    rnd[3] "  << rnd[3]<< endl;
	int kaccept = 0;
	for (int i = 0; i < NIT; i++)
	{
		r = 1.0;
		nx = int(Nx*rnd[i]); ny = int(Ny*rnd[i+NIT]); nz = int(Nz*rnd[i+2*NIT]);
		myz = (ny + nz) % 2;
		//if (i==0){
			//cout << "    rnd[0] "  << rnd[i]<< endl;
			//cout << "    nx="  << nx << endl;
		//}
		//	cout << "    nx="  << nx;
		//	cout << "    rnd[i] "  << rnd[i]<< endl;
		//	cout << "    ny="  << ny;
		//	cout << "    nz="  << nz <<endl;
		
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
		//cout << xdice[0] << endl;
		//cout <<" xdice is " << xdice[0] << "	";
		//cout <<" r is " << r << endl;
		if (*xdice <= accprob){
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
	return;
}


Observables GlauberLoopNMEASURES(int Spin[], double HA[], bool ZnInd[], int Nx, int Ny, int Nz, int NIT, int Nmeas0, int Nmeas,
			double mubBH, double mubAH, double muA, double muB, double muBohr, double J1, double Jperp, double K,
			double beta, int *seed, double *rnd){
	int i;
	double J12 = 2*beta*J1;
	double Jperp2 = 2*beta*Jperp;
	double ezee = 2*beta*mubBH;
	double eK1 = 0.5*exp(-beta*K);
	double eK34 = exp(-0.75*beta*K);

	Observables results;
	
	double* Mlist; double* LB2list; double* Pzlist; 
	Mlist = new double[Nmeas]; 
	LB2list = new double[Nmeas]; Pzlist = new double[Nmeas];
	double* ARlist= new double[Nmeas];
	double* AR = new double[Nmeas];
	double *xdice;
	xdice = new double[1];	
	
	exchangefields(Spin, HA, Nx, Ny, Nz, J1, Jperp, mubAH,  beta);
	for (i=0; i < Nmeas0; i++){
		GlauberLoop(Spin, HA, ZnInd, Nx, Ny, Nz, NIT,
					eK34,  eK1, J12, Jperp2, ezee, seed, rnd,xdice, AR);
	}
	
	
	for (i=0; i < Nmeas; i++){
		GlauberLoop(Spin, HA, ZnInd, Nx, Ny, Nz, NIT,
					eK34,  eK1, J12, Jperp2, ezee, seed, rnd,xdice, AR);
		
		Mlist[i]  =  calMagnetization(Nx, Ny, Nz, Spin, HA, ZnInd, muA, muB, muBohr, eK34, eK1);
		Pzlist[i] = calPolarization(Nx, Ny, Nz, Spin, HA, ZnInd, eK34, eK1);
		LB2list[i] = calLB2Neel(Nx, Ny, Nz, Spin);
		ARlist [i] = *AR; //AR[0]
		
	}
	
	ofstream ofile;
	//string fileName = "FSdataSpinConf\\ConfigFS_" + to_string(mubBH/muA) + ".txt";
	string fileName = "ConfigFS.txt";
	ofile.open(fileName, ios_base::app);
	int NSITES = Nx*Ny*Nz;
	
	for ( int iSpin = 0; iSpin < NSITES; iSpin++ ){
		bool stateif = Spin[iSpin] == 1;
		ofile << stateif <<  ",";
		//ofile << stateif <<  "    ";
	}
	ofile << endl;
	//ofile << " END OF SPIN LIST" << endl;
	ofile.close();
	
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
	delete [] xdice;
	
	return results;
}


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

//****************************************************************************
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
