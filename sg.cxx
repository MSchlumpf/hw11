#include <iostream>
#include <fstream>
#include <complex>
#include <sstream>
#include <cmath>
//-----------------------------------
using namespace std;
//-----------------------------------
typedef complex<double> cmplx;
//-----------------------------------
void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin);


void writeToFile(const cmplx* const v, const string s, const double dx,
         const int Nx, const double xmin, const double alpha,
         const double lambda, const double omega, const double t);


void step(cmplx* const psi0, const int Nx, const double xmin, const double dx, const double omega, const double dt);
//-----------------------------------
int main(){

	const int Nx = 300;
	const double xmin = -40;
	const double xmax = 40;
	const double Tend = 10*M_PI;
	const double dx = (xmax-xmin)/(Nx-1);
	const double dt = 0.1*dx;
	double t = 0;
	const int Na = 10;
	int Nk = int(Tend / Na / dt + 0.5);

	const double lambda = 10;
	const double omega = 0.2;
	const double alpha = sqrt(omega);

	stringstream strm;

	cmplx* psi0 = new cmplx[Nx];

	init(psi0, alpha, lambda, dx, dt, Nx, xmin);

	writeToFile(psi0,"psi_0", dx,Nx,xmin, alpha, lambda, omega,t);



	for (int i = 1; i <= Na; i++) {
		for (int j = 1; j <= Nk-1; j++) {
			step(psi0, Nx, xmin, dx, omega, dt);
			t+=dt;
		}
		strm.str("");
		strm << "psi_" << i;
		writeToFile(psi0,strm.str(), dx,Nx,xmin, alpha, lambda, omega,t);
	}
	cout << "t = " << t << endl;
	delete[] psi0;
	return 0;
}
//-----------------------------------

//-----------------------------------
void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin, const double alpha,
                 const double lambda, const double omega, const double t)
{
	ofstream out(s.c_str());
	double x, xi, xil;
	double h1, h2, h3;
	cmplx ana;
	for(int i=0; i<Nx; i++){
		x = xmin + i * dx;
		xi = alpha * x;
		xil = alpha * lambda;
		h1 = -0.5 * pow(xi - xil*cos(omega*t),2 );
		h2 = omega*t/2 + xi * xil * sin(omega*t);
		h3 =  - 0.25 * xil*xil* sin(2*omega*t);
		ana = cmplx( h1 , -h2 - h3  );
		ana = sqrt(alpha / sqrt(M_PI) ) * exp(ana);
		out << x << "\t" << norm(v[i]) << "\t" << v[i].real() << "\t" << v[i].imag()
		<< "\t" << norm(ana) << "\t" << ana.real() << "\t" << ana.imag() <<  endl;
	}
	out.close();
}
//-----------------------------------

void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin)
{
	const double x0 = dx*Nx * 0.5;
	for(int i=0;i<Nx; i++){
		double x = xmin + i*dx ;
		psi0[i] = sqrt(alpha/sqrt(M_PI)) * exp(- pow(alpha*(x-lambda),2)/2 );
	}
}
//-----------------------------------
void step(cmplx* const psi0, const int Nx, const double xmin, const double dx, const double omega, const double dt){
	//------------step1----------------
	double x = xmin;
	double k = pow(omega, 2);
	// calculate psi0[0]
	cmplx a = cmplx(0, dt/(4*dx*dx));
	cmplx dd = cmplx(1, -(dt/(2*dx*dx)+dt/4.0*k*x*x));
	cmplx temp1 = psi0[0];
	psi0[0] = dd*psi0[0] + a*psi0[1];
	//calculate psi0[1] to psi0[N-2]
	for(int i=1; i<Nx-1; i++){
		x += dx;
		dd = cmplx(1, -(dt/(2*dx*dx)+dt/4.0*k*x*x));
		cmplx temp2 = psi0[i];
		psi0[i] = a*temp1+dd*psi0[i]+a*psi0[i+1];
		temp1 = temp2;
	}
	// calculate psi0[N-1]
	x += dx;
	dd = cmplx(1, -(dt/(2*dx*dx)+dt/4.0*k*x*x));
	psi0[Nx-1] = a*temp1 + dd*psi0[Nx-1];
	
	//------------step2----------------
	a = cmplx(0, -dt/(4*dx*dx));
	cmplx* d = new cmplx[Nx];
	cmplx* u = new cmplx[Nx];
	cmplx* l = new cmplx[Nx];
	// define diagonal, lower and upper elements
	for(int i=0; i<Nx; i++){
		x = xmin +i*dx;
		d[i] = cmplx(1, dt/(2*dx*dx)+dt/4.0*k*x*x);
		u[i] = a;
		l[i] = a;
	}
	// forward substitution to make A upper trigagonal
	for(int i=0; i<Nx-1; i++){
		d[i+1]  -= l[i+1]/d[i]*u[i];
		psi0[i+1] -= l[i+1]/d[i]*psi0[i];
	}
	// solve with backward substitution
	psi0[Nx-1] = psi0[Nx-1]/d[Nx-1];
	for(int i=Nx-2; i>=0; i--){
		psi0[i] = (psi0[i]-u[i]*psi0[i+1])/d[i];
	}
	delete[] d;
	delete[] u;
	delete[] l;

}