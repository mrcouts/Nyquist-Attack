#include <cmath>
#include <iostream>
#include "MatchedFilters.h"

using namespace std;



void SBF1(double *u, double *y, int N, double f_before, double f_now, double BW_before, double BW_now, double T, double *U_1, double *U_2, double *Y_1, double *Y_2)
{
	double K;
	double A;
	double B;
	double f;
	double BW;
	double wc;
	double dw;
	double wcT;
	double dwT;

	double cos_wcT;
	
	double y_1 = Y_1[0];
	double u_1 = U_1[0];
	double y_2 = Y_2[0];
	double u_2 = U_2[0];
	
	/*
	              z² - 2cos(Twc)z + 1
	 G(z) = K ------------------------
	                z² + Az + B 
	*/
	
	f = f_before;
	BW = BW_before;
	wc = 2*M_PI*f;
	dw = 2*M_PI*BW;
	wcT = wc*T;
	dwT = dw*T;
	A = -2*exp(-dwT/2)*cos(0.5*sqrt(4*wcT*wcT - dwT*dwT) );
	B = exp(-dwT);
	cos_wcT = cos(wcT);
	K = (1+A+B)/(2*(1-cos_wcT));

	cout << "A =" << A << " B = " << B << " K = " << K <<  " cos_wcT = " << cos_wcT << " dwT = " << dwT << " wcT = " << wcT <<"\n" ;
	
	y[0] = -A*y_1 - B*y_2 + K*(u[0] -2*cos_wcT*u_1 + u_2);
	
	f = f_before + ((f_now - f_before)/(N-1))*(1);
	BW = BW_before + ((BW_now - BW_before)/(N-1))*(1);
	wc = 2*M_PI*f;
	dw = 2*M_PI*BW;
	wcT = wc*T;
	dwT = dw*T;
	A = -2*exp(-dwT/2)*cos(0.5*sqrt(4*wcT*wcT - dwT*dwT) );
	B = exp(-dwT);
	cos_wcT = cos(wcT);
	K = (1+A+B)/(2*(1-cos_wcT));
	
	y[1] = -A*y[0] - B*y_1 + K*(u[1] -2*cos_wcT*u[0] + u_1);
		
	for (int i=2; i<=N-1; i++)
	{
		f = f_before + ((f_now - f_before)/(N-1))*(i);
		BW = BW_before + ((BW_now - BW_before)/(N-1))*(i);
		wc = 2*M_PI*f;
		dw = 2*M_PI*BW;
		wcT = wc*T;
		dwT = dw*T;
		A = -2*exp(-dwT/2)*cos(0.5*sqrt(4*wcT*wcT - dwT*dwT) );
		B = exp(-dwT);
		cos_wcT = cos(wcT);
		K = (1+A+B)/(2*(1-cos_wcT));

		y[i] = -A*y[i-1] - B*y[i-2] + K*(u[i] -2*cos_wcT*u[i-1] + u[i-2]);
	}
	
	U_1[0] = u[N-1];
	Y_1[0] = y[N-1];
	U_2[0] = u[N-2];
	Y_2[0] = y[N-2];
	
}