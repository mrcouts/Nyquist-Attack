#include <cmath>
#include "FilterClass.h"

FilterClass::FilterClass( double samplerate, int N ) //Construtor da classe pai
{
	this->N = N;
	u = new double[N];
	y = new double[N];

	u_1 = 0;
	u_2 = 0;

	y_1 = 0;
	y_2 = 0;

	wc_1 = 2*M_PI;

	SampleRate = samplerate;
	T = 1/SampleRate;
}

FilterClass::~FilterClass()
{
	delete[] u;
	delete[] y;
}

void FilterClass::LPF1_Bilinear(double f) //Método da classe pai
{
	/*        wc
	G(s) = --------
	        s + wc
	*/

	double w;

	wc = 2*M_PI*f;

	w = wc_1;
	c = w/tan(wc*T/2);
	b_1 = 0;
	b_0 = w;
	a_1 = 1;
	a_0 = w;
	B_0 = b_0 + b_1*c;
	B_1 = b_0 - b_1*c;
	A_0 = a_0 + a_1*c;
	A_1 = a_0 - a_1*c;

	y[0] = (-A_1*y_1 + B_0*u[0] + B_1*u_1)/A_0;

	for (int i=1; i < N; i++)
	{
		w = wc_1 + ((wc - wc_1)/(N-1))*(i);
		c = w/tan(wc*T/2);
		//b_1 = 0;
		b_0 = w;
		//a_1 = 1;
		a_0 = w;
		B_0 = b_0 + b_1*c;
		B_1 = b_0 - b_1*c;
		A_0 = a_0 + a_1*c;
		A_1 = a_0 - a_1*c;

		y[i] = (-A_1*y[i-1] + B_0*u[i] + B_1*u[i-1])/A_0;
	}

	wc_1 = wc;
	u_1 = u[N-1];
	y_1 = y[N-1];
}

BandFilterClass::BandFilterClass(double samplerate, int N):FilterClass(samplerate, N) //Construtor da classe filha
{
	bw_1 = M_PI;
}
