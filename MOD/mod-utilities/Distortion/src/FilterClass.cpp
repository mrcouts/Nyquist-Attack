#include "FilterClass.h"
#include <iostream>

FilterClass::FilterClass( double samplerate, int N ) //Construtor da classe pai
{
	this->N = N;
	y.zeros(N);
	vwc.zeros(N);
	c.zeros(N);
	b0.zeros(N);
	b1.zeros(N);
	a0.zeros(N);
	a1.zeros(N);
	B0.zeros(N);
	B1.zeros(N);
	A0.zeros(N);
	A1.zeros(N);

	u_1 = 0;
	y_1 = 0;

	f_1 = 1;
	wc_1 = 2*M_PI*f_1;

	_wc_ = 8.0/N; //Frequencia de corte em rad/amostra

	_a1_ = -2 + cos(_wc_) + 0.5*sqrt(-4 +  pow(4-2*cos(_wc_),2) );
	_b0_ = 1 + _a1_;


	SampleRate = samplerate;
	T = 1/SampleRate;
}

FilterClass::~FilterClass()
{
	y.clear();
	vwc.clear();
	c.clear();
	b0.clear();
	b1.clear();
	a0.clear();
	a1.clear();
	B0.clear();
	B1.clear();
	A0.clear();
	A1.clear();
}

void FilterClass::ComputeCoef(float f)
{
	/*        wc
	G(s) = --------
	        s + wc
	*/

	wc = 2*M_PI*f;

	vwc(0) = -_a1_*wc_1 + _b0_*wc;
	c(0) = vwc(0)/tan(vwc(0)*T/2);
	for(int i = 1; i < N; i++)
	{
		vwc(i) = -_a1_*vwc(i-1) + _b0_*wc;
		c(i) = vwc(i)/tan(vwc(i)*T/2);
	}
	wc_1 = wc;

	b1.zeros();
	b0 = vwc;
	a1.ones();
	a0 = vwc;

	B0 = b0 + b1 % c;
	B1 = b0 - b1 % c;
	A0 = a0 + a1 % c;
	A1 = a0 - a1 % c;

	B0 = B0/A0;
	B1 = B1/A0;
	A1 = A1/A0;
}

void FilterClass::LPF1(double f, vec *u)
{
	if(f != f_1)
	{
		ComputeCoef(f);
		y(0) = -A1(0)*y_1 + B0(0)*u[0](0) + B1(0)*u_1;
		for (int i=1; i < N; i++) y(i) = -A1(i)*y(i-1) + B0(i)*u[0](i) + B1(i)*u[0](i-1);
		y_1 = y(N-1);
		u_1 = u[0](N-1);
	}
	else
	{
		y(0) = -A1(N-1)*y_1 + B0(N-1)*u[0](0) + B1(N-1)*u_1;
		for (int i=1; i < N; i++) y(i) = -A1(N-1)*y(i-1) + B0(N-1)*u[0](i) + B1(N-1)*u[0](i-1);
		y_1 = y(N-1);
		u_1 = u[0](N-1);
	}
	f_1 = f;
}