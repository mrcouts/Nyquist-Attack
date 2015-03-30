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

double FilterClass::BilinearConstant(float wc, double T)
{
	return wc/tan(wc*T/2);
}

void FilterClass::SoftInput(float f)
{
	wc = 2*M_PI*f;

	vwc(0) = -_a1_*wc_1 + _b0_*wc;
	c(0) = BilinearConstant(vwc(0), T);
	for(int i = 1; i < N; i++)
	{
		vwc(i) = -_a1_*vwc(i-1) + _b0_*wc;
		c(i) = BilinearConstant(vwc(i), T);
	}
	wc_1 = wc;
}

void FilterClass::LPcoef()
{
	/*        wc
	G(s) = --------
	        s + wc
	*/

	b1.zeros();
	b0 = vwc;
	a1.ones();
	a0 = vwc;
}

void FilterClass::HPcoef()
{
	/*        s
	G(s) = --------
	        s + wc
	*/

	b1.ones();
	b0.zeros();
	a1.ones();
	a0 = vwc;
}

void FilterClass::Bilinear1()
{
	B0 = b0 + b1 % c;
	B1 = b0 - b1 % c;
	A0 = a0 + a1 % c;
	A1 = a0 - a1 % c;

	B0 = B0/A0;
	B1 = B1/A0;
	A1 = A1/A0;
}


void FilterClass::LPComputeCoef(float f)
{
	SoftInput(f);
	LPcoef();
	Bilinear1();
}

void FilterClass::HPComputeCoef(float f)
{
	SoftInput(f);
	HPcoef();
	Bilinear1();
}

void FilterClass::F1(double f, vec *u, bool lp)
{
	if(f != f_1)
	{
		if(lp) LPComputeCoef(f);
		else   HPComputeCoef(f);
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

void FilterClass::LPF1(double f, vec *u)
{
	F1(f, u, true);
}

void FilterClass::HPF1(double f, vec *u)
{
	F1(f, u, false);
}

FilterClass2::FilterClass2(double samplerate, int N):FilterClass(samplerate, N) //Construtor da classe filha
{
	b2.zeros(N);
	a2.zeros(N);
	B2.zeros(N);
	A2.zeros(N);

	u_2 = 0;
	y_2 = 0;
}

FilterClass2::~FilterClass2() //Destrutor da classe filha
{
	b2.clear();
	a2.clear();
	B2.clear();
	A2.clear();
}

void FilterClass2::LP2coef()
{
	/*                 wc²
	G(s) = ------------------------
	        s² + sqrt(2) wc s + wc²
	*/

	b2.zeros();
	b1.zeros();
	b0 = vwc % vwc;
	a2.ones();
	a1 = sqrt(2)*vwc;
	a0 = b0;
}

void FilterClass2::HP2coef()
{
	/*                 s²
	G(s) = ------------------------
	        s² + sqrt(2) wc s + wc²
	*/

	b2.ones();
	b1.zeros();
	b0.zeros();
	a2.ones();
	a1 = sqrt(2)*vwc;
	a0 = vwc % vwc;
}

void FilterClass2::F3coef()
{
	/*           wc s
	G(s) = -----------------
	        s² + wc s + wc²
	*/

	b2.zeros();
	b1 = vwc;
	b0.zeros();
	a2.ones();
	a1 = vwc;
	a0 = vwc%vwc;
}

void FilterClass2::Bilinear2()
{
	B0 =   b0 + b1%c +   b2%c%c;
	B1 = 2*b0        - 2*b2%c%c;
	B2 =   b0 - b1%c +   b2%c%c;
	A0 =   a0 + a1%c +   a2%c%c;
	A1 = 2*a0        - 2*a2%c%c;
	A2 =   a0 - a1%c +   a2%c%c;

	B0 = B0/A0;
	B1 = B1/A0;
	B2 = B2/A0;
	A1 = A1/A0;
	A2 = A2/A0;
}

void FilterClass2::LP2ComputeCoef(float f)
{
	SoftInput(f);
	LP2coef();
	Bilinear2();
}

void FilterClass2::HP2ComputeCoef(float f)
{
	SoftInput(f);
	HP2coef();
	Bilinear2();
}

void FilterClass2::F3ComputeCoef(float f)
{
	SoftInput(f);
	F3coef();
	Bilinear2();
}

void FilterClass2::F2(double f, vec *u, bool lp)
{
	if(f != f_1)
	{
		if(lp) LP2ComputeCoef(f);
		else   HP2ComputeCoef(f);
		y(0) = -A1(0)*y_1    - A2(0)*y_2    + B0(0)*u[0](0) + B1(0)*u_1       + B2(0)*u_2;
		y(1) = -A1(1)*y(0)   - A2(1)*y_1    + B0(1)*u[0](1) + B1(1)*u[0](0)   + B2(1)*u_1;
		for (int i=2; i < N; i++)
		y(i) = -A1(i)*y(i-1) - A2(i)*y(i-2) + B0(i)*u[0](i) + B1(i)*u[0](i-1) + B2(i)*u[0](i-2);
		y_1 = y(N-1);
		y_2 = y(N-2);
		u_1 = u[0](N-1);
		u_2 = u[0](N-2);
	}
	else
	{
		y(0) = -A1(N-1)*y_1    - A2(N-1)*y_2    + B0(N-1)*u[0](0) + B1(N-1)*u_1       + B2(N-1)*u_2;
		y(1) = -A1(N-1)*y(0)   - A2(N-1)*y_1    + B0(N-1)*u[0](1) + B1(N-1)*u[0](0)   + B2(N-1)*u_1;
		for (int i=2; i < N; i++)
		y(i) = -A1(N-1)*y(i-1) - A2(N-1)*y(i-2) + B0(N-1)*u[0](i) + B1(N-1)*u[0](i-1) + B2(N-1)*u[0](i-2);
		y_1 = y(N-1);
		y_2 = y(N-2);
		u_1 = u[0](N-1);
		u_2 = u[0](N-2);
	}
	f_1 = f;
}

void FilterClass2::LPF2(double f, vec *u)
{
	F2(f, u, true);
}

void FilterClass2::HPF2(double f, vec *u)
{
	F2(f, u, false);
}

FilterClass3::FilterClass3( double samplerate, int N ) //Construtor da classe pai
{
	this->N = N;
	y.zeros(N);
	
	filter1 = new FilterClass(samplerate, N);
	filter2 = new FilterClass2(samplerate, N);

	A1.zeros(2,N);
	A2.zeros(2,N);
	B0.zeros(2,N);
	B1.zeros(2,N);
	B2.zeros(2,N);

	Y.zeros(2,N);
	Y_1.zeros(2);
	Y_2.zeros(2);

	u_1 = 0;
	u_2 = 0;

	f_1 = 1;
}

FilterClass3::~FilterClass3()
{
	y.clear();
	delete filter1;
	delete filter2;
}

void FilterClass3::_P3ComputeCoef(float f, bool lp)
{
	if(lp) filter1->LPComputeCoef(f);
	else   filter1->HPComputeCoef(f);
	filter2->F3ComputeCoef(f);

	A1.row(0) = filter1->A1.t();
	A1.row(1) = filter2->A1.t();
	A2.row(1) = filter2->A2.t();
	B0.row(0) = filter1->B0.t();
	B0.row(1) = filter2->B0.t();
	B1.row(0) = filter1->B1.t();
	B1.row(1) = filter2->B1.t();
	B2.row(1) = filter2->B2.t();
}

void FilterClass3::F3(double f, vec *u, bool lp)
{
	if(f != f_1)
	{
		_P3ComputeCoef(f, lp);
		Y.col(0) = -A1.col(0)%Y_1        - A2.col(0)%Y_2        + B0.col(0)*u[0](0) + B1.col(0)*u_1       + B2.col(0)*u_2;
		Y.col(1) = -A1.col(1)%Y.col(0)   - A2.col(1)%Y_1        + B0.col(1)*u[0](1) + B1.col(1)*u[0](0)   + B2.col(1)*u_1;
		for (int i=2; i < N; i++)
		Y.col(i) = -A1.col(i)%Y.col(i-1) - A2.col(i)%Y.col(i-2) + B0.col(i)*u[0](i) + B1.col(i)*u[0](i-1) + B2.col(i)*u[0](i-2);
		y = Y.row(0).t() - Y.row(1).t();
		Y_1 = Y.col(N-1);
		Y_2 = Y.col(N-2);
		u_1 = u[0](N-1);
		u_2 = u[0](N-2);
	}
	else
	{
		Y.col(0) = -A1.col(N-1)%Y_1        - A2.col(N-1)%Y_2        + B0.col(N-1)*u[0](0) + B1.col(N-1)*u_1       + B2.col(N-1)*u_2;
		Y.col(1) = -A1.col(N-1)%Y.col(0)   - A2.col(N-1)%Y_1        + B0.col(N-1)*u[0](1) + B1.col(N-1)*u[0](0)   + B2.col(N-1)*u_1;
		for (int i=2; i < N; i++)
		Y.col(i) = -A1.col(N-1)%Y.col(i-1) - A2.col(N-1)%Y.col(i-2) + B0.col(N-1)*u[0](i) + B1.col(N-1)*u[0](i-1) + B2.col(N-1)*u[0](i-2);
		y = Y.row(0).t() - Y.row(1).t();
		Y_1 = Y.col(N-1);
		Y_2 = Y.col(N-2);
		u_1 = u[0](N-1);
		u_2 = u[0](N-2);
	}
	f_1 = f;
}

void FilterClass3::LPF3(double f, vec *u)
{
	F3(f, u, true);
}

void FilterClass3::HPF3(double f, vec *u)
{
	F3(f, u, false);
}