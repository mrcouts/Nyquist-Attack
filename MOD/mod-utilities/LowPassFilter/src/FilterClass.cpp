#include "FilterClass.h"

FilterClass::FilterClass( double samplerate, int N ) //Construtor da classe pai
{
	this->N = N;
	u = new double[N];
	u2= new double[N];
	y = new double[N];

	u_1 = 0;
	u_2 = 0;

	u2_1 = 0;
	u2_2 = 0;

	y_1 = 0;
	y_2 = 0;

	wc_1 = 2*M_PI;

	SampleRate = samplerate;
	T = 1/SampleRate;
}

FilterClass::~FilterClass()
{
	delete[] u;
	delete[] u2;
	delete[] y;
}

void FilterClass::SetInput(float *in)
{
    copy_n(in, N, u);
}

void FilterClass::SetInput(double *in)
{
    copy_n(in, N, u);
}

void FilterClass::CopyOutput(float *out)
{
    copy_n(y, N, out);
}

void FilterClass::CopyOutput(double *out)
{
    copy_n(y, N, out);
}

bool FilterClass::SizeHasChanged(int N)
{
	return this->N != N;
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
	//c = w/tan(w*T/2);
	c = _C_(w*T)/T;

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
		c = _C_(w*T)/T;
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
	u_2 = u[N-2];
	y_2 = y[N-2];
}

void FilterClass::LPF2_Bilinear(double f) //Método da classe pai
{
	/*                 wc²
	G(s) = ------------------------
	        s² + sqrt(2) wc s + wc²
	*/

	double w;

	wc = 2*M_PI*f;

	w = wc_1;
	//c = w/tan(w*T/2);
	c = _C_(w*T)/T;
	c2 = c*c;

	b_2 = 0;
	b_1 = 0;
	b_0 = w*w;
	a_2 = 1;
	a_1 = sqrt(2)*w;
	a_0 = w*w;

	B_0 =   b_0 + b_1*c + b_2*c2;
	B_1 = 2*b_0       - 2*b_2*c2;
	B_2 =   b_0 - b_1*c + b_2*c2;
	A_0 =   a_0 + a_1*c + a_2*c2;
	A_1 = 2*a_0       - 2*a_2*c2;
	A_2 =   a_0 - a_1*c + a_2*c2;

	y[0] = (-A_1*y_1 - A_2*y_2 + B_0*u[0] + B_1*u_1 + B_2*u_2)/A_0;

	w = wc_1 + ((wc - wc_1)/(N-1))*(1);
	c = _C_(w*T)/T;
	c2 = c*c;

	//b_2 = 0;
	//b_1 = 0;
	b_0 = w*w;
	//a_2 = 1;
	a_1 = sqrt(2)*w;
	a_0 = w*w;

	B_0 =   b_0 + b_1*c + b_2*c2;
	B_1 = 2*b_0       - 2*b_2*c2;
	B_2 =   b_0 - b_1*c + b_2*c2;
	A_0 =   a_0 + a_1*c + a_2*c2;
	A_1 = 2*a_0       - 2*a_2*c2;
	A_2 =   a_0 - a_1*c + a_2*c2;

	y[1] = (-A_1*y[0] - A_2*y_1 + B_0*u[1] + B_1*u[0] + B_2*u_1)/A_0;

	for (int i=2; i < N; i++)
	{
		w = wc_1 + ((wc - wc_1)/(N-1))*(i);
		c = _C_(w*T)/T;
		c2 = c*c;

		//b_2 = 0;
		//b_1 = 0;
		b_0 = w*w;
		//a_2 = 1
		a_1 = sqrt(2)*w;
		a_0 = w*w;

		B_0 =   b_0 + b_1*c + b_2*c2;
		B_1 = 2*b_0       - 2*b_2*c2;
		B_2 =   b_0 - b_1*c + b_2*c2;
		A_0 =   a_0 + a_1*c + a_2*c2;
		A_1 = 2*a_0       - 2*a_2*c2;
		A_2 =   a_0 - a_1*c + a_2*c2;

		y[i] = (-A_1*y[i-1] - A_2*y[i-2] + B_0*u[i] + B_1*u[i-1] + B_2*u[i-2])/A_0;
	}

	wc_1 = wc;
	u_1 = u[N-1];
	y_1 = y[N-1];
	u_2 = u[N-2];
	y_2 = y[N-2];
}

void FilterClass::LPF3_Bilinear(double f) //Método da classe pai
{
	/*       wc           wc²
	G(s) = ------  ------------------
	       s + wc  s² +   wc s + wc²
	*/

	//Parte 1

	double w;

	wc = 2*M_PI*f;

	w = wc_1;
	//c = w/tan(w*T/2);
	c = _C_(w*T)/T;

	b_1 = 0;
	b_0 = w;
	a_1 = 1;
	a_0 = w;

	B_0 = b_0 + b_1*c;
	B_1 = b_0 - b_1*c;
	A_0 = a_0 + a_1*c;
	A_1 = a_0 - a_1*c;

	u2[0] = (-A_1*u2_1 + B_0*u[0] + B_1*u_1)/A_0;

	for (int i=1; i < N; i++)
	{
		w = wc_1 + ((wc - wc_1)/(N-1))*(i);
		c = _C_(w*T)/T;
		//b_1 = 0;
		b_0 = w;
		//a_1 = 1;
		a_0 = w;
		B_0 = b_0 + b_1*c;
		B_1 = b_0 - b_1*c;
		A_0 = a_0 + a_1*c;
		A_1 = a_0 - a_1*c;

		u2[i] = (-A_1*u2[i-1] + B_0*u[i] + B_1*u[i-1])/A_0;
	}

	//Parte 2

	w = wc_1;
	//c = w/tan(w*T/2);
	c = _C_(w*T)/T;
	c2 = c*c;

	b_2 = 0;
	b_1 = 0;
	b_0 = w*w;
	a_2 = 1;
	a_1 = w;
	a_0 = w*w;

	B_0 =   b_0 + b_1*c + b_2*c2;
	B_1 = 2*b_0       - 2*b_2*c2;
	B_2 =   b_0 - b_1*c + b_2*c2;
	A_0 =   a_0 + a_1*c + a_2*c2;
	A_1 = 2*a_0       - 2*a_2*c2;
	A_2 =   a_0 - a_1*c + a_2*c2;

	y[0] = (-A_1*y_1 - A_2*y_2 + B_0*u2[0] + B_1*u2_1 + B_2*u2_2)/A_0;

	w = wc_1 + ((wc - wc_1)/(N-1))*(1);
	c = _C_(w*T)/T;
	c2 = c*c;

	//b_2 = 0;
	//b_1 = 0;
	b_0 = w*w;
	//a_2 = 1;
	a_1 = w;
	a_0 = w*w;

	B_0 =   b_0 + b_1*c + b_2*c2;
	B_1 = 2*b_0       - 2*b_2*c2;
	B_2 =   b_0 - b_1*c + b_2*c2;
	A_0 =   a_0 + a_1*c + a_2*c2;
	A_1 = 2*a_0       - 2*a_2*c2;
	A_2 =   a_0 - a_1*c + a_2*c2;

	y[1] = (-A_1*y[0] - A_2*y_1 + B_0*u2[1] + B_1*u2[0] + B_2*u2_1)/A_0;

	for (int i=2; i < N; i++)
	{
		w = wc_1 + ((wc - wc_1)/(N-1))*(i);
		c = _C_(w*T)/T;
		c2 = c*c;

		//b_2 = 0;
		//b_1 = 0;
		b_0 = w*w;
		//a_2 = 1
		a_1 = w;
		a_0 = w*w;

		B_0 =   b_0 + b_1*c + b_2*c2;
		B_1 = 2*b_0       - 2*b_2*c2;
		B_2 =   b_0 - b_1*c + b_2*c2;
		A_0 =   a_0 + a_1*c + a_2*c2;
		A_1 = 2*a_0       - 2*a_2*c2;
		A_2 =   a_0 - a_1*c + a_2*c2;

		y[i] = (-A_1*y[i-1] - A_2*y[i-2] + B_0*u2[i] + B_1*u2[i-1] + B_2*u2[i-2])/A_0;
	}

	wc_1 =  wc;
	u_1  =  u[N-1];
	u2_1 = u2[N-1];
	y_1  =  y[N-1];
	u_2  =  u[N-2];
	u2_2 = u2[N-2];
	y_2  =  y[N-2];
}

BandFilterClass::BandFilterClass(double samplerate, int N):FilterClass(samplerate, N) //Construtor da classe filha
{
	bw_1 = M_PI;
}