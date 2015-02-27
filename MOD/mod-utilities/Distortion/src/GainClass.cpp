#include "GainClass.h"

GainClass::GainClass(int N) //Construtor da classe pai
{
	this->N = N;
	y.zeros(N);
	g.zeros(N);

	g0 = 1;
	gaindB0 = 0;
	wc = 6.0/N; //Frequencia de corte em rad/amostra

	a1 = -2 + cos(wc) + 0.5*sqrt(-4 +  pow(4-2*cos(wc),2) );
	b0 = 1 + a1;
}

GainClass::~GainClass()
{
	y.clear();
	g.clear();
}

void GainClass::dBconv(float gaindB)
{
	gain = pow(10.0,gaindB/20.0);
}

void GainClass::ComputeGain()
{
	g(0) = -a1*g0 + b0*gain;
	for(int i = 1; i < N; i++) g(i) = -a1*g(i-1) + b0*gain;
	g0 = g(N-1);
}

void GainClass::Gain(float gaindB, vec *u)
{
	dBconv(gaindB);
	if(gaindB != gaindB0)
	{
		ComputeGain();
		y = g % u[0];
	}
	else
	{
		g0 = gain;
		y = gain * u[0];
	}
	gaindB0 = gaindB;
}