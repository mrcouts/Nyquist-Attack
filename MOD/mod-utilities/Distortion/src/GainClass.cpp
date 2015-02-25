#include "GainClass.h"

GainClass::GainClass(int N) //Construtor da classe pai
{
	this->N = N;
	y.zeros(N);
	g.zeros(N);

	g0 = 0;
	wc = 6.0/N; //Frequencia de corte em rad/amostra

	a1 = -2 + cos(wc) + 0.5*sqrt(-4 +  pow(4-2*cos(wc),2) );
	b0 = 1 + a1;
}

GainClass::~GainClass()
{
	g.clear();
}

void GainClass::ComputeGain(float gaindB)
{
	double gain = pow(10.0,gaindB/20.0);
	g(0) = -a1*g0 + b0*gain;
	for(int i = 1;i<N;i++) g(i) = -a1*g(i-1) + b0*gain;
	g0 = g(N-1);
}

void GainClass::ApplyGain(vec *u)
{
	y = g % u[0];
}

void GainClass::Gain(float gaindB, vec *u)
{
	ComputeGain(gaindB);
	ApplyGain(u);
}