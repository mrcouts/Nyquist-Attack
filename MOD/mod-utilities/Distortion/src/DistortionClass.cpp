#include "DistortionClass.h"

DistortionClass::DistortionClass(int N) //Construtor da classe pai
{
	this->N = N;
	y.zeros(8*N);
}

DistortionClass::~DistortionClass()
{
	y.clear();
}

void DistortionClass::SoftClip(vec *u)
{
	y = u[0]/(abs(u[0]) + ones(8*N));
}

void DistortionClass::ArcTg(vec *u)
{
	y = (2/M_PI)*atan(M_PI*u[0]/2);
}

void DistortionClass::TgH(vec *u)
{
	y = tanh(u[0]);
}  

void DistortionClass::HardClip(vec *u)
{
	for(int i = 0; i < 8*N; i++)
	{
		if(abs(u[0](i))<1) y(i) = u[0](i);
		else y(i) = (u[0](i)>0 ? 1: -1);
	}
}