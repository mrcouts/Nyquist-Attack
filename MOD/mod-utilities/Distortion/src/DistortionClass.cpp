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

void DistortionClass::TgH(vec *u)
{
	y = tanh(u[0]);
}