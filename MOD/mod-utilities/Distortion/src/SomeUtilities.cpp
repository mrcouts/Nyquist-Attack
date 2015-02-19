#include "SomeUtilities.h"

float VectorNorm1(float *vec, int N)
{
	float Norm = 0;
	for (int i = 0; i<N; i++)
		Norm += abs(vec[i]);
	return Norm;
}
