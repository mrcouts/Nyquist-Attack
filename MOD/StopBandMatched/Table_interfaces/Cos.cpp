#include <cmath>
#include "Cosseno.h"

using namespace std;

double Cos( double x)
{	
	double aux;
	
	if ( (x>COS_fim)||(x<COS_inicio))
	{
		aux = ((x +M_PI)/(2*M_PI));
        aux = floor(aux);
        x = x-aux*(2*M_PI);
	}
	
	int n;
	n = round(((x-COS_inicio)/(COS_fim-COS_inicio))*COS_N);
	
	return Cosseno[n];
}
