#include "C.h"
#include "PreWarping.h"

using namespace std;

double _C_( double x)
{	
	if ( x > C_fim ) return 0;
	if ( x < C_inicio) return 2.0;
		
	int n;
	n = round(((x-C_inicio)/(C_fim-C_inicio))*C_N);
	
	return _c_[n];
}
