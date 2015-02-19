#include "DistortionClass.h"

DistortionClass::DistortionClass( double samplerate, int N ) //Construtor da classe pai
{
	this->N = N;
	u.zeros(N);
	y.zeros(N);

	H = mat(" 0.001882625446606  -0.002710089534866   0.000528221487907   0.116362366467962   0.014596275528642  -0.007805743254415   0.003894924369314  -0.001688091512466;  0.000312944378303   0.001702733514889  -0.010084104330975   0.110845064614024   0.031380665774813  -0.012937335688846   0.005639138698702  -0.002873541303131; -0.001276523364462   0.004971543960620  -0.016918850194126   0.100319469167197   0.049772267353074  -0.017279257145356   0.007023009804036  -0.001614950193983; -0.001720503037984   0.006921079864528  -0.020055571780079   0.085738234601382   0.068393924114215  -0.019939669133575   0.007519008110854  -0.002712118358723; -0.002712118358723   0.007519008110854  -0.019939669133575   0.068393924114215   0.085738234601382  -0.020055571780079   0.006921079864528  -0.001720503037984; -0.001614950193983   0.007023009804036  -0.017279257145356   0.049772267353074   0.100319469167197  -0.016918850194126   0.004971543960620  -0.001276523364462; -0.002873541303131   0.005639138698702  -0.012937335688846   0.031380665774813   0.110845064614024  -0.010084104330975   0.001702733514889   0.000312944378303; -0.001688091512466   0.003894924369314  -0.007805743254415   0.014596275528642   0.116362366467962   0.000528221487907  -0.002710089534866   0.001882625446606");
	h = mat("-0.037011176733667  -0.004861913297324  -0.003833440051794  -0.001918238641928   0.000733427210561   0.003853995956054   0.007107662327022   0.010082652982130   0.012364575556955   0.013567818702727   0.013381583718675   0.011647727671889   0.008384072434459   0.003821516187967  -0.001655560832785  -0.007641913014954  -0.013256069080979  -0.017883565055948  -0.020739173772553  -0.021157433579672  -0.018650352120071  -0.012910831441936  -0.003929134221630   0.008035405375639   0.022451291121282   0.038517749982585   0.055240820419681   0.071506679429859   0.086207335320331   0.098314342691984   0.106899716389173   0.111372317647399   0.111372317647399   0.106899716389173   0.098314342691984   0.086207335320331   0.071506679429859   0.055240820419681   0.038517749982585   0.022451291121282   0.008035405375639  -0.003929134221630  -0.012910831441936  -0.018650352120071  -0.021157433579672  -0.020739173772553  -0.017883565055948  -0.013256069080979  -0.007641913014954  -0.001655560832785   0.003821516187967   0.008384072434459   0.011647727671889   0.013381583718675   0.013567818702727   0.012364575556955   0.010082652982130   0.007107662327022   0.003853995956054   0.000733427210561  -0.001918238641928  -0.003833440051794  -0.004861913297324  -0.037011176733667");

	ic.zeros(7);
	ic2.zeros(56);
	u_ic.zeros(N+7);
	y_over_ic.zeros(8*N+56);
	u_over.zeros(8*N);
	y_over.zeros(8*N);

	SampleRate = samplerate;
	T = 1/SampleRate;
}

DistortionClass::~DistortionClass()
{
	u.clear();
	y.clear();
}

void DistortionClass::SetInput(float *input)
{
	for(int i = 0; i < N; i++) u(i) = input[i];
}

void DistortionClass::CopyOutput(float *output)
{
	for(int i = 0; i < N; i++) output[i] = y(i);
}

void DistortionClass::Wire()
{
	y = u;
}

void DistortionClass::OverWire()
{
	y_over = atan(100*u_over)/10;
}

void DistortionClass::Oversample8x()
{
	u_ic = 8*join_vert(ic,u);

	for(int i = 0; i < N;i++)
	{
		u_over.subvec(8*i,8*i+7) = H*u_ic.subvec(i,i+7);
	}
	ic = u.subvec(N-7,N-1);
}

void DistortionClass::Downsample8x()
{
	y_over_ic = join_vert(ic2,y_over);

	for(int i = 0; i < N;i++)
	{
		y.subvec(i,i) = h*y_over_ic.subvec(8*i,8*i+63);
	}
	ic2 = y_over.subvec(8*N-56,8*N-1);
}