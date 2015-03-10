#include <cmath>
#include <algorithm>
#include <armadillo>

using namespace std;
using namespace arma;

class FilterClass
{
public:
    FilterClass( double samplerate, int N );
    ~FilterClass();
    double BilinearConstant(float wc, double T);
    void SoftInput(float f);
    void LPcoef();
    void HPcoef();
    void Bilinear1();
    void LPComputeCoef(float f);
    void HPComputeCoef(float f);
    void LPF1(double f, vec *u);
    void HPF1(double f, vec *u);

    vec y; //Sinal de saída

//protected:

    int N; // Tamanho dos vetores
    double SampleRate; //Frequência de amostragem em Hz
    double T; //Periodo de amostragem em segundos

    //Condições iniciais:

    double u_1;
    double y_1;

    /*Coeficientes de Laplace:

			bn s^n + ... b1 s + b0
	H(s) = ------------------------
			an s^n + ... a1 s + a0
    */

    vec b0;
    vec b1;

    vec a0;
    vec a1;

    /*Coeficientes em Z:

			B_0 + B_1 z^{-1} + ... B_n z^{-n}
	H(z) = ----------------------------------
			A_0 + A_1 z^{-1} + ... A_n z^{-n}
    */

	vec B0;
    vec B1;

    vec A0;
    vec A1;

    vec c; // Constante da transformação bilinear

    double wc; //Frequencia de corte ou central em rad/s
    double wc_1; //Frequencia de corte ou central em rad/s anterior
    vec vwc;

    float f_1; //Entrada anterior

    //Variáveis auxiliares
    double _a1_;
    double _b0_;
    double _wc_;
};

class FilterClass2:public FilterClass
{
public:
    FilterClass2(double samplerate, int N);
    ~FilterClass2();
    void LP2coef();
    void HP2coef();
    void F3coef();
    void Bilinear2();
    void LP2ComputeCoef(float f);
    void HP2ComputeCoef(float f);
    void F3ComputeCoef(float f);
    void LPF2(double f, vec *u);
    void HPF2(double f, vec *u);

//protected:

    //Condições iniciais:
    double u_2;
    double y_2;

    //Coeficientes de Laplace:
    vec b2;
    vec a2;

    //Coeficientes em Z:
    vec B2;
    vec A2;
};

class FilterClass3
{
public:
    FilterClass3(double samplerate, int N);
    ~FilterClass3();
    void LP3ComputeCoef(float f);
    void LPF3(double f, vec *u);

    vec y; //Sinal de saída

//protected:
    int N; // Tamanho dos vetores

    FilterClass  *filter1;
    FilterClass2 *filter2;

    mat A1;
    mat A2;
    mat B0;
    mat B1;
    mat B2;

    mat Y;
    vec Y_1;
    vec Y_2;

    double u_1;
    double u_2;

    float f_1;
};