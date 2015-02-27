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
    void ComputeCoef(float f);
    void LPF1(double f, vec *u);

    vec y; //Sinal de saída

private:

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