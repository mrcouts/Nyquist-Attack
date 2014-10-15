#include <cmath>

class FilterClass
{
public:
    FilterClass( double samplerate, int N );
    ~FilterClass();
    void LPF1_Bilinear(double f);
    void LPF2_Bilinear(double f);
    void LPF3_Bilinear(double f);

    //Variáveis:
    double SampleRate; //Frequência de amostragem em Hz
    double T; //Periodo de amostragem em segundos

    double *u; //Sinal de entrada
    double *u2; // u-->| |-->u2-->| |-->y
    double *y; //Sinal de saída
    int N; // Tamanho dos vetores

    //Condições iniciais:

    double u_1;
    double u_2;

    double u2_1;
    double u2_2;

    double y_1;
    double y_2;

    /*Coeficientes de Laplace:

			b_n s^n + ... b_1 s + b_0
	H(s) = ---------------------------
			a_n s^n + ... a_1 s + a_0
    */

    double b_0;
    double b_1;
    double b_2;

    double a_0;
    double a_1;
    double a_2;

    /*Coeficientes em Z:

			B_0 + B_1 z^{-1} + ... B_n z^{-n}
	H(z) = ----------------------------------
			A_0 + A_1 z^{-1} + ... A_n z^{-n}
    */

	double B_0;
    double B_1;
    double B_2;

    double A_0;
    double A_1;
    double A_2;

    double c; // Constante da transformação bilinear
    double c2; // c*c

    double wc; //Frequencia de corte ou central em rad/s
    double wc_1; //Frequencia de corte ou central em rad/s anterior
};

class BandFilterClass:public FilterClass
{
public:
	BandFilterClass(double samplerate, int N);

	double bw; //Largura de banda em rad/s = w2 - w1
	double bw_1; //Largura de banda em rad/s anterior

};