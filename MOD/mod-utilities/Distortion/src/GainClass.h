#include <cmath>
#include <algorithm>
#include <armadillo>

using namespace std;
using namespace arma;

class GainClass
{
public:
    GainClass(int N);
    ~GainClass();
    void dBconv(float gaindB);
    void ComputeGain();
    void Gain(float gaindB, vec *u);

    vec y; //Sinal de saída

private:
    int N; // Tamanho dos vetores

    vec g; //Vetor de ganho
    float g0; //Condição inicial
    float gaindB0; //Entrada anterior
    float gain; // 10^(gaindB/20)

    //Variáveis auxiliares
    double a1;
    double b0;
    double wc;
};