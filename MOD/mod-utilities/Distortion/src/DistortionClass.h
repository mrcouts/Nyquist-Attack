#include <cmath>
#include <algorithm>
#include <armadillo>

using namespace std;
using namespace arma;

class DistortionClass
{
public:
    DistortionClass( double samplerate, int N );
    ~DistortionClass();
    void SetInput(float *input);
    void Wire();
    void OverWire();
    void Oversample8x();
    void Downsample8x();
    void CopyOutput(float *output);

    vec u; //Sinal de entrada
    vec y; //Sinal de saída

private:
    int N; // Tamanho dos vetores
    double SampleRate; //Frequência de amostragem em Hz
    double T; //Periodo de amostragem em segundos

    mat H; //Matrix de oversample
    mat h; //Downsample matrix
    vec ic; //Condições iniciais do oversample
    vec ic2;//Condições iniciais do downsample
    vec u_ic; // join_vert(ic,u)
    vec y_over_ic; // join_vert(ic2,y_over_ic)
    vec u_over; //u oversampled
    vec y_over; //y oversampled
};