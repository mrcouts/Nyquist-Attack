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
    void ComputeGain(float gaindB);
    void ApplyGain(vec *signal);

private:
    int N; // Tamanho dos vetores

    vec g; //Vetor de ganho
    float g0; //Condição inicial

    double a1;
    double b0;
    double wc;
};