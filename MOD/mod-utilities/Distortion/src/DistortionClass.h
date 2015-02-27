#include <cmath>
#include <algorithm>
#include <armadillo>

using namespace std;
using namespace arma;

class DistortionClass
{
public:
    DistortionClass(int N);
    ~DistortionClass();
    void TgH(vec *u);

    vec y; //Sinal de saída

private:
    int N; // Tamanho dos vetores
};