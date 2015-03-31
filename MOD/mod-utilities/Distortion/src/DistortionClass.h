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
    void SoftClip(vec *u);
    void ArcTg(vec *u);
    void TgH(vec *u);
    void HardClip(vec *u);

    vec y; //Sinal de saída

private:
    int N; // Tamanho dos vetores
};