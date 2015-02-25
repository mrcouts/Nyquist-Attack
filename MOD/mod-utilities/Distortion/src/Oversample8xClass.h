#include <cmath>
#include <algorithm>
#include <armadillo>

using namespace std;
using namespace arma;

class Oversample8xClass
{
public:
    Oversample8xClass(int N);
    ~Oversample8xClass();
    void Oversample8x(vec *u);

    vec y; //Sinal de saída

private:
    int N; // Tamanho dos vetores

    mat H; //Matrix de oversample
    vec ic; //Condições iniciais do oversample
    vec u_ic; // join_vert(ic,u)
};