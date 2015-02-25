#include <cmath>
#include <algorithm>
#include <armadillo>

using namespace std;
using namespace arma;

class Downsample8xClass
{
public:
    Downsample8xClass(int N);
    ~Downsample8xClass();
    void Downsample8x(vec *u);

    vec y; //Sinal de saída

private:
    int N; // Tamanho dos vetores

    mat h; //Downsample matrix
    vec ic; //Condições iniciais do downsample
    vec u_ic; // join_vert(ic,u)
};