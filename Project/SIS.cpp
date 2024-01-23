#include <iostream>
#include <armadillo>
#include <complex>
#include <time.h>
#include <cmath>

using namespace std;
using namespace arma;
using namespace std::complex_literals;

int main() {
    bool debug = true;
    srand(time(NULL));

    int n = 4;
    int m = 3;
    int inst[m][3];

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < 3; j++) {
            int var = rand() % n + 1;
            int sign = pow(-1, rand() % 2);
            inst[i][j] = sign * var;
        }
    }

    if (debug) {
        for (int i = 0; i < m; i++) {
            cout << "Clause " << i << ": ";
            for (int j = 0; j < 3; j++) {
                cout << inst[i][j] << " ";
            }
            cout << endl;
        }
    }

    cx_mat X = {{0, 1.0 + 0i}, {1.0 + 0i, 0}};
    cx_mat Y = {{0, -1i}, {1i, 0}};
    cx_mat Z = {{1.0 + 0i, 0}, {0, -1.0 + 0i}};
    cx_mat I(2, 2, fill::eye);

    //Initialising Hamiltonians
    cx_mat HC(pow(2, n), pow(2, n));
    HC.fill(0.0 + 0.0i);

    cx_mat H_array[m];
    cx_mat H;
    cx_mat H_m_i = -(X+Z+I);
    cx_mat S[3] = {Z, X, I};

/*  cx_mat H_1 = (Z*I*I*I + I*X*I*I + I*I*I*Z + 2*I*I*I*I);
    cx_mat H_m_1 = -(I*X+I*Z+I*I);
     H_1.brief_print();
    H_m_1.brief_print();
    cx_mat H_F_1 = kron(H_1,kron(H_m_1,I));
    H_F_1.brief_print();

    cx_mat H_2 = (Z*I*I*I + I*I*Z*I + I*I*I*X + 2*I*I*I*I);
    cx_mat H_m_2 = -(I*X+I*Z+I*I);
    cx_mat H_F_2 = kron(H_2,kron(H_m_2,I));
    cx_mat H_F = H_F_1 + H_F_2;

    H_F.brief_print();*/
    cx_mat O = {{-2.0 + 0i, -1.0 + 0i}, {-1.0 + 0i, 0}};
    cx_mat P = {{1.0 + 0i, 2.0 + 0i},{ 2.0 + 0i, -1.0 + 0i}};
    cx_mat OP = O*P;
    OP.brief_print();


}