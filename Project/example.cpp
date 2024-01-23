#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;
#include <iostream>
#include <armadillo>
#include <complex>
#include <time.h>
#include <cmath>
using namespace std;
using namespace arma;
using namespace std::complex_literals;

/* main function*/
int main() {
    bool debug = true;

    /*initialize random seed*/
    srand(time(NULL));

    /*Generate 3-SAT instances in this part of the code*/
    int m = 2;
    int n = 4;
    int inst[m][3];
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < 3; j++) {
            int var = rand() % n + 1;
            int sign = pow(-1, rand() % 2);
            inst[i][j] = sign * var;
        }
    }

    /* print out 3-sat instance*/
    if (debug) {
        for (int i = 0; i < m; i++) {
            cout << "Clause " << i << ": ";
            for (int j = 0; j < 3; j++) {
                cout << inst[i][j] << " ";
            }
            cout << endl;
        }
    }

    /*Generate Hamiltonians in this part of the code*/
    /*Pauli Matrices: need to make sure they're dense matrices*/
    cx_mat X = {{0, 1.0 + 0i}, {1.0 + 0i, 0}};
    cx_mat Y = {{0, -1i}, {1i, 0}};
    cx_mat Z = {{1.0 + 0i, 0}, {0, -1.0 + 0i}};
    cx_mat I(2, 2, fill::eye);

    //Initialising Hamiltonians
    cx_mat HC(pow(2, n), pow(2, n));
    HC.fill(0.0 + 0.0i);

    cx_mat H_array[m];
    int g = 3;
    cx_mat H_arr[n];
    cx_mat H;
    cx_mat H_m;
    cx_mat L[2];
    cx_mat S[g];
    cx_mat H_k;
    cx_mat H_temp;

    /*Generate the Hamiltonians for each clause*/
    for (int k = 0; k < m; k++) {
        H_m = 1.0 + 0i;
        for (int j = 0; j < 3; j++) {
            int var = abs(inst[k][j])-1;
            int sign = inst[k][j] > 0 ? 1 : 0;
            H = 0.5 + 0i;
            for (int i = 0; i < n; i++) {
                if (i == var) {
                    H = kron(H, S[sign]);
                } else {
                    H = kron(H, I);
                }
            }
            H_m = H_m * H;
        }
        H_array[k] = H_m;
        H_array[k].brief_print();
    }

    /*Add the Hamiltonians for each clause to get the full Hamiltonian*/
    for (int k = 0; k < m; k++) {
        HC += H_array[k];
    }

    /*Print the full Hamiltonian*/
    HC.brief_print();

    return 0;
}
