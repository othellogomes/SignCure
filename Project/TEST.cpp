
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
    int n = 3;
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
  int dim = pow(2, n);
    cx_mat HC(dim, dim);
    HC.fill(0.0 + 0.0i);

    cx_mat H_array[m];
    cx_mat H;
    cx_mat H_m_i = -(X+Z+I);
    cx_mat H_m;
    cx_mat S[3] = {Z, X, I};

   for (int i = 0; i < n; i++) {
        cx_mat H_i;
        H_i = 1.0 + 0i;
        for (int j = 0; j < n; j++) {
            if (i == j) {
                H_i = kron(H_i, Z);
            } else {
                H_i = kron(H_i, X);
            }
        }
        cout << "Variable Hamiltonian" << endl;
        HC += H_i;
        HC.brief_print();
   }
    for (int k = 0; k < m; k++){
        for (int o = 0; o < m; o++){
            if (o == k){
                H_array[o] = H_m_i;
            }
            else{
                H_array[o] = I;
            }
        }
        H = H_array[m-1];
        for (int j = m-2; j >= 0; j--){
            H = kron(H_array[j], H);
        }
        H_m = H;
        cout << "Clause Hamiltonian" << endl; 
        H_m.brief_print();
    }

    cx_mat H_F = kron(H_m,HC);
    cout << "Full Hamiltonian" << endl;

    H_F.brief_print();
int eH = expmat(-H_F);

int avgSign = ((trace(eH))/trace(abs(eH)));

cout << eH << endl;
cout << avgSign << endl;
}

