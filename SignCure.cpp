/* Othello D. Gomes
Sign Curing Numerics
Use c++14
*/
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
    cx_mat HC(pow(2, n + m), pow(2, n + m));
    HC.fill(0.0 + 0.0i);

    cx_mat H_m = - (X + Z + I);
    cx_mat S[3];
    cx_mat H;
    cx_mat H_k;
    cx_mat H_temp;
    for (int i = 0; i < m; i++){
        /*Construct H_C */
        /*order qubits: 0....n-1 and n....n+m-1 */
        for (int s = 0; s < 3; s++) {
            int variable = inst[i][s];
            if (abs(variable)== 1){
            if (variable < 0) {
                H_temp = X;
            }
            else if(variable > 0){
                H_temp = Z;
            }
            else {
                H_temp = I;
            }
        }
         if (abs(variable)== 2){
            if (variable < 0) {
                H_temp = X;
            }
            else if(variable > 0){
                H_temp = Z;
            }
            else {
                H_temp = I;
            }
        }
        if (abs(variable)== 3){
            if (variable < 0) {
                H_temp = X;
            }
            else if(variable > 0){
                H_temp = Z;
            }
            else {
                H_temp = I;
            }
        }
            if (abs(variable) < n){
            if (variable < 0) {
                H_temp = X;
            }
            else if(variable > 0){
                H_temp = Z;
            }
            else {
                H_temp = I;
            }
        }
        }
                if (inst[i][0]== -1) {
                    H_k = S[0];
                }
                else if (inst[i][1]== -1) {
                    H_k = S[1];
                }
                else if(inst[i][2]== -1) {
                    H_k = S[2];
                }
                else{
                    H_k = I;
                }
      }
    for (int j = 0; j < n; j++){ 
    for (int i = 0; i < m; i++){
                   if (inst[i][0] - 1 == j) {
                    H_k = kron(H_k,S[0]);
                }
                else if (inst[i][1] - 1 == j) {
                    H_k = kron(H_k,S[1]);
                }
                else if(inst[i][2] - 1 == j){
                    H_k = kron(H,S[2]);
                }
                else{
                    H_k = kron(H_k,I);
                }
    }
    }
    for (int j = 0; j < n; j++) {
        bool found_qubit = false;
        for (int i = 0; i < m; i++) {
            if (abs(inst[i][0]) - 1 == j || abs(inst[i][1]) - 1 == j || abs(inst[i][2]) - 1 == j) {
                found_qubit = true;
                break;
            }
        }
        if (found_qubit) {
            H_k = kron(H_k, H_m);
        }
        else {
            H_k = kron(H_k, I);
        }
        cx_mat H = HC + H_k;
    }
 H.print();
    return 0;
}
/*Intialise S for the if and else statements*/

/* Compute Avg. Sign <S>*/

//int eH = expmat(-H);

//int avgSign = (trace(eH)/trace(abs(eH)));
