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
            bool flag=true;
            int var;
            while(flag){
                var = rand() % n + 1;
                flag=false;
                for(int k=0; k<j; k++){
                    if(var==abs(inst[i][k])){
                        flag=true;
                    }
                }
            }
            int sign = pow(-1, rand() % 2);
            inst[i][j] = sign * var;
        }
    }

    /* print out 3-sat instance*/
    if (debug) {
        for (int i = 0; i < m; i++) {
            cout << "Clause " << i << ": ";    
            for (int k = 0; k < m; k++){               
                for (int j = 0; j < 3; j++) {
                    cout << inst[i][j] << " ";
                }
            cout << endl;
        }
    }
    }

    /* Generate Hamiltonians in this part of the code */
/* Pauli Matrices: need to make sure they're dense matrices */
cx_mat X = {{0, 1.0 + 0i}, {1.0 + 0i, 0}};
cx_mat Y = {{0, -1i}, {1i, 0}};
cx_mat Z = {{1.0 + 0i, 0}, {0, -1.0 + 0i}};
cx_mat I(2, 2, fill::eye);

// Initialize Hamiltonians
cx_mat HC(pow(2, n), pow(2, n));
HC.fill(0.0 + 0.0i);

cx_mat H_array[m];
cx_mat H_arr[3][n];
cx_mat H;
cx_mat H_m_i = -(X + Z + I);
cx_mat H_m;
cx_mat S[3]; // Variable Hamiltonians S1, S2, S3

// Generate the clause Hamiltonian using the for loop for tensoring with I or the clause side of H_m.
// The code is done backwards, meaning that the tensor conditionals are done first and the main for loop is done later on.
// This makes no difference whatsoever, as how the code is looped (forward or backward) will make no difference.

// Create Hamiltonians acting on variable qubits
 for (int k = 0; k < m; k++) {
for (int idx = 0; idx < 3; idx++) {
    H = 1.0 + 0i;
        // For each clause, decide if I'm acting with I, X, or Z based on the clause k.
        cx_mat H_n;
        int curr_var = inst[k][idx];
        for (int i = 0; i < n; i++) {
            if (abs(curr_var) == i + 1 && curr_var > 0) {
                H_arr[idx][i] = Z;
            } else if (abs(curr_var) == i + 1 && curr_var < 0) {
                H_arr[idx][i] = X;
            } else {
                H_arr[idx][i] = I;
            }
            H_n = H_arr[idx][i];
            H = kron(H, H_n);
        }
    S[idx] = H;
}
}
cx_mat S1 = S[0];
cout << "S1:" << endl;
S1.brief_print();
cx_mat S2 = S[1];
cout << "S2:" << endl;
S2.brief_print();
cx_mat S3 = S[2];
cout << "S3:" << endl;
S3.brief_print();
cout << "Variable Hamiltonians:" << endl;

cx_mat H_v = S1 + S2 + S3 + 2 * kron(I, kron(I,I)); 
H_v.brief_print();
return 0;
}