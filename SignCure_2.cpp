/* Othello D. Gomes
UMD - Institute for Advanced Computer Studies
QuICS
05/15/2023
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
    cx_mat H;
    cx_mat H_n;
    cx_mat H_arr[n];
    cx_mat H_m_i = -(X+Z+I);
    cx_mat H_m;
    cx_mat S[3] = {Z, X, I};

    //This generates the clause Hamiltonian using the for loop for tensoring with I or Clause side of H_m.
    //The code is done backwards, meaning that the tensor conditionals are done first and the main for loop is done later on.
    //This makes no difference whatsoever, as in how the code is looped (forward or backward) will make no difference.
    for (int k = 0; k < m; k++){
    cx_mat H_i;
        H_i = 1.0 + 0i;
        for (int i = 0; i < n; i++) {         
            for(int idx=0; idx < 3; idx++){
                int curr_var = inst[k][idx];
                cout << curr_var << " Curr Var" << endl;
                if(abs(curr_var) == i && curr_var > k){
                    H_arr[i] = Z;
                }
                else if(abs(curr_var) == i && curr_var < k){
                    H_arr[i] = X;
                }
                else{
                    H_arr[i] = I;
                }
                H_n = H_arr[i];
                H_n.brief_print();
                }
                //write conditional statements based on curr_var to determine what Hamiltonian acts on qubit i. Will be one of I,X,Z
                //want to break out of this for loop if abs(curr_var)=i once you've set the Hamiltonian appropriately 
            //H_n = H_arr[i];
        cout << "Variable Hamiltonians." << endl;
        cout << "END OF Variable Hamiltonians." << endl;   
            //for each variable qubit decide if I'm acting with I, X, or Z based on the clause k
                    
        }
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
        H_m.brief_print();

        //for this clause k you need to generate the corresponding hamiltonian on variable qubits INSIDE THIS FOR LOOP. This clause has 3 non-trivial variables which will be acted on by either X or Z. All other variables are acted on by identity. Tensor these together then tensor the variable Hamiltonian with the clause Hamiltonian and THAT is one term of the full Hamiltonian which is a sum over clauses.

        //Structure is H=sum_{clauses} H(on clause qubits corresponding to clause) tensor H(on variable qubits corresponding to clause)
        
        //HC += H_i;
        //HC.print();
        //tensor H_i and H_m together to create term for this clause
        //add this tensor product to the full Hamiltonian. 
    }
}

/* Below is the code I've made to generate the variable Hamiltonians: so far, it has yielded some reliable results but the signs for the full Hamiltonian, which are tensored with the Clause Hamiltonians in the form of
kron(H_m,HC) appear to be negative -- which shouldn't be the case since the Full Hamiltonian is supposed to be invariant under conjugation.*/

//below here isn't right yet
// for (int j = 0; j < n; j++) { //this for loop doesn't make sense.
//     if (i == j) {
//         H_i = kron(H_i, Z);
//     } else {
//         H_i = kron(H_i, X);
//     }    
                                             

/*Intialise S for the if and else statements*/

/* Compute Avg. Sign <S>*/

//int eH = expmat(-H);

//int avgSign = (trace(eH)/trace(abs(eH)));
