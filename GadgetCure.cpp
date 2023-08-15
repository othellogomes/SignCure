/* Othello D. Gomes
UMD - Institute for Advanced Computer Studies
QuICS
08/03/2023
Sign Curing Numerics w/ Gadgets
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
    int m = 1;
    int n = 3;
    int inst[m][3];
     //complex<double> avgSign(1.0, 0.0);
    //while (avgSign == complex<double>(1.0, 0.0)){
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
    cx_mat I_kron (pow(2,n),pow(2,n),fill::eye);

    //Initialising Hamiltonians
    cx_mat HC(pow(2, n+m), pow(2, n+m));
    HC.fill(0.0 + 0.0i);

    cx_mat H_array[m];
    cx_mat H_arr[3][n];
    cx_mat H_arra[n];
    cx_mat H;
    cx_mat H_m_i = -(X+Z+I);
    cx_mat H_m;
    cx_mat S[3];
    cx_mat H_1;
//cx_mat Test_H = kron(I,I);
//cx_mat Test_H2 = kron(Test_H, Z);
//Test_H2.brief_print();
//cout << "Test Kron Prod." << endl;
    //This generates the clause Hamiltonian using the for loop for tensoring with I or Clause side of H_m.
    //The code is done backwards, meaning that the tensor conditionals are done first and the main for loop is done later on.
    //This makes no difference whatsoever, as in how the code is looped (forward or backward) will make no difference.
    
    for (int k = 0; k < m; k++){
        // Create Hamiltonian acting on variable qubits
        for (int idx = 0; idx < 3; idx++) {
            H_1 = 1.0 + 0i;
            cx_mat H_n;
            // For each variable qubit, decide if I'm acting with I, X, or Z based on the clause k.
            int curr_var = inst[k][idx];
            for (int i = 0; i < n; i++) {
                if (abs(curr_var) == i+1 && curr_var > 0) {
                    H_arr[idx][i] = Z;
                } else if (abs(curr_var) == i+1 && curr_var < 0) {
                    H_arr[idx][i] = X;
                } else {
                    H_arr[idx][i] = I;
                }
                H_n = H_arr[idx][i];
                H_1 = kron(H_1, H_n);
            }

        S[idx] = H_1;
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

cx_mat H_v = S1 + S2 + S3 + 2 * I_kron; 
    
        //Create Hamiltonian acting on clause qubits
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
        H_m = H; // Clause Side of HC Tensor Product.
        H_m.brief_print();
        cout << "Clause Hamiltonians" << endl;
    cx_mat H_C = kron(H_m, H_v);
    H_C.brief_print();
    cout << "Full Hamiltonian" << endl; 
    HC += H_C;
    // Summation of Hamiltonian over m clauses.

    //Adding gadget
    cx_mat I_kr(pow(2,3*(n+m)),pow(2,3*(n+m)),fill::eye);
    cx_mat Hk = kron(HC, I_kr);
    Hk.brief_print();
   for (int u = 0; u < n+m-1; u++){
    int a_u = (n+m-1) + 3*u + 1;
    int b_u = (n+m-1) + 3*u + 2;
    int c_u = (n+m-1) + 3*u + 3;
    //I can now see why this approach would be messy....
    //Better to just create new identities and tensor them...
    cx_mat I_au(pow(2,a_u),pow(2,a_u));
    cx_mat I_bu(pow(2,b_u),pow(2,b_u));
    cx_mat I_cu(pow(2,c_u),pow(2,c_u));
    cx_mat I_u(pow(2,u),pow(2,u));

    cx_mat G1 =-kron((I_cu),kron(X,I_u))-kron(I_cu,kron(Z,I_u));
    cx_mat G2 =-kron(I,kron(X,kron(I-u,kron(X,I_au))));
    cx_mat G3 =-kron(I,kron(Y,kron(I-u,kron(Y,I_au))));
    cx_mat G4 =-kron(I,kron(Z,kron(I-u,kron(Z,I_au))));
    cx_mat G5 = -3*kron(I_au,kron(X,kron(X,I_bu)));
    cx_mat G6 = -kron(I_au,kron(Y,kron(Y,I_bu)));
    cx_mat G7 = -2*kron(I_au,kron(Z,kron(Z,I_bu)));
    cx_mat G8 = -kron(I_bu,kron(X,kron(X,I_cu)));
    cx_mat G9 = -kron(I_bu,kron(Y,kron(Y,I_cu)));
    cx_mat G10 = -kron(I_bu,kron(Z,kron(Z,I_cu)));

    cx_mat G = G1 + G2 + G3 + G4 + G5 + G6 + G7 + G8 + G9 + G10; 
    }
    //Create loops for the ancillary qubits corresponding to the matrices within the gadget.
    //This will definitely require conditionals similar to the variables of the Hamiltonian in line 97.
    cx_mat G_u[4];
    //For loop for the ancillary qubits, going from u, m+n+3(u-1) + 1.... m+n+3(u-1)+3
    for(int l = 0; l < u; l++){
    //probably have to add for loops for the ancillary qubits.
        G_u[0] = -(X+Z);
        G_u[1] = -(X*X + Y*Y + Z*Z);
        G_u[2] = -(3*X*X + Y*Y + 2*Z*Z);
        G_u[3] = -(X*X + Y*Y + Z*Z);
    }
}

    cout << "Summed H." << endl;

     //cx_mat H_2;
    cx_mat H_stoq = HC;
   for(int i = 0; i < n+m; i++){
    for (int j = 0; j < n+m; j++){
       // H_2 = -1.0 + 0i;
        //H_stoq.brief_print();
        if(i!=j){
            H_stoq(i,j) = -abs(HC(i,j));
        }
    }
   }
    cx_mat eH = expmat(-HC);
    complex<double> eH_1 = trace(eH);
    complex<double> eH_2 = trace(expmat(-H_stoq));
    complex<double> avgSign = eH_1/eH_2;
    cout << "Average Sign" << endl;
    cout << avgSign << endl;
  
   //cout << "Average Sign is not 1. " << endl;
   
    return 0;
}
// Add Gadgets - ancillary qubits and Hamiltonians in that Lemma!!!!

        //for this clause k you need to generate the corresponding hamiltonian on variable qubits INSIDE THIS FOR LOOP. This clause has 3 non-trivial variables which will be acted on by either X or Z. All other variables are acted on by identity. Tensor these together then tensor the variable Hamiltonian with the clause Hamiltonian and THAT is one term of the full Hamiltonian which is a sum over clauses.

        //Structure is H=sum_{clauses} H(on clause qubits corresponding to clause) tensor H(on variable qubits corresponding to clause)
        /*cx_mat H_i;
        H_i = 1.0 + 0i;
        for (int i = 0; i < n; i++) {            
            //for each variable qubit decide if I'm acting with I, X, or Z based on the clause k
                    for(int k = 0; k < n; i++){
                        if (k > i){
                            H_i = Z;
                        }
                        else if (k < i){
                            H_i = X;
                        }
                        else{
                            H_i = I;
                        }
                    }
           
        }*/
         //cout << "Variable Hamiltonian" << endl;
        //HC += H_i;
        //HC.print();
        //tensor H_i and H_m together to create term for this clause
        //add this tensor product to the full Hamiltonian. 

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
