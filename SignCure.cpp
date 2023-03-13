/* Othello D. Gomes
Sign Curing Numerics

Use c++14

*/
#include <iostream>
#include <armadillo> 
#include <complex>  
using namespace std;
using namespace arma;
using namespace std::complex_literals;

/*Generate 3-SAT instances in this part of the code*/

int main()
{
int m = 2; 
int n = 4;  
int inst [m][3];
for (int i = 0; i < m-1; i++){
  for (int j = 0; j <= 2; j++){
    int var = rand() % n + 1;
    int sign = rand() % 4 - 4;
    inst [i][j] = sign * var;
  }
}

/*Generate Hamiltonians in this part of the code*/
/*Pauli Matrices: need to make sure they're dense matrices*/

cx_mat X = {{0, 1.0+0i},
    {1.0+0i, 0}};
cx_mat Y = { {0, -1i},
            {1i, 0} };
cx_mat Z = {{1.0+0i, 0},
           {0, -1.0+0i}};
cx_mat I(2, 2, fill::eye);



cx_mat ZI = arma::kron(Z,I);

ZI.print();
  
cx_mat H_1= arma::zeros(2^(n+m),2^(n+m)); //Initialising Hamiltonians
cout << H_1 << endl;

}
/*mat X = {{0, 1},{1, 0}};
mat Y = {{0, eye},{eye, 0}};          
mat Z = {{1,0},{0,-1}};
mat I = {{1,0},{0,1}};
*/        

/*int H_1 = zeros(2^[n+m]); */ /*Initialising Hamiltonian*/

/*Intialise S for the if and else statements*/

// int H = zeros({2^[n+m]}{2^[n+m]}); /*Initialising Hamiltonian*/
//   for (int i = 0; i < m-1; i++){
//     /*Construct H_C */
//     /*order qubits: 0....n-1 and n....n+m-1 */
    
//     H_m = - (X + Z + I); 
//     if sign(inst[i][0] == -1){
//       S[0] = X;
//    }
//   else{
//         S[0] = Z;
//       }
//   }
    
//   for (int i = 0; i < m-1; i++){
//     int H_m = - (X + Z + I); 
//     if sign(inst[i][1] == -1){
//       S[1] = X;
//    }
//       else{
//         S[1] = Z;
//     }
//   }
// for (int i = 0; i < m-1; i++){
//     H_m = - (X + Z + I); 
//     if sign(inst[i][2] == -1){
//       S[2] = X;
//    }
//       else{
//         S[2] = Z;
//     }
//   }
// int qubits = abs(inst[1][:])-1;
//   if (qubits = 0){
//     H_k = S[]
//   }
//   else{
//     H_k = I; /* */
//   }
//   for (int j = 1; j < n-1; j++){
//     if(j=qubits){
//       H_k = kron(H_k,S)  
//     }
//     else{
//       H_k = kron(H_k, I)
//     }
//   }
// for (int j = n; j < n-m+1; j++){
//   if(j-n == i){
//     H_k = kron(H_k, H_m)
//   }
//   else{
//     H_k = kron(H_k, I)
//   }
// }
// H = H + H_k; 
// /*Intialise S for the if and else statements*/

// /* Compute Avg. Sign <S>*/

// int eH = expmat(-H);

// int avgSign = (trace(eH)/trace(abs(eH)));

// /* Compute Avg. Sign <S>*/


// return 0;
// }
