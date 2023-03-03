/* Othello D. Gomes
Sign Curing Numerics
*/
#include <iostream>
#include <armadillo> 
#include <complex>  
using namespace std;
using namespace arma;

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
std::complex<double> imaginary(1);

cx_mat X = {{0, 1},
          {1, 0}};
cx_mat Y = { {0, -imaginary},
          {imaginary, 0} };
cx_mat Z = {{1, 0},
          {0, -1}};
cx_mat I(2, 2, fill::eye);

cout << X << endl;

/*mat X = {{0, 1},{1, 0}};
mat Y = {{0, eye},{eye, 0}};          
mat Z = {{1,0},{0,-1}};
mat I = {{1,0},{0,1}};
*/        

/*int H_1 = zeros(2^[n+m]); */ /*Initialising Hamiltonian*/

/*Intialise S for the if and else statements*/

/* Compute Avg. Sign <S>*/

return 0;
}
