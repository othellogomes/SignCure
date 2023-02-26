#include <iostream>
#include <armadillo> 

using namespace std;
using namespace arma;

/*Generate 3-SAT instances in this part of the code*/

int main()
{
m = 2; 
n = 4;  
inst = [m][3];
for (int i = 0,1; i < m-1; i++) {
  for (int j = 0,1,2; j++){
  var = rand() % n + 1;
  sign = rand() % 4 - 4;
  inst [i][j] = sign * var;
}
}






/*Generate Hamiltonians in this part of the code*/

/*Pauli Matrices: need to make sure they're dense matrices*/

mat X = {{0, 1}{1, 0}};
mat Y = {{0, -i}{i, 0}};          
mat Z = {{1,0}{0,-1}};
mat I = {{1,0}{0,1}}         

Hamiltonian = zeros({2^{n+m}}{2^{n+m}}); /*Initialising Hamiltonian*/



/* Compute Avg. Sign <S>*/









}
