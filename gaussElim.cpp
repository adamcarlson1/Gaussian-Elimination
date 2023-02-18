#include <iostream>
#include <string>
#include "stdio.h"
#include "math.h"
#include "string.h"
#include <cstdlib>
using namespace std;
using std :: cout;
using std :: endl;
using std :: string;


// =========================================================================
// Small Number Filters
// =========================================================================

double filter(double val)
{
  double smallVal = 1.e-10;

  if ( fabs(val) < smallVal ) return 0.;
  return val;
}


bool isZero(double val)
{
  double smallVal = 1.e-10;
  if ( fabs(val) < smallVal ) return true;
  return false;
  
}

bool isOne(double val)
{
  double smallVal = 1.e-10;
  if ( fabs(val-1.0) < smallVal ) return true;
  return false;
  
}

// =========================================================================
// printSystem
// =========================================================================

void printSystem(string message, int m , int n , double *A , double *b)
{
  
  cout << endl;
  cout << "--------------------------------------" << endl;
  cout << message << endl;
  cout << "--------------------------------------" << endl;
  cout << endl;

  for ( int i = 0 ; i < m ; ++i )
    {
      for ( int j = 0 ; j < n ; ++j )
	{
	  double out = A[i*n + j];
	  printf("%15.5e ", filter(out));
	}
      printf("      |    %15.5e ", filter(b[i]) );
      cout << endl;
    }
}

// =========================================================================
// swapRows
// =========================================================================

void swapRows( int m , int n , double *A , double *b, int row1, int row2 )
{
  double tmp[n];

  for ( int j = 0 ; j < n ; ++j )
    {
      tmp [          j ] = A   [ row1*n+j ];
      A   [ row1*n + j ] = A   [ row2*n+j ];
      A   [ row2*n + j ] = tmp [        j ];

    }
  
  double tmp_b = b[row1];
  b[row1] = b[row2];
  b[row2] = tmp_b; 
}


// =========================================================================
// findPivot
// =========================================================================

int findPivot( int m , int n , double *A , double *b , int row , int col )
{

  for ( int colTry = col ; colTry < n ; ++colTry )
    {
      // Find largest pivot in this column, not counting rows above "row"

      int largest = row;
  
      for ( int i = row ; i < m ; ++i )
	  if ( fabs( A[i*n + colTry] ) > fabs( A[ largest*n + colTry ] ) ) largest = i;

      // Swap rows

      swapRows(m,n,A,b,row,largest);

      if ( !isZero( A[row*n + colTry] ) ) return colTry;

    }

  return -1;

}


// =========================================================================
// M A I N   P R O G R A M
// =========================================================================

int main(int argc, char *argv[])
{
  int m; // number of rows
  int n; // number of columns

  #include "setCase.h"

  // --------------------------
  // Save Original
  // --------------------------

  double Asave[m][n];
  double bsave[n]   ;
  int    pivotCol[m];

  for ( int i = 0 ; i < m ; ++i )
    {
    bsave[i] = b[i];
    for ( int j = 0 ; j < n ; ++j ) Asave[i][j] = A[i][j];
    }

  printSystem("Original Matrix" , m , n , &A[0][0] , &b[0] );

  // --------------------------
  // Forward Phase
  // --------------------------

  // Work from the top down, i.e., starting with pivotRow = 0,
  // eliminating under the diagonal for each row

  for ( int pivotRow = 0 ; pivotRow < m ; ++pivotRow )
    {
      pivotCol[pivotRow] = findPivot(m,n, &A[0][0], &b[0], pivotRow, pivotRow);

      for ( int elimRow = pivotRow + 1 ; elimRow < m ; ++ elimRow )
	{
	  double fac = -A[elimRow][pivotCol[pivotRow]]/A[pivotRow][pivotCol[pivotRow]];

	  for ( int col = pivotCol[pivotRow] ; col < n ; ++col ) A[elimRow][col] += A[pivotRow][col] * fac;

	  b[elimRow] += b[pivotRow] * fac;
	}

      printSystem("After elimination on row " + std::to_string(pivotRow) + " pivot Col = " + std::to_string(pivotCol[pivotRow]) , m , n , &A[0][0]  , &b[0] );

    }

  // --------------------------
  // Print Echelon Form
  // --------------------------

  printSystem("Echelon Form" , m , n , &A[0][0]  , &b[0]  );

  // --------------------------
   // Backward Phase
  // --------------------------

  for ( int pivotRow = m-1 ; pivotRow >= 0 ; --pivotRow )
    {
      if ( ! isZero(A[pivotRow][pivotCol[pivotRow]] ) )
      for ( int elimRow = pivotRow - 1 ; elimRow >= 0 ; -- elimRow )
	{

	  double fac = -A[elimRow][pivotCol[pivotRow]]/A[pivotRow][pivotCol[pivotRow]];

	  for ( int col = 0 ; col < n ; ++col ) A[elimRow][col] += A[pivotRow][col] * fac;

	  b[elimRow] += b[pivotRow] * fac;
	}
    }

  printSystem("Semi-reduced echelon: "  , m , n , &A[0][0]  , &b[0] );

  // Normalize pivots

  for ( int i = 0 ; i < m ; ++i )
    {
      double denom = A[i][pivotCol[i] ];
      if ( ! isZero(denom) )
      {
	b[i] /= denom;
	for ( int j = 0 ; j < n ; ++j ) A[i][j] /= denom;
      }
    }

  printSystem("Reduced echelon: " , m , n , &A[0][0]  , &b[0] );

  // --------------------------
  // Find free variables
  // --------------------------

  // Programming Assignment (part 1) --------------------- Begin
  //
  //
  //
  // Write lines of code that compute freeVars[].
  //
  //
  //
  // Programming Assignment (part 1) --------------------- Begin
  int freeVars[n]; // freeVars[i] = 1 if variable i is free, 0 if the variable i is not free
  for ( int i = 0 ; i < n ; ++i ) {
      freeVars[i] = 1; // mark all variables as free (1
  }
  for ( int i = 0 ; i < m ; ++i ){ // loop through the pivot columns
        freeVars[pivotCol[i]] = 0; // if the variable is a pivot column, mark it as not free (0)
  }
  // Programming Assignment (part 1) --------------------- End

    // --------------------------
    // Print free variables
    // --------------------------

  cout << "\n";
  cout << "\n";
  cout << "Free Variables: \n";
  cout << "\n";
  cout << "\n";


  for ( int i = 0 ; i < n ; ++i ) cout << "freeVars[" << i << "] = " << freeVars[i] << endl; // print the free variables (if any are found)

  // --------------------------
  // Print space
  // --------------------------

  // Programming Assignment (part 2) --------------------- Begin
  //
  //
  //
  // Write lines of code that print the general solution, which
  // might be unique or in terms of the free variables.
  //
  //
  //
  // Programming Assignment (part 2) --------------------- Begin

  cout << "\n";
  cout << "\n";
  cout << "General Solution: \n";
  cout << "\n";
  cout << "\n";

    int freeVarCount = 0; // number of free variables
    for (int i = 0; i < n; i++) { // loop through the free variables
        if(freeVars[i]){ // if the variable is free
            freeVarCount++; // increment the number of free variables
        }
    }
    cout << "General Solution in terms of " << freeVarCount << " free variables: \n";
    cout << "\n";
    cout << "\n";
    // print the general solution in terms of the free variables (if there are any)
    for (int x = 0; x < m; x++) { // loop through the all the rows of A (the pivot columns)
        if (!freeVars[x]) { // if the variable is not free, then print the solution for each pivot column
            cout << "x[" << x << "]" << " = "; // print the variable name
            for (int y = 0; y < n; y++) { // loop through the columns of A (which are the free variables)
                if (freeVars[y]) { // if the variable is free
                    cout << " " << "( " << 0 - A[x][y] << " )" << " + x[" << y // output the solution in terms of the free variables
                        << "]" << " * ";
                }
            }
            cout << "( " << b[x] << " )" << endl; // for each pivot column, print the solution for that column
        }
    }
}
// Programming Assignment (part 2) --------------------- End
