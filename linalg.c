#include "sieve.h"
#include <unistd.h>
#define SIZE_T_BITS (sizeof(size_t) * 8)

Matrix initMatrix(size_t rows, size_t columns)
{
  Matrix Matrix;
  Matrix.numRows = rows;
  Matrix.numCols = columns;
  Matrix.mem = calloc(rows * ((columns - 1) / SIZE_T_BITS + 1), sizeof(size_t));
  return Matrix;
}

void freeMatrix(Matrix Matrix)
{
  free(Matrix.mem);
}

/**
 * Given the sieved values and factorization table, returns a Matrix with the appropriate values
 */
Matrix loadMatrix(BigNumList currValues, FactorizationTable table)
{
  // determine how many values sieved down to 1
  int numRelations, i, j;
  for (i = 0; i < currValues.len; i++)
  {
    numRelations += (mpz_get_ui(currValues.nums[i]) == 1);
  }

  Matrix Matrix = initMatrix(table.numPrimes, numRelations);
  j = 0;
  for (i = 0; i < currValues.len; i++)
  {
    if (mpz_get_ui(currValues.nums[i]) == 1)
    {
      // update the ith column
      j++;
    }
  }
  return Matrix;
}

void set1(Matrix Matrix, size_t row, size_t col)
{
  size_t rowSize = (Matrix.numCols - 1) / SIZE_T_BITS + 1; // number of size_t's in a row
  size_t *data = Matrix.mem + row * rowSize + col / SIZE_T_BITS;
  //printf("%p\n", (void *)data);
  *data |= ((size_t)1) << (col % SIZE_T_BITS);
}

int get(Matrix Matrix, size_t row, size_t col)
{
  size_t rowSize = (Matrix.numCols - 1) / SIZE_T_BITS + 1; // number of size_t's in a row
  size_t *data = Matrix.mem + row * rowSize + col / SIZE_T_BITS;
  //printf("%lu / %lu = %lu\n", col, SIZE_T_BITS, col / SIZE_T_BITS);
  //printf("data = %zx\n", *data);
  return (*data & (((size_t)1) << (col % SIZE_T_BITS))) > 0;
}

void printMatrix(Matrix m)
{
  int i,j;
  for (i = 0; i < m.numRows; i++)
  {
    for (j = 0; j < m.numCols; j++)
    {
      printf("%d", get(m, i, j));
    }
    printf("\n");
  }
}

/**
 * Sets row2 to row1 + row2 (mod 2).
 */
void addRows(Matrix m, size_t row1, size_t row2)
{
  if (row1 == row2) return;
  size_t rowSize = (m.numCols - 1) / SIZE_T_BITS + 1; // number of size_t's in a row
  size_t *row1ptr = m.mem + row1 * rowSize;
  size_t *row2ptr = m.mem + row2 * rowSize;
  size_t i;
  for (i = 0; i < rowSize; i++)
  {
    *(row2ptr + i) ^= *(row1ptr + i); // use XOR for efficient addition mod 2
  }
}

void swapRows(Matrix m, size_t row1, size_t row2)
{
  if (row1 == row2) return;
  size_t rowSize = (m.numCols - 1) / SIZE_T_BITS + 1; // number of size_t's in a row
  size_t *row1ptr = m.mem + row1 * rowSize;
  size_t *row2ptr = m.mem + row2 * rowSize;
  size_t i;
  for (i = 0; i < rowSize; i++)
  {
    *(row2ptr + i) ^= *(row1ptr + i);
    *(row1ptr + i) ^= *(row2ptr + i);
    *(row2ptr + i) ^= *(row1ptr + i);
  }
}

void rref(Matrix m)
{
  size_t lead, r, i;
  size_t rowCount = m.numRows;
  size_t colCount = m.numCols;
  
  lead = 0;
  for (r=0; r<rowCount; r++) {
      if (lead >= colCount)
          return;
      i = r;
      while (!get(m, i,lead)) {
          i++;
          if (i == rowCount) {
              i = r;
              lead++;
              if (lead == colCount)
                  return;
          }
      }
      swapRows(m, i, r );
      for (i=0; i<rowCount; i++) {
          if ( i != r  && get(m,i,lead)) {
              addRows(m,r, i) ;
          }
      }
      lead++;
  }
}



typedef struct {
    Matrix m;
    size_t numVecs;
} PreppedMatr;

// Returns a matrix which makes it easier to find the basis of the kernel
PreppedMatr prepareMatrix(Matrix m)
{

  rref(m);


  size_t n = m.numCols;
  Matrix prepped = initMatrix(n,n);

  size_t i, numVecs = 0;
  for ( i = 0; i < m.numRows; i++)
  {
    size_t j = 0;
    while (!get(m,i,j) && j < n)
    {
        j++;
    }
    if (j == n)
        continue;

    size_t k;
    for (k = 0; k < n; k++)
    {
        if (get(m,i,k))
        {
            set1(prepped, j, k);
            numVecs++;
        }
    }
  }

  PreppedMatr mtr;
  mtr.m = prepped;
  mtr.numVecs = numVecs;
  return mtr;
}

Matrix kernelBasis(Matrix m)
{
  PreppedMatr prepped = prepareMatrix(m);
  freeMatrix(m);
  Matrix newMat = prepped.m;
  size_t n = newMat.numCols;
  size_t numVecs = prepped.numVecs;

  Matrix basis = initMatrix(numVecs, n);
  size_t c,v = 0;
  for (c = 0; c < n; c++)
  {
    if (!get(newMat,c,c))
    {
        int i;
        for (i=0; i < n; i++)
        {
            size_t val = get(newMat,i,c) | i==c;
            if (val)
            {
                set1(basis,v, i);
            }
        }
        v++;
    }
  }
  freeMatrix(newMat);

  return basis;
}


/*
int main()
{
  Matrix m = initMatrix(5,100);
  printMatrix(m);
  set1(m, 0, 2);
  set1(m, 0, 3);
  set1(m, 0, 32);
  set1(m, 0, 70);
  set1(m, 2, 30);
  set1(m, 2, 70);
  printf("__________\n");
  printMatrix(m);
  addRows(m, 0, 1);
  printf("__________\n");
  printMatrix(m);
  addRows(m, 0, 2);
  printf("__________\n");
  printMatrix(m);
  freeMatrix(m);
}*/
