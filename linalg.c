#include "sieve.h"
#define SIZE_T_BITS (sizeof(size_t) * 8)

Matrix initMatrix(size_t rows, size_t columns)
{
  Matrix matrix;
  matrix.numRows = rows;
  matrix.numCols = columns;
  matrix.mem = calloc(rows * ((columns - 1) / SIZE_T_BITS + 1), sizeof(size_t));
  return matrix;
}

void freeMatrix(Matrix matrix)
{
  free(matrix.mem);
}

/**
 * Given the sieved values and factorization table, returns a matrix with the appropriate values
 */
Matrix loadMatrix(BigNumList currValues, FactorizationTable table)
{
  // determine how many values sieved down to 1
  int numRelations, i, j;
  for (i = 0; i < currValues.len; i++)
  {
    numRelations += (mpz_get_ui(currValues.nums[i]) == 1);
  }

  Matrix matrix = initMatrix(table.numPrimes, numRelations);
  j = 0;
  for (i = 0; i < currValues.len; i++)
  {
    if (mpz_get_ui(currValues.nums[i]) == 1)
    {
      // update the ith column
      j++;
    }
  }
  return matrix;
}

void set1(Matrix matrix, size_t row, size_t col)
{
  size_t rowSize = (matrix.numCols - 1) / SIZE_T_BITS + 1; // number of size_t's in a row
  size_t *data = matrix.mem + row * rowSize + col / SIZE_T_BITS;
  //printf("%p\n", (void *)data);
  *data |= ((size_t)1) << (col % SIZE_T_BITS);
}

int get(Matrix matrix, size_t row, size_t col)
{
  size_t rowSize = (matrix.numCols - 1) / SIZE_T_BITS + 1; // number of size_t's in a row
  size_t *data = matrix.mem + row * rowSize + col / SIZE_T_BITS;
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
// TODO unfinished, Yoni idk if you wanna use this in the null space thing
void rref(Matrix m)
{
  size_t i, j, leadingEntryIndex;
  leadingEntryIndex = 0; // k will be index of leading entry of current row
  for (i = 0; i < m.numCols; i++)
  {
    if (get(m, i, leadingEntryIndex))
    {

    } else
    {
      for (j = i + 1; j < m.numRows; j++)
      {
        if (get(m, j, leadingEntryIndex))
        {

        }
      }
    }

  }
}

/**
 * Sets row2 to row1 + row2 (mod 2).
 */
void addRows(Matrix m, size_t row1, size_t row2)
{
  size_t rowSize = (m.numCols - 1) / SIZE_T_BITS + 1; // number of size_t's in a row
  size_t *row1ptr = m.mem + row1 * rowSize;
  size_t *row2ptr = m.mem + row2 * rowSize;
  size_t i;
  for (i = 0; i < rowSize; i++)
  {
    *(row2ptr + i) ^= *(row1ptr + i); // use XOR for efficient addition mod 2
  }
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
