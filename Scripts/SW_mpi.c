#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdlib.h>
#include <mpi.h>

#define MAX_SEQ 50

#define CHECK_NULL(_check) {\
   if ((_check)==NULL) \
      fprintf(stderr, "Null Pointer allocating memory\n");\
   }


#define AA 20           // number of amino acids
#define MAX2(x,y)     ((x)<(y) ? (y) : (x))
#define MAX3(x,y,z)   (MAX2(x,y)<(z) ? (z) : MAX2(x,y))

// function prototypes
void error(char *);   /** error handling */

int char2AAmem[256];
int AA2charmem[AA];
void initChar2AATranslation(void);

int getRowCount(int rowsTotal, int mpiRank, int mpiSize) {
    /* Adjust slack of rows in case rowsTotal is not exactly divisible */
    return (rowsTotal / mpiSize) + (rowsTotal % mpiSize > mpiRank);
}

int getColCount(int colsTotal, int index, int numBlocks) {
    /* Adjust slack of rows in case rowsTotal is not exactly divisible */
    return (colsTotal / numBlocks) + (colsTotal % numBlocks > index);
}

int main(int argc, char *argv[]) {

  /* DECLARE VARIABLES */
  FILE * in1, *in2, *pam;
  char ch;
  int temp;
  int i,j,diag,down,right,DELTA;
  char *aout,*bout;
  int max, Max;

  short *a, *b;
  int sim[AA][AA];		// PAM similarity matrix
  int N;
  int nc;

  int rank, size;
  int *hptrLocal;
  int **hLocal;
  int rowsLocal;
  int numBlocks;


  /*Error handling for input file */
  if ((( argc != 5) && (argc!=6)) && (argc!=7)) {
    fprintf(stderr,"%s protein1 protein2 PAM gapPenalty [N]\n",argv[0]);
    exit(1);
  } else if (argc==5)  /* Maximum size of the proteins, they should have equal size */
  {
    /***** Initialization of input file and pair array **/
    in1   = fopen(argv[1],"r");
    in2   = fopen(argv[2],"r");
    pam   = fopen(argv[3],"r");
    DELTA = atoi(argv[4]);
    N = MAX_SEQ;
  } else if (argc==6){
    in1   = fopen(argv[1],"r");
    in2   = fopen(argv[2],"r");
    pam   = fopen(argv[3],"r");
    DELTA = atoi(argv[4]);
    N     = atoi(argv[5]);
  } else {
    in1   = fopen(argv[1],"r");
    in2   = fopen(argv[2],"r");
    pam   = fopen(argv[3],"r");
    DELTA = atoi(argv[4]);
    N     = atoi(argv[5]);
    numBlocks = atoi(argv[6]);
  }

  // Star MPI world
  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );



  /* How many rows has each process? */

  int rows_rank[size];
  for (int i=0; i<size; i++) rows_rank[i] = getRowCount(N+1, i, size);
  rowsLocal = rows_rank[rank];
  if (rank!=0) rowsLocal++;

  MPI_Barrier(MPI_COMM_WORLD);
  // Alocate the memory
  // Test the allocation is correct

  CHECK_NULL((aout = (char *) malloc(sizeof(char)*2*N)));
  CHECK_NULL((bout = (char *) malloc(sizeof(char)*2*N)));
  CHECK_NULL((a = (short *) malloc(sizeof(short)*(N+1))));
  CHECK_NULL((b = (short *) malloc(sizeof(short)*(N+1))));
  CHECK_NULL((hptrLocal = (int *) malloc(sizeof(int)*(rowsLocal)*(N+1))));
  CHECK_NULL((hLocal = (int **) malloc(sizeof(int*)*(rowsLocal))));
  /* Mount h[N][N] */

  for(i=0;i<rowsLocal;i++) hLocal[i]=hptrLocal+i*(N+1);

  for (i=0;i<rowsLocal;i++) hLocal[i][0]=0;
  for (j=0;j<=N;j++) hLocal[0][j]=0;

  /* Reading the files */

  initChar2AATranslation();
  /* end  AMPP */

  /** read PAM250 similarity matrix **/
  /* begin AMPP */
  fscanf(pam,"%*s");
  /* end  AMPP */
  for (i=0;i<AA;i++)
    for (j=0;j<=i;j++) {
        if (fscanf(pam, "%d ", &temp) == EOF) {
          printf(stderr, "PAM file empty\n");
          fclose(pam);
          exit(1);
        }
        sim[i][j]=temp;
      }
  fclose(pam);
  for (i=0;i<AA;i++)
    for (j=i+1;j<AA;j++)
      sim[i][j]=sim[j][i]; // symmetrify

  /* begin AMPP */
  /** read first file in array "a" **/
  i=0;
  do {
    nc=fscanf(in1,"%c",&ch);
    if (nc>0 && char2AAmem[ch]>=0) {
      a[++i] = char2AAmem[ch];
    }
  } while (nc>0 && (i<N));
  a[0]=i;
  fclose(in1);

  /** read second file in array "b" **/
  i=0;
  do {
    nc=fscanf(in2,"%c",&ch);
    if (nc>0 && char2AAmem[ch]>=0) {
      b[++i] = char2AAmem[ch];
    }
  } while (nc>0 && (i<N));
  b[0]=i;
  fclose(in2);
  /* Files done! */

  // Set uo de blocks Sizes
  int col_block[numBlocks];
  for (int bl=0; bl<numBlocks; bl++) col_block[bl] = getColCount(N, bl, numBlocks);

  /** Tha algorithm **/
  Max=0;
  MPI_Status status;
  int iA;
  int prevRows = -1;
  for (int i=0; i<rank; i++) prevRows = prevRows + rows_rank[i];
  if ( rank == 0) prevRows = 0;
  for (int bl=0; bl<numBlocks; bl++){

      int sum = 1;
      for (int i=0; i<bl; i++) sum = sum + col_block[i];

      if (rank != 0) MPI_Recv(&hLocal[0][sum], col_block[bl], MPI_INT, rank-1, 1, MPI_COMM_WORLD, &status);
      for (i=1; i<rowsLocal; i++)
        for (j=sum;j<(col_block[bl]+sum);j++) {
          iA = i + prevRows;
          diag    = hLocal[i-1][j-1] + sim[a[iA]][b[j]];
          down    = hLocal[i-1][j] + DELTA;
          right   = hLocal[i][j-1] + DELTA;
          max=MAX3(diag,down,right);
          if (max <= 0)         hLocal[i][j]=0;
          else if (max == diag) hLocal[i][j]=diag;
          else if (max == down) hLocal[i][j]=down;
          else                  hLocal[i][j]=right;
          if (max > Max) Max=max;
        }
      if (rank != (size-1)) MPI_Send(&hLocal[rowsLocal-1][sum], col_block[bl], MPI_INT, rank+1, 1, MPI_COMM_WORLD);
  }
  int score;
  MPI_Reduce(&Max, &score, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
  if (rank==0) printf("The score of the aligment is %d.\n", score);

  MPI_Finalize();
}





void error(char * s) {
  fprintf(stderr,"%s\n",s);
  exit(1);
}

void initChar2AATranslation(void) {
    int i;
    for(i=0; i<256; i++) char2AAmem[i]=-1;
    char2AAmem['c']=char2AAmem['C']=0;
    AA2charmem[0]='c';
    char2AAmem['g']=char2AAmem['G']=1;
    AA2charmem[1]='g';
    char2AAmem['p']=char2AAmem['P']=2;
    AA2charmem[2]='p';
    char2AAmem['s']=char2AAmem['S']=3;
    AA2charmem[3]='s';
    char2AAmem['a']=char2AAmem['A']=4;
    AA2charmem[4]='a';
    char2AAmem['t']=char2AAmem['T']=5;
    AA2charmem[5]='t';
    char2AAmem['d']=char2AAmem['D']=6;
    AA2charmem[6]='d';
    char2AAmem['e']=char2AAmem['E']=7;
    AA2charmem[7]='e';
    char2AAmem['n']=char2AAmem['N']=8;
    AA2charmem[8]='n';
    char2AAmem['q']=char2AAmem['Q']=9;
    AA2charmem[9]='q';
    char2AAmem['h']=char2AAmem['H']=10;
    AA2charmem[10]='h';
    char2AAmem['k']=char2AAmem['K']=11;
    AA2charmem[11]='k';
    char2AAmem['r']=char2AAmem['R']=12;
    AA2charmem[12]='r';
    char2AAmem['v']=char2AAmem['V']=13;
    AA2charmem[13]='v';
    char2AAmem['m']=char2AAmem['M']=14;
    AA2charmem[14]='m';
    char2AAmem['i']=char2AAmem['I']=15;
    AA2charmem[15]='i';
    char2AAmem['l']=char2AAmem['L']=16;
    AA2charmem[16]='l';
    char2AAmem['f']=char2AAmem['F']=17;
    AA2charmem[17]='L';
    char2AAmem['y']=char2AAmem['Y']=18;
    AA2charmem[18]='y';
    char2AAmem['w']=char2AAmem['W']=19;
    AA2charmem[19]='w';
}
