lcg_rand <- function(int N, double SEED, double* DUNIF, int NDIM)
{
        int    i;
        double  DZ, DOVER, DZ1, DZ2, DOVER1, DOVER2;
        double  DTWO31, DMDLS, DA1, DA2;

        DTWO31 = 2147483648.0; /* DTWO31=2**31  */
        DMDLS  = 2147483647.0; /* DMDLS=2**31-1 */
        DA1 = 41160.0;       /* DA1=950706376 MOD 2**16 */
        DA2 = 950665216.0;   /* DA2=950706376-DA1 */

        DZ = SEED;
        if ( N > NDIM )
                N = NDIM;
        for ( i=1; i<=N; i++ ) {
                DZ = floor(DZ);
                DZ1 = DZ*DA1;
                DZ2 = DZ*DA2;
                DOVER1 = floor(DZ1/DTWO31);
                DOVER2 = floor(DZ2/DTWO31);
                DZ1 = DZ1-DOVER1*DTWO31;
                DZ2 = DZ2-DOVER2*DTWO31;
                DZ = DZ1+DZ2+DOVER1+DOVER2;
                DOVER = floor(DZ/DMDLS);
                DZ = DZ-DOVER*DMDLS;
                DUNIF[i-1] = DZ/DMDLS;
                SEED = DZ;
        }

        return SEED;
}

lcg <- function()
{
        double  *DUNIF, SEED;
        int             i, counter;
        unsigned bit;
        int             num_0s, num_1s, v, bitsRead;

        SEED = 23482349.0;
        if ( ((epsilon = (BitSequence *) calloc(tp.n, sizeof(BitSequence))) == NULL) ||
             ((DUNIF = (double*)calloc(tp.n, sizeof(double))) == NULL) ) {
                printf("Insufficient memory available.\n");
                exit(1);
        }
        counter = 1;
 
        for ( v=0; v<tp.numOfBitStreams; v++ ) {
                num_0s = 0;
                num_1s = 0;
                bitsRead = 0;
                SEED = lcg_rand(tp.n, SEED, DUNIF, tp.n);
                for ( i=0; i<tp.n; i++ ) {
                        if ( DUNIF[i] < 0.5 ) {
                                bit = 0;
                                num_0s++;
                        }
                        else {
                                bit = 1;
                                num_1s++;
                        }
                        bitsRead++;
                        epsilon[i] = bit;
                }
                fprintf(freqfp, "\t\tBITSREAD = %d 0s = %d 1s = %d\n", bitsRead, num_0s, num_1s); fflush(freqfp);
                nist_test_suite();
                }
        free(DUNIF);
        free(epsilon);
}