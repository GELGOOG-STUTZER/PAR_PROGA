#include<mpi.h>
#include<stdio.h>

int main(int argc, char *argv[]) {
	int commsize, my_rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &commsize);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	long double sum=0, sum0;
	int N;
	MPI_Status status;
	if(my_rank == 0) {
        	scanf("%d", &N);
		int send_status;
		for(int to_thread=1; to_thread<commsize; to_thread++) {
			MPI_Send(&N, 1, MPI_INT, to_thread, 0, MPI_COMM_WORLD);
		}
		for(int to_thread=1; to_thread<commsize; to_thread++) {
			MPI_Recv(&sum0, 1, MPI_LONG_DOUBLE, to_thread, 0, MPI_COMM_WORLD, &status);
		sum = sum + sum0;
		}
		printf("SUMMA = %Lf\n", sum);
		
	}

	if(my_rank != 0) {
		MPI_Recv(&N, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		int i;
		for(i=0; (i+1)*(commsize-1)<=N; i++) {
			sum = sum + 1.0/(i*(commsize-1) + my_rank);
		}
		if(N%(commsize-1) >= my_rank) {
			sum = sum + 1.0/(N - (N%(commsize-1)) + my_rank);
		}
		MPI_Send(&sum, 1, MPI_LONG_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}


	MPI_Finalize();
}
