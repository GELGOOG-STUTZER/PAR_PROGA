#include<mpi.h>
#include<stdio.h>

int main(int argc, char *argv[]) {
	int commsize, my_rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &commsize);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	int N;
	MPI_Status status;
	if(my_rank == 0) {
        	scanf("%d", &N);
		int send_status;
		printf("RANK = %d   N = %d\n", my_rank, N);
		MPI_Send(&N, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
		MPI_Recv(&N, 1, MPI_INT, commsize-1, 0, MPI_COMM_WORLD, &status);
		printf("RESULT = %d\n", N);
		
	}

	if(my_rank != 0) {
		MPI_Recv(&N, 1, MPI_INT, my_rank-1, 0, MPI_COMM_WORLD, &status);
		N = N + 2;
		printf("RANK = %d   N = %d\n", my_rank, N);
		MPI_Send(&N, 1, MPI_INT, (my_rank+1)%commsize, 0, MPI_COMM_WORLD);
	}


	MPI_Finalize();
}
