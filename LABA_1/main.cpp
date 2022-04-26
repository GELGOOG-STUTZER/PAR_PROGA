#include <mpi.h>
#include <iostream>
#include <vector>
#include <array>
#include <cassert>
#include <fstream>
#include <math.h>


constexpr size_t X_max = 100;
constexpr size_t T_max = 100;
constexpr double dt = 0.01;
constexpr double dx = 0.01;
constexpr size_t X_arr_size = X_max / dx;
constexpr size_t T_arr_size = T_max / dt;

MPI_Status status;

struct Start_Cond {
	std::array<double, X_arr_size> x0;
	std::array<double, T_arr_size> t0;
};

Start_Cond get_start_cond() {
	Start_Cond I;

	for(size_t i = 0; i != X_arr_size; ++i) {
		I.x0[i] = sin(static_cast<double>(i) * 2 * dx);
	}

	for(size_t i = 0; i != T_arr_size; ++i) {
		I.t0[i] = sin(static_cast<double>(i) * dt);
	}

	return I;
}

double f(double x, double t) {
	return cos(x) + cos(t);
}

double single_thread() {
	double time0 = MPI_Wtime();
	std::array<std::array<double, X_arr_size>, T_arr_size> data;
	Start_Cond I = get_start_cond();

	for(size_t i = 0; i != X_arr_size; ++i) {
		data[0][i] = I.x0[i];
	}

	for(size_t i = 1; i != T_arr_size; ++i) {
		data[i][0] = I.t0[i];
	}

	const double k0 = 2 * dt;
	const double k1 = 2 * dx;
	const double k2 = (dt + dx) / (2 * dx * dt);
	double time1 = MPI_Wtime();

	for(int t = 1; t != T_arr_size; ++t) {
		for(int x = 1; x != X_arr_size; ++x) {
			data[t][x] = (f((x + 0.5) * dx, (t - 0.5) * dt) -
				     (data[t][x - 1] - data[t - 1][x - 1] - data[t - 1][x]) / k0 -
				     (data[t - 1][x] - data[t - 1][x - 1] - data[t][x - 1]) / k1) / k2;
		}
	}

	double time2 = MPI_Wtime();
	std::cout << "OVERALL RUN TIME = " << time2 - time0 << std::endl;
	std::cout << "COMPUTING TIME = " << time2 - time1 << std::endl;
	std::ofstream out;
	out.open("Output.txt");

	for(size_t i = 0; i != X_arr_size; ++i) {

		for(size_t j = 0; j != T_arr_size; ++j) {
			out << data[i][j] << " ";
		}

		out << std::endl;

	}

	return data[T_arr_size - 1][X_arr_size - 1];
}

void multi_thread() {
	double time0 = MPI_Wtime();
	int t_num, rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &t_num);

//	if(t_num == 1) {
//		double res;
//		res = single_thread();
//		return;
//	}
	int thread_size = (X_arr_size - 1) / t_num;
	int offset = (X_arr_size - 1) % t_num;
	
	if(rank < offset) {
		thread_size++;
	}

	int x_start;

	if(rank < offset) {
		x_start = rank * thread_size + 1;
	}
	else {
		x_start = rank * thread_size + offset + 1;
	}

	std::vector<std::vector<double>> data(T_arr_size);

	for(int i = 0; i != T_arr_size; ++i) {
		data[i] = std::vector<double>(thread_size);
	}

	const double k0 = 2 * dt;
        const double k1 = 2 * dx;
        const double k2 = (dt + dx) / (2 * dx * dt);
	const Start_Cond I = get_start_cond();

	for(int i = 0; i != thread_size; ++i) {
		data[0][i] = I.t0[x_start + i];
	}

	double time1 = MPI_Wtime();

	for(int t = 1; t != T_arr_size; ++t) {
		double left_prev, left_cur;
		
		if(rank == 0) {
			left_prev = I.x0[t - 1];
			left_cur = I.x0[t];
		}
		else {
			if(t == 1) {
				left_prev = I.t0[x_start - 1];
			}
			else {
				left_prev = left_cur;
			}
			assert(MPI_Recv(&left_cur, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &status) == MPI_SUCCESS);
		}

		double lcur = left_cur, lpr = left_prev;

		for(int n = 0; n != thread_size; ++n) {
			data[t][n] = (f((n+x_start +0.5) * dx, (t - 0.5) * dt) -
				     (lcur - lpr - data[t - 1][n]) / k0 -
				     (data[t - 1][n] - lcur - lpr) / k1) / k2;
			lpr = data[t - 1][n];
			lcur = data[t][n];
		}

		if(rank != t_num - 1) {
			assert(MPI_Send(&(data[t][thread_size - 1]), 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD) == MPI_SUCCESS);
		}

	}

	if(rank == t_num - 1) {
		double time2 = MPI_Wtime();
		std::cout << "OVERALL RUN TIME = " << time2 - time0 << std::endl;
        	std::cout << "COMPUTING TIME = " << time2 - time1 << std::endl;
		std::cout << data[T_arr_size - 1][thread_size - 1] << std::endl;
	}
}

int main(int argc, char *argv[]) {
	int rank;
	MPI_Init(&argc, &argv);
	multi_thread();
	MPI_Finalize();
	return 0;
}

//10000 X 10000
//n = 1 TIME = 7.24s
//n = 2 TIME = 3.76s
//n = 3 TIME = 2.51s
//n = 4 TIME = 1.95s
//
//n = 2
//S = 7.24 / 3.76 = 1.93
//E = 1.93 / 2    = 0.96
//
//n = 3
//S = 7.24 / 2.51 = 2.88
//E = 2.88 / 3    = 0.96
//
//n = 4
//S = 7.24 / 1.95 = 3.71
//E = 3.71 / 4    = 0.93
