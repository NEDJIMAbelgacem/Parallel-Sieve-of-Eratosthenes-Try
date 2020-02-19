// Sieve-of-Eratosthenes-parallel.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <mpi.h>
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <array>
#include <random> 
#include <algorithm>

// LIMIT_N should be multiple of world_size
#define LIMIT_N 700000l
#define NB_NODES 8
#define ROOT_RANK 0
#define NUMBERS_BUFF_TAG 1
#define CALCULATION_STATE_TAG 2
#define SEED 0
#define NB_NODE_PRIMES 10

#define root_get_marker_v(n) (marker_array[(n % world_size) * (LIMIT_N / world_size) + n / world_size])

// n should be in process data range which means n % world_size == process_rank
#define node_get_marker_v(n) (marker_array[n / world_size])

int* marker_array;

int world_size, process_rank;
std::vector<int> global_primes;

int prime_search_index;

void root_initialize_marker() {
	marker_array = new int[LIMIT_N];
	for (int i = 0; i < LIMIT_N / world_size; ++i) marker_array[i] = false;
	// world_size should be prime so marker_array[1] which correspond to 7 * 1 is prime
	marker_array[1] = true;

	for (int i = LIMIT_N / world_size; i < LIMIT_N; ++i) marker_array[i] = true;
}

void node_initialize_marker() {
	marker_array = new int[LIMIT_N / world_size];
	for (int i = 0; i < LIMIT_N / world_size; ++i) marker_array[i] = true;
	if (process_rank == 1) marker_array[0] = false;
}

bool root_send_calculation_end_notification(int ended, int dest_node) {
	MPI_Send(&ended, 1, MPI_INT, dest_node, CALCULATION_STATE_TAG, MPI_COMM_WORLD);
}

bool node_receive_calculation_end_notification() {
	MPI_Status status;
	int calculation_ended;
	MPI_Recv(&calculation_ended, 1, MPI_INT, ROOT_RANK, CALCULATION_STATE_TAG, MPI_COMM_WORLD, &status);
	return calculation_ended;
}

bool is_divisible_by(int n, std::vector<int>& primes) {
	for (int p : primes) {
		if (n % p == 0) return true;
	}
	return false;
}

int last_prime = 2;
std::vector<int> to_be_sent_primes;
void root_get_next_primes(int nb_primes) {
	to_be_sent_primes.clear();
	while (last_prime < LIMIT_N && to_be_sent_primes.size() < nb_primes) {
		if (root_get_marker_v(last_prime) && !is_divisible_by(last_prime, to_be_sent_primes)) {
			to_be_sent_primes.push_back(last_prime);
			global_primes.push_back(last_prime);
		}
		last_prime++;
	}
}

void root_send_primes() {
	int* primes_to_send = new int[to_be_sent_primes.size()];
	for (int i = 0; i < to_be_sent_primes.size(); ++i) primes_to_send[i] = to_be_sent_primes[i];
	for (int node = 1; node < world_size;  ++node) {
		MPI_Send(primes_to_send, to_be_sent_primes.size(), MPI_INT, node, NUMBERS_BUFF_TAG, MPI_COMM_WORLD);
	}
	delete[] primes_to_send;
}

int* node_received_primes;
int node_received_primes_count;
void node_receive_primes() {
	MPI_Status status;
	MPI_Recv(node_received_primes, world_size, MPI_INT, ROOT_RANK, NUMBERS_BUFF_TAG, MPI_COMM_WORLD, &status);
	MPI_Get_count(&status, MPI_INT, &node_received_primes_count);
}

// apply sieve algorithm on process data using numbers in node_received_primes
void node_apply_sieve() {
	for (int i = 0; i < node_received_primes_count; ++i) {
		int& p = node_received_primes[i];
		for (int n = 2 * p; n < LIMIT_N; n += p) {
			if (n % world_size != process_rank) continue;
			node_get_marker_v(n) = false;
		}
	}
}

void root_receive_data() {
	MPI_Status status;
	for (int i = 1; i < world_size; ++i) {
		MPI_Recv(marker_array + i * LIMIT_N / world_size, LIMIT_N / world_size, MPI_INT, i, NUMBERS_BUFF_TAG, MPI_COMM_WORLD, &status);
		//int count;
		//MPI_Get_count(&status, MPI_INT, &count);
	}
}

int node_primes[NB_NODE_PRIMES + 1];
void node_find_primes() {
	int& nb_primes = node_primes[0];

}

void node_initialize_data() {
	node_primes[0] = 1;
	node_primes[1] = process_rank * LIMIT_N / world_size;
}

void node_send_data() {
	MPI_Send(marker_array, LIMIT_N / world_size, MPI_INT, ROOT_RANK, NUMBERS_BUFF_TAG, MPI_COMM_WORLD);
}

int main(int argc, char* argv[])
{
	double t1, t2, broadcast_duration = 0.0;
	MPI_Init(&argc, &argv);
	MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	// Get the rank of the process
	MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);

	node_received_primes = new int[world_size];
	if (process_rank == ROOT_RANK) {
		double t1, t2;
		double data_transmission_delay = 0.0;
		t1 = MPI_Wtime();
		root_initialize_marker();
		while (last_prime < LIMIT_N) {
			for (int node = 1; node < world_size; ++node) 
				root_send_calculation_end_notification(false, node);
			root_get_next_primes(world_size - 1);
			root_send_primes();
			root_receive_data();
		}
		for (int node = 1; node < world_size; ++node)
			root_send_calculation_end_notification(true, node);
		t2 = MPI_Wtime();
		std::cout << "parallel computation time : " << (t2 - t1) << std::endl;
		for (int i = 0; i < 100; ++i) std::cout << global_primes[i] << std::endl;
	}
	else {
		node_initialize_marker();
		double data_receive_time = 0.0;
		double data_send_time = 0.0;
		while (!node_receive_calculation_end_notification()) {
			data_receive_time -= MPI_Wtime();
			node_receive_primes();
			data_receive_time += MPI_Wtime();
			//std::cout << "NODE" << process_rank << " : received primes" << std::endl;
			node_apply_sieve();
			data_send_time -= MPI_Wtime();
			node_send_data();
			data_send_time += MPI_Wtime();
		}
		std::cout << "NODE" << process_rank << " data receive time " << data_receive_time << std::endl;
		std::cout << "NODE" << process_rank << " data send time " << data_send_time << std::endl;
	}
	MPI_Finalize();
}
