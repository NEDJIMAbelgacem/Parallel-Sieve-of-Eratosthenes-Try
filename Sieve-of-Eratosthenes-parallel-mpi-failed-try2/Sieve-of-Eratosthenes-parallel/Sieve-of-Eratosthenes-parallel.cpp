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

#define LIMIT_N 10000000l
#define ROOT_RANK 0
#define NUMBERS_BUFF_TAG 1
#define CALCULATION_STATE_TAG 2
#define SEED 0

bool* marker_array;

int world_size, process_rank;

int* numbers = nullptr;
int numbers_count = 0;

double node_to_root_communication_delays = 0.0;
double root_to_nodes_communication_delays = 0.0;
double local_sieve_computation_cost = 0.0;

void run_sieve() {
	/*
	std::map<int, bool> is_prime;
	//std::set<int> ordered_primes;
	for (int i = 0; i < numbers_count; ++i) {
		is_prime[numbers[i]] = true;
		//ordered_primes.insert(numbers[i]);
	}
	sieve_result.clear();
	for (int j = 0; j < numbers_count; ++j) {
		int& p = numbers[j];
		if (is_prime[p]) {
			sieve_result.push_back(p);
			for (int i = 2 * p; i < LIMIT_N; i += p) {
				if (is_prime.find(i) != is_prime.end()) is_prime[i] = false;
			}
		}
	}*/
	
	//std::vector<int> v;
	//sieve_result.clear();
	int k = 0;
	for (int i = 0; i < numbers_count; ++i) {
		if (!marker_array[numbers[i]]) continue;
		for (int j = numbers[i] * 2; j < LIMIT_N; j += numbers[i]) {
			marker_array[j] = false;
		}
		numbers[k] = numbers[i];
		++k;
	}
	numbers_count = k;
}

void feed_calculations_to_root() {
	//for (int j = 0; j < sieve_result.size(); ++j) numbers[j] = sieve_result[j];
	MPI_Send(numbers, numbers_count, MPI_INT, ROOT_RANK, NUMBERS_BUFF_TAG, MPI_COMM_WORLD);
}

int* distribution_arr;
void initialize_numbers_array() {
	numbers[0] = 5;
	numbers_count = 1;
	for (int k = 1; (6 * k + 5) < LIMIT_N; ++k) {
		numbers[numbers_count] = 6 * k + 1;
		numbers[numbers_count + 1] = 6 * k + 5;
		numbers_count += 2;
	}
	distribution_arr = new int[numbers_count];
	//distribution_arr = new int[(numbers_count + world_size - 1) / world_size + 1];
}

void distribute_numbers() {
	//MPI_Scatter(numbers, numbers_count, MPI_INT, nullptr, 0, MPI_INT, ROOT_RANK, MPI_COMM_WORLD);
	//std::sort(numbers, numbers + numbers_count);
	for (int i = 1; i < world_size; ++i) {
		int k = 0;
		for (int j = i - 1; j < numbers_count; j += world_size - 1, ++k) distribution_arr[k] = numbers[j];
		MPI_Send(distribution_arr, k, MPI_INT, i, NUMBERS_BUFF_TAG, MPI_COMM_WORLD);
	}
}

void nodes_receive_numbers() {
	MPI_Status numbers_status;
	MPI_Recv(numbers, LIMIT_N / 3 + 1, MPI_INT, ROOT_RANK, NUMBERS_BUFF_TAG, MPI_COMM_WORLD, &numbers_status);
	MPI_Get_count(&numbers_status, MPI_INT, &numbers_count);
}

void root_receive_data() {
	MPI_Status receive_status;
	numbers_count = 0;
	for (int i = 1; i < world_size; ++i) {
		int process_numbers_count;
		MPI_Recv(numbers + numbers_count, LIMIT_N / 3 + 1, MPI_INT, i, NUMBERS_BUFF_TAG, MPI_COMM_WORLD, &receive_status);
		MPI_Get_count(&receive_status, MPI_INT, &process_numbers_count);
		numbers_count += process_numbers_count;
	}
}

void root_notify_calculations_state(int ended) {
	for (int i = 1; i < world_size; ++i) {
		MPI_Send(&ended, 1, MPI_INT, i, CALCULATION_STATE_TAG, MPI_COMM_WORLD);
	}
}

bool nodes_fetch_notification() {
	int notification_data;
	MPI_Recv(&notification_data, 1, MPI_INT, ROOT_RANK, CALCULATION_STATE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	return notification_data;
}

std::vector<int> primes;
void sequential_sieve_algorithm(int n) {
	bool* is_prime = new bool[n];
	for (int i = 0; i < n; ++i) is_prime[i] = i % 2 == 1;
	is_prime[1] = false;
	is_prime[2] = true;
	for (int i = 3; i < n; i += 2) {
		if (is_prime[i]) {
			for (int j = 2 * i; j < n; j += i) {
				is_prime[j] = false;
			}
			primes.push_back(i);
		}
	}
}

void root_last_sieve() {
	int* marker_arr = new int[LIMIT_N];
	for (int i = 0; i < LIMIT_N; ++i) marker_arr[i] = true;
	int k = 0;
	std::sort(numbers, numbers + numbers_count);
	for (int i = 0; i < numbers_count; ++i) {
		if (!marker_arr[numbers[i]]) continue;
		for (int j = numbers[i] * 2; j < LIMIT_N; j += numbers[i]) marker_arr[j] = false;
		numbers[k] = numbers[i];
		++k;
	}
	numbers_count = k;
}

int main(int argc, char* argv[])
{
	double t1, t2, broadcast_duration = 0.0;
	MPI_Init(&argc, &argv);
	MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	// Get the rank of the process
	MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
	numbers = new int[(LIMIT_N + 2) / 3 + 1];

	if (process_rank == ROOT_RANK) {
		t1 = MPI_Wtime();
		initialize_numbers_array();
		std::cout << "ROOT : numbers initialized, " << numbers_count << " initial possible primes" << std::endl;
		int primes_count;
		for (int i = 0; i < 2; ++i) {
			root_notify_calculations_state(true);
			//std::cout << "ROOT : sent notification" << std::endl;
			distribute_numbers();
			//std::cout << "ROOT : numbers distributed" << std::endl;
			primes_count = numbers_count;
			root_receive_data();
			std::cout << "ROOT : calculations result received, reduced to " << numbers_count << " possible primes" << std::endl;
		} //while (primes_count != numbers_count);
		root_notify_calculations_state(false);
		root_last_sieve();

		t2 = MPI_Wtime();
		std::cout << "ROOT : execution ended" << std::endl;
		std::cout << "ROOT : parallel execution time : " << (t2 - t1) << std::endl;

		t1 = MPI_Wtime();
		sequential_sieve_algorithm(LIMIT_N);
		t2 = MPI_Wtime();
		std::cout << "ROOT : sequential execution time : " << (t2 - t1) << std::endl;
		std::cout << "ROOT : primes under " << LIMIT_N << " count is " << numbers_count << std::endl;
		std::cout << "ROOT : the first 100 primes >= 5 of course \n";
		std::sort(numbers, numbers + numbers_count);
		for (int i = 0; i < 100; ++i) std::cout << numbers[i] << "\n";
		//std::set<int> primes;
		//for (int i = 0; i < numbers_count; ++i) primes.insert(numbers[i]);
		//int i = 0;
		//std::cout << "prines count " << primes.size() << std::endl;
		//for (int k : primes) {
		//	++i;
		//	std::cout << k << std::endl;
		//	if (i == 100) break;
		//}
	}
	else {
		marker_array = new bool[LIMIT_N];
		for (int i = 0; i < LIMIT_N; ++i) marker_array[i] = true;
		marker_array[0] = marker_array[1] = false;
		while (nodes_fetch_notification()) {
			root_to_nodes_communication_delays -= MPI_Wtime();
			nodes_receive_numbers();
			root_to_nodes_communication_delays += MPI_Wtime();
			//std::cout << "NODE" << process_rank << " : received " << numbers_count << " numbers" << std::endl;
			local_sieve_computation_cost -= MPI_Wtime();
			run_sieve();
			local_sieve_computation_cost += MPI_Wtime();
			//std::cout << "NODE" << process_rank << " : completed sieve calculation, " << sieve_result.size() << " primes after sieve" << std::endl;
			node_to_root_communication_delays -= MPI_Wtime();
			feed_calculations_to_root();
			node_to_root_communication_delays += MPI_Wtime();
		}
		std::cout << "NODE" << process_rank << " : to root " << node_to_root_communication_delays << ", from ROOT " << root_to_nodes_communication_delays << std::endl;
		std::cout << "NODE" << process_rank << " communication penalty " << (node_to_root_communication_delays + root_to_nodes_communication_delays) << std::endl;
		std::cout << "NODE" << process_rank << " local sieve computation cost : " << local_sieve_computation_cost << std::endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Finalize();
}
