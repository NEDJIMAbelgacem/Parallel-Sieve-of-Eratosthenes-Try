// Sieve-of-Eratosthenes-parallel.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <omp.h>
#include <iostream>
#include <vector>

// kinda failed implementation 1 : from 10 to 30 % improvement

#define LIMIT_N (100000l * 1024)
#define NB_THREADS 4
#define NB_PRIMES 256

// Parallel implementation data
bool is_prime[LIMIT_N];
long long primes[LIMIT_N];
long long primes_count = 0;

long long to_be_sieved_primes[NB_PRIMES];
long long to_be_sieved_primes_count = 0;

// variable for analysing performance
double average_wait_time = 0.0;
double average_sieving_time = 0.0;
double sequential_section_time = 0.0;

int nb_threads = NB_THREADS;

void parallel_sieve()
{
	int thread_id, i, j;
	long long prime = 3, p;
	int start_index, end_index;
	double wait_time = 0.0;
	double sieving_time = 0.0;
	long long local_primes[NB_PRIMES];

	// initialize arrays
	for (long long i = 0; i < LIMIT_N; ++i) is_prime[i] = i % 2;
	for (long long i = 6; i < LIMIT_N; i += 3) is_prime[i] = false;
	for (long long i = 10; i < LIMIT_N; i += 5) is_prime[i] = false;
	is_prime[1] = false;
	is_prime[2] = true;
#pragma omp parallel firstprivate(nb_threads) default(shared) private(thread_id, i, j, p, start_index, end_index, sieving_time, wait_time, local_primes)
	{
		thread_id = omp_get_thread_num();
		start_index = thread_id * LIMIT_N / nb_threads;
		end_index = start_index + LIMIT_N / nb_threads;
		while (true) {

			sieving_time -= omp_get_wtime();
			for (i = thread_id; i < to_be_sieved_primes_count; i += nb_threads) {
				p = to_be_sieved_primes[i];
				for (j = 2 * p; j < LIMIT_N; j += p) is_prime[j] = false;
			}
			sieving_time += omp_get_wtime();

			wait_time -= omp_get_wtime();
#pragma omp barrier
			wait_time += omp_get_wtime();

			// sequential section
#pragma omp master
			{
				sequential_section_time -= omp_get_wtime();
				to_be_sieved_primes_count = 0;
				while (true)
				{
					if (prime >= LIMIT_N) break;
					to_be_sieved_primes[to_be_sieved_primes_count] = prime;
					to_be_sieved_primes_count++;
					if (to_be_sieved_primes_count >= NB_PRIMES) break;
					prime += 2;
					while (prime < LIMIT_N && (!is_prime[prime])) prime += 2;
				}
				sequential_section_time += omp_get_wtime();
			}

			wait_time -= omp_get_wtime();
#pragma omp barrier
			wait_time += omp_get_wtime();
			if (prime >= LIMIT_N) break;
		}

		average_wait_time += wait_time;
		average_sieving_time += sieving_time;
	}

	average_wait_time /= nb_threads;
	average_sieving_time /= nb_threads;

	primes[primes_count++] = 2;
	for (long long i = 1; i < LIMIT_N; i += 2) {
		if (is_prime[i]) primes[primes_count++] = i;
	}
}


long long seq_primes[LIMIT_N];
bool seq_is_prime[LIMIT_N];
long long seq_primes_count = 0;
void sequential_sieve()
{
	for (long long i = 0; i < LIMIT_N; ++i) seq_is_prime[i] = i % 2;
	seq_is_prime[2] = true;
	seq_is_prime[1] = false;
	seq_primes[seq_primes_count] = 2;
	seq_primes_count = 1;
	long long j;
	for (long long p = 3; p < LIMIT_N; p += 2) {
		if (seq_is_prime[p]) {
			for (j = 2 * p; j < LIMIT_N; j += p) {
				seq_is_prime[j] = false;
			}
			seq_primes[seq_primes_count] = p;
			seq_primes_count++;
		}
	}
}



bool test_identical_results() 
{
	if (primes_count != seq_primes_count) return false;
	for (long long i = 0; i < primes_count; ++i) {
		//std::cout << primes[i] << " " << seq_primes[i] << std::endl;
		if (primes[i] != seq_primes[i]) return false;
	}
	return true;
}

int main()
{
	std::cout << "n=" << LIMIT_N << std::endl;
	double sequential_time = -omp_get_wtime();
	sequential_sieve();
	sequential_time += omp_get_wtime();
	std::cout << "sequential sieve time : " << sequential_time << std::endl;

	omp_set_num_threads(NB_THREADS);
	double parallel_time = -omp_get_wtime();
	parallel_sieve();
	parallel_time += omp_get_wtime();


	std::cout << "parallel algorithm execution time : " << parallel_time << std::endl;
	std::cout << "speed up : " << (sequential_time / parallel_time) << std::endl;
	std::cout << "sequential section time " << sequential_section_time << std::endl;
	std::cout << "average wait time " << average_wait_time << std::endl;
	std::cout << "average sieving time " << average_sieving_time << std::endl;

	if (test_identical_results()) std::cout << "identical results between sequential and parallel algorithms" << std::endl;
	else {
		std::cout << "primes count : " << primes_count << std::endl;
		std::cout << "sequential primes count : " << seq_primes_count << std::endl;
		for (int i = 0; i < 20; ++i) {
			std::cout << primes[i] << " " << seq_primes[i] << std::endl;
		}
	}
	std::cin.get();
}
