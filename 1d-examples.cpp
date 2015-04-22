// Simple 1d examples for the Ticketed Diffusion
// Monte Carlo (TDMC) algorithm.
// (c) James Farris Dama 2015

#include "tdmc.h"
#include <iostream>
#include <vector>
#include <random>

using namespace std;

// Make a particle drift deterministically in 1D.
void drift_dynamics(double &x) {
	x += .1;
}

// Make a particle undergo pure Brownian dynamics in 1D.
constexpr double timestep = .0001;
void brownian_dynamics(double &x) {
	static random_device rd;
	static mt19937_64 gen(rd());
	static normal_distribution<> brownian_increment(0, sqrt( 2 * timestep));
	x += brownian_increment(gen);
}

// A rare event sampling type weight.
constexpr double beta = -1.25;
double dmc_chi(const double &x_new, const double &x_old) {
	return  beta * (x_new - x_old);
}

int main(void) {
	// Run with drift dynamics for a single step.
	auto final_walker_data = run_tdmc<double, drift_dynamics, dmc_chi>(0.0, 1, 5);
	cout << "Returned walker states:";
	for (auto walker_data : final_walker_data) {
		cout << " " << walker_data.first;
	}
	cout << endl;
	// Run with Brownian dynamics for 1 time unit n_replicates times.
	int n_replicates = 10000;
	int n_timesteps = (int) (1 / timestep);
	auto final_walker_data2 = run_tdmc<double, brownian_dynamics, dmc_chi>(0.0, n_timesteps, n_replicates);
	cout << "Final walker number: " << final_walker_data2.size() << endl;
	cout << "Average final walkers per initial walker: " << (double) final_walker_data2.size() / n_replicates << endl;
	return 0;
}