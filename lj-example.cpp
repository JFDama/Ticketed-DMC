// 2D LJ7 cluster example for the Ticketed Diffusion
// Monte Carlo (TDMC) algorithm.
// (c) James Farris Dama 2015

#include "tdmc.h"
#include <iostream>
#include <array>
#include <vector>
#include <random>

using namespace std;

//--------------------------------------------------------------
// Initial definitions and declarations
//--------------------------------------------------------------

// Define input variables as used in the paper.
constexpr double paper_gamma = 0.4;
constexpr double paper_lambda = 1.9;
constexpr double paper_epsilon = 0.001;
constexpr int paper_nsteps = (int) (2 / paper_epsilon);
// Choose a number of replicates to perform.
constexpr int n_runs = 10000;

// Define the number of particles and problem dimension.
// This code has been written so as to make adjusting 
// this example to study the analogous rare event in 
// an LJ 13 cluster in 3D a natural exercise.
constexpr int n_particles = 7;
constexpr int dimension = 2;
constexpr int n_vars = n_particles * dimension;
// Define the state type for these simulations.
using lj_state_vec = std::array<double, n_vars>;
// Simple function for printing these states as a line.
void print_lj_state(lj_state_vec &x) {for (double real : x) cout << real << " "; cout << endl;}

// Two system physics convenience functions defined in full
// below main.
lj_state_vec calc_lj_force(const lj_state_vec &x);
lj_state_vec calc_random_displacement();

// Dynamics propagation function. Advance state according to
// Brownian dynamics on the LJ potential surface for one time step.
void propagate_state(lj_state_vec &x) {
	lj_state_vec conservative_force = calc_lj_force(x);
	lj_state_vec random_displacement = calc_random_displacement();
	for (int i = 0; i < n_vars; i++) x[i] += conservative_force[i] * paper_epsilon + random_displacement[i];
}

// Convenience function for calculating a birth/death process
// intermediate variable defined in full below main.
double paper_V(const lj_state_vec &x);

// DMC walker birth and death determining function.
double dmc_chi(const lj_state_vec &x_new, const lj_state_vec &x_old) {
	return paper_V(x_new) - paper_V(x_old);
}

// Brute-force-as-DMC birth and death determining function.
double trivial_dmc_chi(const lj_state_vec &x_new, const lj_state_vec &x_old) {
	return 0;
}

// Convenience functions for analysis: calculating the 
// rearrangement coordinate, centering a cluster, and
// calculating which particle is closest to the center
// of the cluster.
double rearrangement_coord(const lj_state_vec &x);
void center_cluster(lj_state_vec &x);
int central_particle(const lj_state_vec &x);

//--------------------------------------------------------------
// Main program
//--------------------------------------------------------------

// Run the 2D LJ7 example using TDMC.
int main(void) {
	cout << "Reproducing the Hairer and Weare seven particle LJ cluster example." << endl;
	cout << "Using parameters γ = " << paper_gamma << ", λ = " << paper_lambda << ", and ε = " << paper_epsilon << "." << endl;
	// Initialize 7 particles, 6 in a hexagon around 
	// a single central particle.
	lj_state_vec initial_state = {{0.0, 0.0, -1.0, 0.0, 1.0, 0.0, -0.5, 0.866025, 0.5, 0.866025, -0.5, -0.866025, 0.5, -0.866025}};
	// Scale so that the particles are close to the 
	// sixfold symmetric energy minimum.
	for (int i = 0; i < n_vars; i++) {
		initial_state[i] *= 1.11846;
	}
	// Print the initial state and the initial LJ forces on the particles.
	cout << "Initial state:" << endl;
	print_lj_state(initial_state);
	lj_state_vec f = calc_lj_force(initial_state);
	cout << "Initial forces:" << endl;
	print_lj_state(f);

	// Run TDMC and brute force simulation n_runs times each.
	auto final_walker_data = run_tdmc<lj_state_vec, propagate_state, dmc_chi>(initial_state, paper_nsteps, n_runs);
	auto final_bf_walker_data = run_tdmc<lj_state_vec, propagate_state, trivial_dmc_chi>(initial_state, paper_nsteps, n_runs);

	// Calculate the estimators for TDMC.
	// Calculate means and variances of indicators
	// of existence and occupation of regions B and D.
	// Region B: a particle other than the initial central
	// particle is within .1σ of the cluster center.
	// Region D: a particle other than the initial central
	// particle is closest to the cluster center.
	// Region D contains Region B for all practical purposes.
	double I_B_estimator_sum = 0.0;
	double I_D_estimator_sum = 0.0;
	double one_estimator_sum = 0.0;
	double I_B_estimator_sumsq = 0.0;
	double I_D_estimator_sumsq = 0.0;
	double one_estimator_sumsq = 0.0;

	int curr_replicate = 0;
	double I_B_estimator = 0.0;
	double I_D_estimator = 0.0;
	double one_estimator = 0.0;
	for (auto walker_datum :  final_walker_data) {
		// If we've run out of walkers for the current replicate,
		// gather the statistics for that replicate and move on
		// to the next.
		if (walker_datum.second != curr_replicate) {
			I_B_estimator_sum += I_B_estimator;
	    	I_D_estimator_sum += I_D_estimator;
	    	one_estimator_sum += one_estimator;
	    	I_B_estimator_sumsq += I_B_estimator * I_B_estimator;
	    	I_D_estimator_sumsq += I_D_estimator * I_D_estimator;
	    	one_estimator_sumsq += one_estimator * one_estimator;
			curr_replicate = walker_datum.second;
			I_B_estimator = 0.0;
			I_D_estimator = 0.0;
			one_estimator = 0.0;
		}
		// Calculate the expectation of indicator functions by
		// weighting the walker by its final state bias as described
		// in the paper, then multiplying by the indicator of interest.
		// Sum over all the walkers for each replicate.
		double weight = exp(dmc_chi(walker_datum.first, initial_state));
		one_estimator += weight;
		if (rearrangement_coord(walker_datum.first) < 0.1) {
			I_B_estimator += weight;
		}
		if (central_particle(walker_datum.first) != 0) {
			I_D_estimator += weight;
		}
	}
		
	// Calculate the estimators for brute force simulation.
	// The expectation for existence is one in this case, and
	// an indicator value is equal to its square so don't bother
	// gathering both separately.
	double bf_I_B_estimator_sum = 0.0;
	double bf_I_D_estimator_sum = 0.0;
	for (auto bf_walker_datum :  final_bf_walker_data) {
		// Calculate the indicators.
		if (rearrangement_coord(bf_walker_datum.first) < 0.1) {
			bf_I_B_estimator_sum += 1.0;
		}
		if (central_particle(bf_walker_datum.first) != 0) {
			bf_I_D_estimator_sum += 1.0;
		}
	}

	// Print the results.
	cout << "Final estimate for E[1]: " << one_estimator_sum / n_runs  << ", variance " << (one_estimator_sumsq / n_runs - (one_estimator_sum / n_runs) * (one_estimator_sum / n_runs)) << endl;
	cout << "Final estimate for E[I_B]: " << I_B_estimator_sum / n_runs  << ", variance " << (I_B_estimator_sumsq / n_runs - (I_B_estimator_sum / n_runs) * (I_B_estimator_sum / n_runs)) << endl;
	cout << "Final estimate for E[bf_I_B]: " << bf_I_B_estimator_sum / n_runs  << ", variance " << (bf_I_B_estimator_sum / n_runs - (bf_I_B_estimator_sum / n_runs) * (bf_I_B_estimator_sum / n_runs)) << endl;
	cout << "Final estimate for E[I_D]: " << I_D_estimator_sum / n_runs  << ", variance " << (I_D_estimator_sumsq / n_runs - (I_D_estimator_sum / n_runs) * (I_D_estimator_sum / n_runs)) << endl;
	cout << "Final estimate for E[bf_I_D]: " << bf_I_D_estimator_sum / n_runs  << ", variance " << (bf_I_D_estimator_sum / n_runs - (bf_I_D_estimator_sum / n_runs) * (bf_I_D_estimator_sum / n_runs)) << endl;

	return 0;
}

//--------------------------------------------------------------
// Dynamics and paper_V convenience function definitions
//--------------------------------------------------------------

// Calculate the total forces experienced by each particle
// in each direction as a sum of LJ pair forces (in LJ units).
lj_state_vec calc_lj_force(const lj_state_vec &x) {
	// Start from having zero forces.
	lj_state_vec force = {};
	// Calculate the pair forces.
	for (int i = 0; i < n_particles; i++) {
		for (int j = i + 1; j < n_particles; j++) {
			// Calculate the pair displacement and squared distance.
			double distsq = 0;
			array<double, dimension> disp = {};
			for (int k = 0; k < dimension; k++) {
				disp[k] = x[dimension * i + k] - x[dimension * j + k];
				distsq += disp[k] * disp[k];
			}
			// Use the displacement and distance to calculate forces.
			for (int k = 0; k < dimension; k++) {
				force[dimension * i + k] += (48.0 * pow(distsq, -7) - 24.0 * pow(distsq, -4)) * disp[k];
				force[dimension * j + k] += (48.0 * pow(distsq, -7) - 24.0 * pow(distsq, -4)) * -disp[k];
			}
		}
	}
	return force;
}

// Calculate a random displacement given the timestep 
// and temperature.
lj_state_vec calc_random_displacement() {
	// Seed the RNG and set up a Gaussian number generator
	// prior to the first call of this function. 
	// Static variables are used in lieu of a functor class
	// to make the code a little less indirect for educational
	// purposes.
	static random_device rd;
	static mt19937_64 gen(rd());
	static normal_distribution<> brownian_increment(0, sqrt(2 * paper_gamma * paper_epsilon));
	
	// Start from zero.
	lj_state_vec disp = {};
	for (int i = 0; i < n_vars; i++) {
		disp[i] = brownian_increment(gen);
	}
	return disp;
}

// Calculate the sampling-biasing function defined for this
// system.
double paper_V(const lj_state_vec &x) {
	return (paper_lambda / paper_gamma) * rearrangement_coord(x);
}

//--------------------------------------------------------------
// Analysis convenience function definitions
//--------------------------------------------------------------

// Calculate the minimum of the distances between the cluster
// center of mass and all of the particles that did not begin
// as the central particle of the cluster.
double rearrangement_coord(const lj_state_vec &x) {
	// Calculate the center position of the cluster.
	array<double, dimension> center_position = {};
	for (int i = 0; i < n_particles; i++) {
		for (int j = 0; j < dimension; j++) {
			center_position[j] += x[dimension * i + j];
		}
	}
	for (int i = 0; i < dimension; i++) {
		center_position[i] /= n_particles;
	}
	// Find the minimum squared distance from one of the 
	// initially peripheral particles to the cluster center.
	// This differs from the paper only in that here the
	// particles are indexed from 0 rather than from 1.
	double min_center_dist_sq;
	for (int i = 1; i < n_particles; i++) {
		double dist_sq = 0;
		for (int j = 0; j < dimension; j++) {
			double disp = (x[dimension * i + j] - center_position[j]);
			dist_sq += disp * disp;
		}
		if (i == 1) {
			min_center_dist_sq = dist_sq;
		} else {
			min_center_dist_sq = min(min_center_dist_sq, dist_sq);	
		}
	}
	return sqrt(min_center_dist_sq);
}

// Center a cluster state vector.
void center_cluster(lj_state_vec &x) {
	// Calculate the center position of the cluster.
	array<double, dimension> center_position = {};
	for (int i = 0; i < n_particles; i++) {
		for (int j = 0; j < dimension; j++) {
			center_position[j] += x[dimension * i + j];
		}
	}
	for (int i = 0; i < dimension; i++) {
		center_position[i] /= n_particles;
	}
	// Subtract that from every particle's position.
	for (int i = 0; i < n_particles; i++) {
		for (int j = 0; j < dimension; j++) {
			x[dimension * i + j] -= center_position[j];
		}
	}
}

// Return the index of the particle closest to a cluster center.
int central_particle(const lj_state_vec &x) {
	// Make a local copy.
	auto y = x;
	// Center the copied cluster.
	center_cluster(y);
	// Find the particle closest to (0, 0).
	double min_center_dist_sq;
	int center_particle_index;
	for (int i = 0; i < n_particles; i++) {
		// Calculate the squared distance from zero.
		double dist_sq = 0;
		for (int j = 0; j < dimension; j++) {
			dist_sq += y[dimension * i + j] * y[dimension * i + j];
		}
		// If this is the first or the lowest squared distance
		// seen so far, set it as the new minimum and
		// record that this particle is closest to the center.
		if (i == 0 || dist_sq < min_center_dist_sq) {
			min_center_dist_sq = dist_sq;
			center_particle_index = i;
		}
	}
	return center_particle_index;
}