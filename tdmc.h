// A reasonably generic implementation of the Ticketed
// Diffusion Monte Carlo (TDMC) algorithm introduced by
// Hairer and Weare (DOI: 10.1002/cpa.21526).
// Meant for educational purposes.
// If you use this code please cite it using its
// DOI on Github and Zenodo; find it in the README 
// of www.github.com/jfdama/ticketed-dmc .
// (c) James Farris Dama 2015

#ifndef TDMC_H
#define TDMC_H

#include <list>
#include <vector>
#include <utility>
#include <tuple>

#include <random>
#include <cmath>
#include <algorithm>

#include <iostream>
#include <ctime>
#include <cstdint>

// The TDMC algorithm is written in terms of an arbitrary
// dynamical sampling function, P in the paper and
// "sample_dynamics" here, and an arbitrary stepwise
// reweighting function, χ in the paper and "chi" here. 

// These templates allow one to define TDMC on walkers with any
// type of state T, sampling dynamics P:T->T or P:(T, timestep)->T, 
// and χ:(T, T)->double or χ:(T, T, timestep)->double that are 
// already defined at compile time. The specialized function then 
// takes as arguments an initial walker state, a number of steps to 
// run, and a number of walkers to initialize with, and after running 
// for the desired number of steps it returns a vector of pairs of 
// the final samples and the integer ids of the replicates they 
// correspond to.

// Implementation for P:(T, timestep)->T and χ:(T, T, timestep)->double.
template <typename T, void (*propagate_sample)(T&, int), double (*chi)(const T&, const T&, int)>
std::vector< std::pair<T, int> > run_tdmc(T initial_walker_state, int n_dmc_steps, int n_initial_walkers) {

	// Initialize clock and random number generator.
	std::cout << "Beginning TDMC." << std::endl;
    clock_t start_clock = clock();
    uint64_t n_dynamics_evaluations = 0;
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<> standard_uniform_dist(0, 1);

    // Define a walker data structure: walkers have a state,
    // a ticket value, and correspond to one of the initial
    // walkers given by 'replicate_id'.
    struct walker_data {
    	T state;
    	double ticket;
    	int replicate_id;
    };

	// Create a list of walkers with uniformly distributed tickets.
	auto walker_list = std::list< walker_data >();
	for (int i = 0; i < n_initial_walkers; i++) {
		walker_data new_walker_data = {initial_walker_state, standard_uniform_dist(gen), i};
		walker_list.push_back(new_walker_data);
	}

	// Perform ticketed DMC sampling for n_dmc_steps steps.
	for (int i = 0; i < n_dmc_steps; i++) {
		// At every step, advance and adjust the copy number of each walker.
		// Iterate over all walkers starting from the first.
		auto curr_walker = walker_list.begin();
		while (curr_walker != walker_list.end()) {
			// Calculate a new state, saving the old & new outside the
			// list for just a moment.
			T initial_walker_state = curr_walker->state;
			propagate_sample(curr_walker->state, i); n_dynamics_evaluations++;
			// Calculate the generalized DMC weight for the proposed step.
			double step_weight = exp(-chi(curr_walker->state, initial_walker_state, i));
			// If the weight is lower than the walker's ticket, delete the walker
			// and advance to the next.
			if (step_weight < curr_walker->ticket) {
				curr_walker = walker_list.erase(curr_walker);
			// Otherwise update the walker's ticket value and clone it as needed.
			} else {
				int replicate_id = curr_walker->replicate_id;
				int n_clones_needed = std::max(1, (int) (step_weight + standard_uniform_dist(gen)));
				// Update the ticket value by dividing it by the step weight.
				curr_walker->ticket = curr_walker->ticket / step_weight;
				// Add new clones as needed.
				int curr_n_clones = 1;
				while (curr_n_clones < n_clones_needed) {
					std::uniform_real_distribution<> new_ticket_dist(1.0 / step_weight, 1);
					walker_data clone_walker_data = {curr_walker->state, new_ticket_dist(gen), replicate_id};
					curr_walker = walker_list.insert(++curr_walker, clone_walker_data);
					curr_n_clones++;
				}
				// Advance to the next walker.
				++curr_walker;
			}
		}
	}

	// Print completion log message.
	std::cout << "Completed TDMC in " << (double)(clock() - start_clock)/CLOCKS_PER_SEC << " seconds using " << n_dynamics_evaluations << " evaluations of the sampling dynamics." << std::endl;

	// Convert to the return type and finish.
	std::vector< std::pair<T, int> > final_walker_states;
	final_walker_states.reserve(walker_list.size());
	for (auto walker : walker_list) {
		final_walker_states.push_back(std::make_pair(walker.state, walker.replicate_id));
	}
	return final_walker_states;
}

// Template adaptors for making sampling dynamics P:T->T and
// birth-death weight functions χ:(T, T)->double feed into the
// P:(T, timestep)->T and χ:(T, T, timestep)->double template 
// just above.

// Dress up P:T->T as a function P:(T, timestep)->T with no 
// actual timestep dependence.
template <typename T, void (*propagate_sample)(T&) > void int_dressed_propagate_sample(T &x, int step) {propagate_sample(x);}
// Dress up χ:(T, T)->double as a function
// χ:(T, T, timestep)->double with no actual timestep dependence.
template <typename T, double (*chi)(const T&, const T&)> double int_dressed_chi(const T &x, const T &y, int step) {return chi(x, y);}

// If neither dynamics nor birth-death process depends on input, 
// this template is used.
template <typename T, void (*propagate_sample)(T&), double (*chi)(const T&, const T&)>
std::vector< std::pair<T, int> > run_tdmc(T initial_walker_state, int n_dmc_steps, int n_initial_walkers) {
	return run_tdmc<T, int_dressed_propagate_sample<T, propagate_sample>, int_dressed_chi<T, chi> >(initial_walker_state, n_dmc_steps, n_initial_walkers);
}

// If only the dynamics depends on timestep, this template is used.
template <typename T, void (*propagate_sample)(T&, int), double (*chi)(const T&, const T&)>
std::vector< std::pair<T, int> > run_tdmc(T initial_walker_state, int n_dmc_steps, int n_initial_walkers) {
	return run_tdmc<T, propagate_sample, int_dressed_chi<T, chi> >(initial_walker_state, n_dmc_steps, n_initial_walkers);
}

// If only the birth-death process depends on timestep, this template 
// is used.
template <typename T, void (*propagate_sample)(T&), double (*chi)(const T&, const T&, int)>
std::vector< std::pair<T, int> > run_tdmc(T initial_walker_state, int n_dmc_steps, int n_initial_walkers) {
	return run_tdmc<T, int_dressed_propagate_sample<T, propagate_sample>, chi >(initial_walker_state, n_dmc_steps, n_initial_walkers);
}

#endif