//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef Age_hpp
#define Age_hpp

#include <memory>
#include <vector>

#include "Decimal.h"

#include "Gen.hpp"
#include "GenMarker.hpp"
#include "GenGrid.hpp"

#include "IBD.hpp"


namespace Age
{
	// Individual chromosome
//	struct Chromo
//	{
//		Gen::Sample::Key     individual;
//		mutable Gen::ChrType chromosome;
//		
//		using Pair = std::pair< Chromo, Chromo >;
//	};
	
	
	static constexpr int n_clocks = 3;
	
	enum ClockType : int
	{
		MUT_CLOCK = 0,
		REC_CLOCK = 1,
		CMB_CLOCK = 2
	};
	
	
	struct Gamete
	{
		Gen::Sample::Key individual;
		Gen::ChrType     chromosome;
		
		// construct
		Gamete();
		Gamete(const Gen::Sample::Key &, const Gen::ChrType);
		Gamete(const Gamete &);
		Gamete(Gamete &&);
		
		// assign
		Gamete & operator = (const Gamete &);
		Gamete & operator = (Gamete &&);
		
		// sort
		bool operator <  (const Gamete &) const;
		bool operator  > (const Gamete &) const;
		bool operator == (const Gamete &) const;
		bool operator != (const Gamete &) const;
		
		using Pair = std::pair<Gamete, Gamete>;
	};
	
	
	// Cumulative coalescent function
	struct CCF
	{
		// construct
		CCF();
		CCF(const CCF &);
		CCF(CCF &&);
		
		// assign
		CCF & operator = (const CCF &);
		CCF & operator = (CCF &&);
		
		size_t shape;
		decimal_t rate;
		
		decimal_t q25;
		decimal_t q50;
		decimal_t q75;
		
		decimal_vector_t d;
		
		bool good;
		bool pass;
	};
	
	
	// Composite likelihood estimate
	struct CLE
	{
		// construct
		CLE();
		
		size_t n_shared;
		size_t n_others;
		
		decimal_t mean;
		decimal_t mode;
		decimal_t median;
		decimal_t ci_025;
		decimal_t ci_975;
		decimal_t estim;
		decimal_t lower;
		decimal_t upper;
		
		decimal_vector_t grid;
		
		bool good;
	};
	
	
	// Estimation parameters
	struct Param
	{
		using Data = std::shared_ptr< Param >;
		
		
		// construct
		Param(const Gen::Grid::Data, const size_t = 10000, const decimal_t = 1e-08, const size_t = 1024, const decimal_t = 40.0, const bool = true);
		
		// set theta manually
		void set_theta(const double);
		
		// estimate theta (Watterson estimator)
		void est_theta(const Gen::Grid::Data);
		
		
		const size_t    nt; // length of times prior
		const size_t    Ng; // sample size, genotypes
		const size_t    Nh; // sample size, haplotypes
		const size_t    Nm; // Number of markers
		const decimal_t Ne; // effective population size
		const decimal_t Mr; // mutation rate
		
		decimal_t theta; // theta per site
		bool theta_est;  // flag theta is Watterson estimator
		bool theta_man;  // flag theta is manually set
		
		
		// settings
		
		bool include_prior; // coalescent times prior
		// bool use_mut_clock; // mutation clock
		// bool use_rec_clock; // recombination clock
		bool use_post_prob; // include posterior probabilities (HMM)
		bool use_hard_brks; // assume uncertain recombination breakpoints
		
		bool run_mut_clock;
		bool run_rec_clock;
		bool run_cmb_clock;
		
		//bool apply_filter_fixed;
		//bool apply_filter_detect;
		bool apply_nearest_neighb;
		bool relax_nearest_neighb;
		bool use_tree_consistency; // count pairwise differences only if allele count <= focal fk
		
		size_t limit_sharers; // max. number of sharer pairs, randomly subsamples
		size_t outgroup_size; // max. number of other individuals/chromosomes that do not share the focal mutation
		
		size_t limit_sharers_max;
		size_t outgroup_size_max;
		
		size_t breakpt_range; // max number of considered sites after breakpoint
		
		size_t nearest_range; // max. range LHS/RHS to focal site
		
		size_t threads;
		
		
		// variables
		
		IBD::Segment     boundary;
		decimal_vector_t position;
		decimal_vector_t distance;
		decimal_vector_t frequency;
		decimal_vector_t log_het;
		decimal_vector_t log_hom;
		decimal_vector_t cum_log_hom;
		
		decimal_vector_t prior; // defined time grid
		decimal_vector_t log_prior; // log scale
	};
}


#endif /* Age_hpp */

