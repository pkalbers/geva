//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef IBD_HMM_hpp
#define IBD_HMM_hpp

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <list>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

#include "Decimal.h"

#include "IBD.hpp"

#include "Gen.hpp"
#include "GenMarker.hpp"


namespace IBD
{
	// Hidden Markov Model
	namespace HMM
	{
		static const std::string name = "Hidden Markov Model (HMM)";
		static const std::string abbr = "HMM";


		// Possible hidden states

		static constexpr unsigned hidden_state_n = 2;

		enum HiddenState : uint_fast8_t
		{
			NON_STATE = 0,
			IBD_STATE = 1,
			STATE_UNDEF = 2
		};

		using path_vector_t = std::vector< HiddenState >;
		using path_distr_t = Side< path_vector_t >;

		using prob_t = std::array< decimal_t, hidden_state_n >;
		using prob_vector_t = std::vector< prob_t >;
		using prob_distr_t = Side< prob_vector_t >;

		using distr_n = Side< size_t >;


		inline distr_n size_distr(const Gen::Marker::Key & focal, const size_t & length)
		{
			return distr_n(focal.value + 1, length - focal.value);
		}

		inline path_distr_t init_path_distr(const distr_n & size)
		{
			return path_distr_t(path_vector_t(size[LHS], STATE_UNDEF), path_vector_t(size[RHS], STATE_UNDEF));
		}

		inline prob_distr_t init_prob_distr(const distr_n & size, const decimal_t value = decimal_nil)
		{
			const prob_t init = {value, value};
			return prob_distr_t(prob_vector_t(size[LHS], init), prob_vector_t(size[RHS], init));
		}


		// Observed genotype pair

		// static constexpr unsigned obs_gen_pair_n = 6;
		static constexpr unsigned obs_hap_pair_n = 3;

		// enum ObsGenPair : uint_fast8_t
		// {
		// 	G00 = 0,
		// 	G01 = 1,
		// 	G02 = 2,
		// 	G11 = 3,
		// 	G12 = 4,
		// 	G22 = 5,
		// 	G__ = 6
		// };

		enum ObsHapPair : uint_fast8_t
		{
			H00 = 0,
			H01 = 1,
			H11 = 2,
			H__ = 3
		};

		// using obs_vector_t = std::vector< ObsGenPair >;
		using obs_vector_t = std::vector< ObsHapPair >;


		// template < ObsGenPair O >
		// constexpr bool is_obstype(const ObsGenPair o)
		// {
		// 	return (O == o);
		// }

		template < ObsHapPair O >
		constexpr bool is_obstype(const ObsHapPair o)
		{
			return (O == o);
		}


		// constexpr ObsGenPair obs_gen_pair(const Gen::gen_t g0, const Gen::gen_t g1)
		// {
		// 	using namespace Gen;
		//
		// 	return
		// 	(is_genotype<G0>(g0)) ? ((is_genotype<G0>(g1)) ? G00: (is_genotype<G1>(g1)) ? G01: (is_genotype<G2>(g1)) ? G02: G__):
		// 	(is_genotype<G1>(g0)) ? ((is_genotype<G0>(g1)) ? G01: (is_genotype<G1>(g1)) ? G11: (is_genotype<G2>(g1)) ? G12: G__):
		// 	(is_genotype<G2>(g0)) ? ((is_genotype<G0>(g1)) ? G02: (is_genotype<G1>(g1)) ? G12: (is_genotype<G2>(g1)) ? G22: G__):
		// 	G__;
		// }

		constexpr ObsHapPair obs_hap_pair(const Gen::hap_t h0, const Gen::hap_t h1)
		{
			using namespace Gen;

			return
			(is_haplotype<H0>(h0)) ? ( (is_haplotype<H0>(h1)) ? H00: (is_haplotype<H1>(h1)) ? H01: H__ ):
			(is_haplotype<H1>(h0)) ? ( (is_haplotype<H0>(h1)) ? H01: (is_haplotype<H1>(h1)) ? H11: H__ ):
			H__;
		}

		// inline obs_vector_t obs_gen_pair_vector(const Gen::Variant::Vector::Data v0, const Gen::Variant::Vector::Data v1)
		// {
		// 	using namespace Gen;
		//
		// 	const size_t size = v0->size();
		//
		// 	if (size == 0 || size != v1->size())
		// 	{
		// 		return obs_vector_t(0);
		// 	}
		//
		// 	obs_vector_t ov(size);
		//
		// 	for (size_t i = 0; i < size; ++i)
		// 	{
		// 		ov[i] = obs_gen_pair(v0->gen(i), v1->gen(i));
		// 	}
		//
		// 	return ov;
		// }

		inline obs_vector_t obs_hap_pair_vector(const Gen::hap_vector_t & h0, const Gen::hap_vector_t & h1)
		{
			using namespace Gen;

			const size_t size = h0.size();

			if (size == 0 || size != h1.size())
			{
				return obs_vector_t(0);
			}

			obs_vector_t ov(size);

			for (size_t i = 0; i < size; ++i)
			{
				ov[i] = obs_hap_pair(h0[i], h1[i]);
			}

			return ov;
		}



		// Model probabilities
		class Model
		{
		public:

			using Data = std::shared_ptr<Model>;

			using inits_type = std::array< decimal_t, hidden_state_n >;
			using emiss_type = std::array< std::array< decimal_t, obs_hap_pair_n >, hidden_state_n >;
			using dists_type = decimal_t;
			using trans_type = std::array< std::array< decimal_t, hidden_state_n >, hidden_state_n >;

			using inits_list = std::vector< inits_type >;
			using emiss_list = std::vector< emiss_type >;
			using dists_list = std::vector< dists_type >;
			using trans_list = std::vector< trans_type >;

			using trans_map = std::unordered_map< size_t, trans_list >;


			// construct
			Model(const size_t &, const size_t &, inits_list &&, inits_list &&, emiss_list &&, dists_list &&);

			// generate transitions for target
			void prepare_transition(const size_t &);

			// return probabilities
			inits_list const & initial_con();
			inits_list const & initial_dis();
			emiss_list const & emission();
			trans_list const & transition(const size_t &);

			// return genetic distance
			decimal_t distance(const size_t &) const;

			// calculate transition matrix
			static decimal_t  calc_expected_age(const size_t &, const size_t &);
			static trans_type calc_trans_matrix(const size_t &, const size_t &, const size_t &, const decimal_t &);


			const size_t Ne; // effective size
			const size_t Nh; // sample size (haplotypes)


			bool do_iterative;
			
		private:

			inits_list inits_con; // concordant initial probabilties
			inits_list inits_dis; // discordant initial probabilties
			emiss_list emiss; // emission probabilties
			dists_list dists; // intervariant distances
			trans_map  trans; // transition probabilties per fk

			std::mutex guard;
		};


		// Hidden Markov Model algorithm
		class Algorithm
		{
		public:

			// construct
			Algorithm(const Gen::hap_vector_t, const Gen::hap_vector_t, const Model::Data, const bool);
			
			// decode segment from Viterbi path
			Segment detect(const size_t &, const Gen::Marker::Key &);
//			Segment detect_iterative(const size_t &, const Gen::Marker::Key &);

			// return Viterbi path
			prob_distr_t const & viterbi() const;
			Distribution const & weights() const;
			path_distr_t const & path() const;

			// return posterior probabilities
			prob_distr_t const & posterior(const bool);
			decimal_vector_t     posterior(const LR, const HiddenState, const bool);
			Distribution const & fwd_weights() const;
			Distribution const & bwd_weights() const;

			// print to stream
			void print(const Gen::Sample::Key::Pair &, const Gen::Marker::Key &, std::ostream & = std::cout) const;
			static void print_header(std::ostream & = std::cout);


		private:

			// execute main algorithm
			void execute_viterbi(const size_t &, const LR);
			void execute_viterbi_switch(const size_t &, const LR);
			void execute_forward(const size_t &, const LR);
			void execute_backward(const size_t &, const LR);
			void execute_posterior(const size_t, const LR, const bool);


			const obs_vector_t obs;  // full observation sequence
			const Model::Data  model;  // shared pointer to model parameters
			const size_t size; // length of sequence


			size_t fk; // current focal allele count
			Gen::Marker::Key focal; // current focal site

			//Side< Chain > chain; // Markov chain

			bool good; // indicate completion of viterbi algorithm
			bool post; // indicate completion of posterior probabilities


			path_distr_t v_path; // Viterbi path
			prob_distr_t v_prob; // Viterbi scaled probabilities
			Distribution v_wght; // Viterbi weights from scaling

			prob_distr_t f_prob; // forward scaled probabilities
			Distribution f_wght; // forward weights from scaling

			prob_distr_t b_prob; // backward scaled probabilities
			Distribution b_wght; // backward weights from scaling

			prob_distr_t p_prob; // posterior probabilities
			
			
//			Gen::hap_vector_t focal_A, other_A;
//			Gen::hap_vector_t focal_B, other_B;
			const bool is_discord;
		};
	}
}


#endif /* IBD_HMM_hpp */
