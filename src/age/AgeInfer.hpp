//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef AgeInfer_hpp
#define AgeInfer_hpp

#include <algorithm>
#include <deque>
#include <iomanip>
#include <iostream>
#include <memory>
#include <mutex>
#include <unordered_set>

#include "Random.h"
#include "Threadpool.h"
#include "Clock.hpp"

#include "Gen.hpp"
#include "GenSample.hpp"
#include "GenMarker.hpp"
#include "GenGrid.hpp"
#include "GenShare.hpp"

#include "IBD.hpp"
#include "IBD_DGT.hpp"
#include "IBD_FGT.hpp"
#include "IBD_HMM.hpp"
#include "IBD_SIM.hpp"

#include "Age.hpp"
#include "AgeDensity.hpp"
#include "AgeEstimate.hpp"


namespace Age
{
	struct Pair;
	struct Site;
	
	
	// Pair of individuals/chromosomes
	struct Pair
	{
		using Data = std::shared_ptr< Pair >;
		using List = std::deque< Data >;
		
		
		// construct
		Pair(const Gamete::Pair, const bool);
		
		// print to file
		static void print_header(std::ostream &, const Param::Data, const bool);
		void print(std::ostream &, const Param::Data, const bool) const;
		
		
		std::weak_ptr<Site> site; // pointer to focal site
		Gamete::Pair        pair; // pair of individuals/chromosomes
		const bool       sharing; // indicate sharing of focal site
		
		decimal_t        missing; // calculated pairwise missing rate
		IBD::Segment     segment; // detected IBD segment
		IBD::SegDiff     segdiff; // mutational differences over segment
		
		std::array<CCF, n_clocks> ccf; // cumulative coalescent function
		
		bool done;
	};
	
	
	
	// Select nearest neighbours
	class Near
	{
	public:
		
		struct Rank
		{
			// construct
			Rank(const Gamete::Pair &, const size_t);
			Rank(const Rank &);
			Rank(Rank &&);
			
			const Gamete::Pair pair; // con/discordant pair of haplotypes
			const size_t       dist; // Hamming distance
			const size_t       rand; // random value for random sorting
			
			// sort
			bool operator < (const Rank &) const;
			bool operator > (const Rank &) const;
			
			
			using List = std::list<Rank>;
		};
		
		
		// construct
		Near(const size_t &, const Gen::Marker::Key &, const Gen::Grid::Data, const Param::Data);
		
		// perform all pairwise comparisons
		bool pairwise(const Param::Data);
		
		// apply filter to pairs
		//bool filter(const Param::Data);
		
		Rank::List concord;
		Rank::List discord;

		
	private:
		
		struct Chunk
		{
			// construct
			Chunk(const Gen::Marker::Key &, const Gamete &, const Gen::hap_vector_t &, const Param::Data); // ranked
			
			const Gamete      chr;
			Gen::hap_vector_t lhs;
			Gen::hap_vector_t rhs;
			
			// get Hamming distance
			size_t dist(const Chunk &) const;
			
			using List = std::vector<Chunk>;
		};
		
		struct Hold
		{
			Hold(Near *, const Gen::Sample::Key &, const Gen::Marker::Key &, const Gen::Grid::Data, const Param::Data);
			
			void run();
			
			Near *                 ptr;
			const Gen::Sample::Key key;
			const Gen::Marker::Key foc;
			const Gen::Grid::Data  grd;
			const Param::Data      par;
		};
		
		Chunk::List ins; // Haplotypes carrying the focal allele
		Chunk::List out; // All other haplotypes
		
		Threadpool< Hold > pool;
		std::mutex         lock;
	};
	
	
	
	// Focal site shred by subset of individuals/chromosomes
	struct Site
	{
		using Data = std::shared_ptr< Site >;
		using List = std::deque< Data >;
		
		// construct
		Site(const size_t &, const Gen::Marker::Key &, const Gen::Sample::Key::Vector &, const Gen::Grid::Data, const Param::Data);
		
		// construct for simulated results
		Site(const Gen::Marker::Key &, const IBD::SIM::Result::Data);
		
		
		// calculate allele frequencies
		void frequency(const Gen::Grid::Data);
		
		// estimate age
		void estimate(const Param::Data);
		
		// filter pairs based on CCF summary metric
		void filter(const ClockType, const Param::Data); // apply minimum pair exclusion threshold
		
		// print to file
		static void print_header(std::ostream &, const Param::Data, const bool);
		void print(std::ostream &, const size_t &, const bool, const bool) const;
		
		
		const size_t fk; // target fk
		const Gen::Marker::Key         focus; // focal site shared by pairs
		const Gen::Sample::Key::Vector share; // individuals sharing site
		
		IBD::FGT::Frequency::Data freq; // sharer allele frequencies
		
		Pair::List list;
		
		// composite likelihood estimates
		std::array<CLE, n_clocks> raw;
		std::array<CLE, n_clocks> adj;
		
		bool done;
		
		std::mutex guard;
	};
	
	
	// Exclude minimum number of pairs based on CCF summary metric
	struct MinExclude
	{
		decimal_t time; // time threshold
		size_t    ncon; // number of excluded concordant pairs
		size_t    ndis; // number of excluded discordant pairs
		decimal_t wsum; // weighted sum
		
		using List = std::vector<MinExclude>;
	};
	
	struct SortNode
	{
		const decimal_t node;
		const size_t    rand;
		const Pair::List::iterator pair;
		
		// construct
		SortNode(const decimal_t &, const Pair::List::iterator &);
		
		// sort
		bool operator <  (const SortNode &) const;
		bool operator  > (const SortNode &) const;
		
		using List = std::list<SortNode>;
	};
	
	
	// Result queue
	class Queue
	{
	public:
		
		struct Hold
		{
			// construct
			Hold(const size_t &, const Gen::Marker::Key &, const Gen::Sample::Key::Vector &);
			
			// construct for simulated results
			Hold(const Gen::Marker::Key &);
			
			size_t fk;
			Gen::Marker::Key site;
			Gen::Sample::Key::Vector share;
			
			using List = std::deque< Hold >;
		};
		
		
		// construct
		Queue(const Gen::Share::Data, const size_t, const Gen::Grid::Data, const Param::Data, const IBD::SIM::Result::Data = nullptr, const bool = false);
		
		// next batch
		size_t next(const IBD::SIM::Result::Data = nullptr);
		
		// total number of pairs
		size_t size() const;
		
		
		Pair::List pairs;
		Site::List sites;
		
		
	private:
		
		const Gen::Grid::Data  source;
		const Param::Data      param;
		const size_t           limit;
		
		size_t     total;
		Hold::List queue;
	};
	
	
	// Exection of segment detection and age inference
	class Infer
	{
	public:
		
		// constructs
		Infer(const Param::Data, const IBD::DetectMethod, const decimal_t, const Pair::Data, const Gen::Grid::Data, const IBD::HMM::Model::Data = nullptr);
		
		// execute detection
		void run();
		
		
	private:
		
		// determine chromosomes
		static Gen::ChrType chr_share(Gen::Variant &&);
		static Gen::ChrType chr_other(Gen::Variant &&);
		
		// determine segment differences
		void detect_segdiff(const Site::Data, const Gen::hap_vector_t, const Gen::hap_vector_t);
		void approx_segdiff(const Site::Data, const Gen::hap_vector_t, const Gen::hap_vector_t);
		//		void approx_segdiff_concordant(const Site::Data, const Gen::Variant::Vector::Data, const Gen::Variant::Vector::Data);
		//		void approx_segdiff_discordant(const Site::Data, const Gen::Variant::Vector::Data, const Gen::Variant::Vector::Data);
		
		// run methods
//		void fgt(const Site::Data, const Gen::Variant::Vector::Data, const Gen::Variant::Vector::Data);
//		void dgt(const Site::Data, const Gen::Variant::Vector::Data, const Gen::Variant::Vector::Data);
		void hmm(const Site::Data, const Gen::hap_vector_t, const Gen::hap_vector_t);
//		void sim(const Site::Data, const Gen::Variant::Vector::Data, const Gen::Variant::Vector::Data);
		
		
		const Param::Data param; // age estimation parameters
		const IBD::DetectMethod method; // chosen method
		const Pair::Data        target; // target pair
		const Gen::Grid::Data   source; // grid data source
		const IBD::HMM::Model::Data model; // HMM model
		const decimal_t max_missing_rate;
	};
}


#endif /* AgeInfer_hpp */

