//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef GenShare_hpp
#define GenShare_hpp

#include <algorithm>
#include <iterator>
#include <map>
#include <memory>
#include <set>
#include <utility>
#include <vector>

#include "Random.h"

#include "GenSample.hpp"
#include "GenMarker.hpp"
#include "GenGrid.hpp"
#include "Binary.hpp"


namespace Gen
{
	// Index of target sites and pairs
	class Share
	{
	public:
		
		using Data = std::shared_ptr< Share >;
		
		using target_t = std::set< size_t >;
		using region_t = std::pair< size_t, size_t >;
		
		static constexpr size_t minimum = 2;
		
		
		// Share index, by pair and site
		struct Index
		{
			using Sites = std::map< Marker::Key, Sample::Key::Vector >;
			using Pairs = std::map< Sample::Key::Pair, Marker::Key::Set >;
			
			
			size_t fk;
			
			Sites sites; // list of focal sites
			Pairs pairs; // list of pairs sharing focal sites
			
			
			using Map = std::map< size_t, Index >;
		};
		
		
		// construct
		Share();
		
		// create new index for target genotype count
		bool detect(const target_t &, const Grid::Data);
		
		// create index for input list of target sites
		size_t select(const target_t &, const Grid::Data);
		
		// access table
		Index::Map const & get() const;
		
		// print size of index to stream
		void print(const size_t &, std::ostream & = std::cout) const;
		
		// load index map from binary file
		bool load(const std::string &, const Grid::Data);
		
		// save index map to binary file
		void save(const std::string &, const Grid::Data) const;
		
		
		region_t region;
		size_t max_sites;
		size_t max_pairs;
		
		
	private:
		
		// detect share across samples
		void detect_share(const Grid::Data);
		
		// make sharer pairs
		static void create_pairs_index(Index &);
		void create_pairs_table();
		
		// re-generate sites from pairs
		static void remake_sites_index(Index &);
		void remake_sites_table();
		
		// random subset of sites
		static void subset_sites_index(Index &, const size_t &);
		void subset_sites_table();
		
		// random subset of pairs
		static void subset_pairs_index(Index &, const size_t &);
		void subset_pairs_table();
		
		// remove empty indices in table
		void clean_table();
		
		
		// load and save functions
		size_t load_header(Binary &, const Grid::Data);
		Index::Sites::value_type load_index_site(Binary &) const;
		Index::Pairs::value_type load_index_pair(Binary &) const;
		void save_header(Binary &, const Grid::Data) const;
		void save_index_site(Binary &, const Index::Sites::const_iterator &) const;
		void save_index_pair(Binary &, const Index::Pairs::const_iterator &) const;
		
		
		Index::Map table;
		bool good;
		
		// constant char sequence as identifier for binary file navigation
		static const char checkpoint[4];
		
		typedef uint8_t  Bin1;
		typedef uint16_t Bin2;
		typedef uint32_t Bin4;
		typedef uint64_t Bin8;
	};
}


#endif /* GenShare_hpp */
