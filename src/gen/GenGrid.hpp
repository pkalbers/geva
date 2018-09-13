//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef GenGrid_hpp
#define GenGrid_hpp


#include <algorithm>
#include <fstream>
#include <iterator>
#include <list>
#include <memory>
#include <mutex>
#include <unordered_map>
#include <sstream>

#include "Random.h"

#include "Binary.hpp"

#include "Gen.hpp"
#include "GenSample.hpp"
#include "GenMarker.hpp"
#include "GenVariant.hpp"


namespace Gen
{
	// Data matrix grid
	class Grid
	{
	public:
		
		using Data = std::shared_ptr< Grid >;
		
		using Type   = value_t; // genotype
		using Vector = value_vector_t; // vector of genotypes
		
		
		// constructs
		Grid(const std::string &);
		Grid(Grid &&); // move
		Grid(const Grid &) = delete; // no copy
		
		
		// read from file
		Vector read(const Sample::Key &, const bool = true);
		
		// fetch from cache
		Variant::Vector::Data get(const Sample::Key &);
		
		// limit cache size
		void cache(const size_t = 0);
		
		
		// get sample/marker information
		Sample::Reference sample(const Sample::Key &) const;
		Marker::Reference marker(const Marker::Key &) const;
		
		// get direct reference to sample/marker vector
		Sample::Vector const & sample() const;
		Marker::Vector const & marker() const;
		
		// return sample/marker size
		size_t sample_size() const { return this->size_sample; }
		size_t marker_size() const { return this->size_marker; }
		
		// return compression
		bool compressed() const { return this->compression; }
		
		// print samples/markers
		void print_sample(std::ostream & = std::cout) const;
		void print_marker(std::ostream & = std::cout) const;
		void print_sample(const std::string &) const;
		void print_marker(const std::string &) const;
		
		
	private:
		
		using interval_t = std::array< size_t, 2 >;
		using guard_t    = std::lock_guard<std::mutex>;
		using buffer_t   = std::unordered_map< size_t, Variant::Vector::Data >;
		
		
		// construct
		Grid(Binary &&);

		// create file index
		void load();
		
		// read sample/marker information
		void load_sample();
		void load_marker();
		
		// prune buffer
		void prune();
		
		// variables read from header
		size_t   size_sample;
		size_t   size_marker;
		interval_t  interval;
		bool     compression;
		
		Binary   source; // binary source file
		buffer_t buffer; // variant data of each individual
		
		size_t buffer_limit; // max cache size
		size_t buffer_count; // current cache count
		
		Sample::Vector sample_list; // vector of sample information
		Marker::Vector marker_list; // vector of marker information
		
		std::mutex guard; // mutex lock
		
		// constant char sequence as identifier for binary file navigation
		static const char checkpoint[4];
		
		
	public:
		
		
		// Make new grid
		class Make
		{
		private:
			
			using vector_t = std::vector< Type >;
			using matrix_t = std::vector< Vector >;
			using source_t = std::vector< Binary >;
			
			const std::string unique; // unique identifier
			const std::string output; // output file
			const bool        comprs; // compression setting
			matrix_t          matrix; // genotype matrix
			size_t        irow, nrow; // number of samples in matrix
			size_t        icol, ncol; // number of markers in vector
			
			interval_t interval; // marker interval in buffer
			
			source_t sources; // vector of binary temp. files
			
			
		public:
			
			// construct
			Make(const std::string &, const size_t &, const size_t &, const bool);
			
			// fill buffer with genotypes
			void insert(const gen_t);
			
			// write buffer to binary file and reset buffer
			void save(const bool = false);
			
			// save overhang and concatenate source files
			void finish(Sample::Vector &, Marker::Vector &);
			
			bool good;
			bool full;
		};
		
		
		// Join multiple grids
		class Join
		{
		private:
			
			// Grid source
			struct Source
			{
				// construct
				Source(const Grid::Data);
				
				// sort, check overlap
				bool operator <  (const Source &) const;
				bool operator  > (const Source &) const;
				bool operator == (const Source &) const;
				
				
				const Grid::Data grid;
				
				const size_t position_min;
				const size_t position_max;
			};
			
			using source_t = std::list< Source >;
			using vector_t = std::vector< Type >;
			
			
			const std::string output; // output file
			source_t          source; // list of grids
			
			bool           compression; // common compression setting
			size_t         marker_size; // total marker size
			size_t         sample_size; // common sample size
			Sample::Vector sample_list; // vector of sample information
			
			bool good;
			
			
		public:
			
			// construct
			Join(const std::string &);
			
			// insert grid
			void insert(const Grid::Data);
			
			// join inserted grids
			void finish();
		};
		
		
		friend class Make;
		friend class Join;
		
		
	private:
		
		typedef uint8_t  Bin1;
		typedef uint16_t Bin2;
		typedef uint32_t Bin4;
		typedef uint64_t Bin8;
		
		
		// save sample/marker information
		static void save_sample(Binary &, const Sample &);
		static void save_marker(Binary &, const Marker &);
	};
}


#endif /* GenGrid_hpp */
