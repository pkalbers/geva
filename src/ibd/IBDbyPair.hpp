//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef IBDbyPair_hpp
#define IBDbyPair_hpp

#include <iomanip>
#include <algorithm>
#include <list>
#include <map>
#include <memory>
#include <vector>

#include "Decimal.h"
#include "Random.h"
#include "Progress.hpp"

#include "IBD.hpp"
#include "IBD_DGT.hpp"
#include "IBD_FGT.hpp"
#include "IBD_HMM.hpp"

#include "Gen.hpp"
#include "GenSample.hpp"
#include "GenMarker.hpp"
#include "GenVariant.hpp"
#include "GenGrid.hpp"
#include "GenShare.hpp"


namespace IBD
{
	namespace Pair
	{
		// Target of detection
		class Target
		{
		public:
			
			using Data = std::shared_ptr< Target >;
			using List = std::list< Data >;
			
			using Segmentation = std::map< Segment, Gen::Marker::Key::Set >;
			
			
			// constructs
			Target(const size_t &, const Gen::Sample::Key::Pair &, const Gen::Marker::Key::Set &);
			
			
			const size_t fk; // target fk
			const Gen::Sample::Key::Pair   pair; // pair of individuals
			const Gen::Marker::Key::Vector site; // focal sites shared by pair
			
			Segmentation map;
			
			decimal_t missing; // calculated pairwise missing rate
			
			
		private:
			
			// shuffle focal sites
			static Gen::Marker::Key::Vector shuffle(const Gen::Marker::Key::Set &);
		};
		
		
		// Result container
		class Result
		{
		public:
			
			using Data = std::shared_ptr< Result >;
			
			
			// construct
			Result(const Gen::Share::Data);
			
			// return number of execution units
			size_t size() const;
			
			// access list
			Target::List const & get() const;
			
			// print to stream
			size_t print(std::ostream & = std::cout) const;
			
			
		private:
			
			Target::List target;
			size_t       length;
		};
		
		
		// Exection of detection
		class Detect
		{
		public:
			
			// construct
			Detect(const DetectMethod, const decimal_t, const Target::Data, const Gen::Grid::Data, const HMM::Model::Data = nullptr);
			
			// execute detection using selected methods
			void run();
			
			
		private:
			
			// run methods
			void dgt(const Gen::Variant::Vector::Data, const Gen::Variant::Vector::Data);
			void fgt(const Gen::Variant::Vector::Data, const Gen::Variant::Vector::Data);
			void hmm(const Gen::Variant::Vector::Data, const Gen::Variant::Vector::Data);
			
			
			const DetectMethod    method;
			const Target::Data    target; // focal share pair
			const Gen::Grid::Data source; // grid data source
			const HMM::Model::Data model; // HMM model
			const decimal_t max_missing_rate;
		};
	}
}


#endif /* IBDbyPair_hpp */

