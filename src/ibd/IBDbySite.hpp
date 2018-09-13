//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef IBDbySite_hpp
#define IBDbySite_hpp

#include <iomanip>
#include <map>
#include <memory>
#include <vector>

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
	namespace Site
	{
		// Target of detection
		struct Target
		{
			using Data = std::shared_ptr< Target >;
			using List = std::list< Data >;
			
			using Segmentation = std::map< Gen::Sample::Key::Pair, std::pair< Segment, decimal_t > >;
			
			
			// construct
			Target(const size_t &, const Gen::Marker::Key &, const Gen::Sample::Key::Vector &);
			
			
			const size_t fk; // target fk
			const Gen::Marker::Key         focus; // focal site shared by pairs
			const Gen::Sample::Key::Vector share; // individuals sharing site
			
			Segmentation map;
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
			
			// print to stream
			void print(std::ostream & = std::cout) const;
			
			
		private:
			
			// run methods
			void dgt(const Gen::Sample::Key::Pair &, const Gen::Variant::Vector::Data, const Gen::Variant::Vector::Data, const decimal_t &);
			void fgt(const Gen::Sample::Key::Pair &, const Gen::Variant::Vector::Data, const Gen::Variant::Vector::Data, const decimal_t &);
			void hmm(const Gen::Sample::Key::Pair &, const Gen::Variant::Vector::Data, const Gen::Variant::Vector::Data, const decimal_t &);
			
			
			const DetectMethod        method;
			const Target::Data        target; // focal share pair
			const Gen::Grid::Data source; // grid data source
			const HMM::Model::Data model; // HMM model
			const decimal_t max_missing_rate;
		};
	}
}


#endif /* IBDbySite_hpp */

