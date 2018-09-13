//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef GenVariant_hpp
#define GenVariant_hpp

#include <memory>
#include <mutex>
#include <stdexcept>

#include "Gen.hpp"
#include "GenMarker.hpp"


namespace Gen
{
	// Variant data
	class Variant
	{
	public:
		
		// Chromosome data
		class Vector
		{
		public:
			
			using Data = std::shared_ptr< Vector >;
			
			
			// construct
			Vector();
			Vector(const gen_vector_t &);
			Vector(gen_vector_t &&);
			Vector(const gen_vector_t &, const bool);
			Vector(gen_vector_t &&, const bool);
			Vector(const Vector &);
			Vector(Vector &&);
			
			// assign
			Vector & operator = (const gen_vector_t &);
			Vector & operator = (gen_vector_t &&);
			
			// return variant
			Variant at(const Marker::Key &) const;
			Variant operator [] (const Marker::Key &) const;
			
			// return genotype
			gen_vector_t const & gen() const;
			gen_t gen(const Marker::Key &) const;
			
			// return haplotype
			hap_vector_t const & hap(const ChrType);
			hap_t hap(const ChrType, const Marker::Key &) const;
			
			// check if phased
			bool is_phased() const;
			
			// enforce phasing
			void is_phased(const bool);
			
			// return size
			size_t size() const;
			
			
		private:
			
			using hap_paired_t = std::array< hap_vector_t, ploidy >;
			
			
			gen_vector_t g;
			hap_paired_t h;
			bool         p;
			size_t       n;
			
			bool good;
			
			std::mutex guard;
		};
		
		
		// construct
		Variant(const gen_t);
		Variant(const Variant &);
		
		// return genotype
		gen_t gen() const;
		
		// return haplotype
		hap_t hap(const ChrType) const;
		
		// check if phased
		bool is_phased() const;
		
		
	private:
		
		const gen_t      g; // genotype
		const hap_pair_t h; // haplotypes
		const bool       p; // flag if genotype is phased
	};
}


#endif /* GenVariant_hpp */

