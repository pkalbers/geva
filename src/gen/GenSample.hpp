//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef GenSample_hpp
#define GenSample_hpp

#include "Identity.h"

#include "Gen.hpp"


namespace Gen
{
	struct Sample
	{
		using Key = Identity<Sample>;
		
		typedef std::vector< Sample >   Vector;
		typedef Vector::const_iterator  Iterator;
		typedef Vector::const_reference Reference;
		
		
		// constructs
		Sample();
		//Sample(const Sample &);
		//Sample(Sample &&);
		
		
		// print to stream
		void print(std::ostream & = std::cout) const;
		
		// return as string
		std::string str() const;
		
		// compare
		bool operator == (const Sample &) const;
		bool operator != (const Sample &) const;
		
		
		Key         index;
		std::string label;
		bool        phase;
		
		static const std::string header; // header for printing
	};
}


#endif /* GenSample_hpp */

