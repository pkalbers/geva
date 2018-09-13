//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef Approx_h
#define Approx_h

#include <limits>


// Approximation through linear interpolattion

template< typename Float >
Float approx(const Float x,
			 const Float x0, const Float x1,
			 const Float y0, const Float y1,
			 const Float lower_bound = std::numeric_limits<Float>::lowest(),
			 const Float upper_bound = std::numeric_limits<Float>::max())
{
	static constexpr Float epsn = std::numeric_limits<Float>::epsilon();
	static constexpr Float half = static_cast<Float>(0.5);
	
	const Float num = x  - x0;
	const Float den = x1 - x0;
	
	// prevent division by zero
	if (den < epsn)
	{
		return y0 + ((y1 - y0) * half);
	}
	
	const Float y = y0 + ((y1 - y0) * (num / den));
	
	// lower bound
	if (y < lower_bound)
	{
		return lower_bound;
	}
	
	// upper bound
	if (y > upper_bound)
	{
		return upper_bound;
	}
	
	return y;
}


#endif /* Approx_h */

