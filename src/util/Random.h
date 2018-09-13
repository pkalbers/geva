//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef Random_h
#define Random_h

#include <algorithm>
#include <chrono>
#include <limits>
#include <random>
#include <set>
#include <string>
#include <vector>


using random_engine_f = std::mt19937;
using random_engine_t = random_engine_f::result_type;

using random_output_t = size_t;
using random_vector_t = std::vector< random_output_t >;


//  random seed
static random_engine_t random_seed = static_cast<random_engine_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count() % std::numeric_limits<random_engine_t>::max());


// random generators
static random_engine_f random_generator(random_seed);


// get and set random seed

inline void set_random_seed(const size_t seed)
{
	random_seed = static_cast<random_engine_t>(seed);
	random_generator = random_engine_f(random_seed);
}

inline random_engine_t get_random_seed()
{
	return random_seed;
}


// random number

inline random_output_t random_number()
{
	return static_cast<random_output_t>(random_generator());
}

inline random_output_t random_number(const size_t upper)
{
	return std::uniform_int_distribution<random_output_t>(0, upper - 1)(random_generator);
}

inline random_output_t random_number(const size_t lower, const size_t upper)
{
	return std::uniform_int_distribution<random_output_t>(lower, upper - 1)(random_generator);
}


// random number vector

inline random_vector_t random_vector(const size_t size, const bool unique = false, const bool sorted = false)
{
	random_vector_t v(size);
	
	if (unique)
	{
		std::set< random_output_t > u;
		
		if (sorted)
		{
			while (u.size() < size)
				u.insert(random_number());
			
			v.assign(u.begin(), u.end());
		}
		else
		{
			for (size_t i = 0; i < size; ++i)
			{
				random_output_t n = random_number();
				u.insert(n);
				
				while (u.size() < i + 1)
				{
					n = random_number();
					u.insert(n);
				}
				
				v[i] = n;
			}
		}
	}
	else
	{
		for (size_t i = 0; i < size; ++i)
			v[i] = random_number();
		
		if (sorted)
		{
			std::sort(v.begin(), v.end());
		}
	}
	
	return v;
}

inline random_vector_t random_vector(const size_t size, const size_t upper, const bool unique = false, const bool sorted = false)
{
	random_vector_t v(size);
	
	if (unique)
	{
		if (upper < size)
			throw std::logic_error("Random number vector size must be larger than limit");
		
		std::set< random_output_t > u;
		
		if (sorted)
		{
			while (u.size() < size)
				u.insert(random_number(upper));
			
			v.assign(u.begin(), u.end());
		}
		else
		{
			for (size_t i = 0; i < size; ++i)
			{
				random_output_t n = random_number(upper);
				u.insert(n);
				
				while (u.size() < i + 1)
				{
					n = random_number(upper);
					u.insert(n);
				}
				
				v[i] = n;
			}
		}
	}
	else
	{
		for (size_t i = 0; i < size; ++i)
			v[i] = random_number(upper);
		
		if (sorted)
		{
			std::sort(v.begin(), v.end());
		}
	}
	
	return v;
}

inline random_vector_t random_vector(const size_t size, const size_t lower, const size_t upper, const bool unique = false, const bool sorted = false)
{
	random_vector_t v(size);
	
	if (unique)
	{
		if (upper - lower < size)
			throw std::logic_error("Random number vector size must be larger than range");
		
		std::set< random_output_t > u;
		
		if (sorted)
		{
			while (u.size() < size)
				u.insert(random_number(lower, upper));
			
			v.assign(u.begin(), u.end());
		}
		else
		{
			for (size_t i = 0; i < size; ++i)
			{
				random_output_t n = random_number(lower, upper);
				u.insert(n);
				
				while (u.size() < i + 1)
				{
					n = random_number(lower, upper);
					u.insert(n);
				}
				
				v[i] = n;
			}
		}
	}
	else
	{
		for (size_t i = 0; i < size; ++i)
			v[i] = random_number(lower, upper);
		
		if (sorted)
		{
			std::sort(v.begin(), v.end());
		}
	}
	
	return v;
}


// coin toss

inline bool random_coin()
{
	std::bernoulli_distribution distr;
	return distr(random_generator);
	
	// static constexpr random_output_t half = static_cast<random_output_t>(((std::numeric_limits<random_engine_t>::max() - 1) / 2) + 1);
	// return (random_number() < half);
}


// random string

inline std::string random_string(const size_t n)
{
	static const char charset[] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
	static const size_t charlen = sizeof(charset) - 1;
	
	std::string str(n, '\0');
	
	for (size_t i = 0; i < n; ++i)
	{
		str[i] = charset[ random_number(charlen) ];
	}
	
	return str;
}


#endif /* Random_h */

