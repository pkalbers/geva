//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef Identity_h
#define Identity_h

#include <functional>
#include <set>
#include <utility>
#include <vector>


template <typename T>
struct Identity
{
	using Type = T;
	
	// constructs
	Identity() {}
	Identity(const size_t i) : value(i) {}
	
	// assign
	Identity<T> & operator = (const Identity<T> & other)
	{
		this->value = other.value;
		return *this;
	}
	Identity<T> & operator = (const size_t id)
	{
		this->value = id;
		return *this;
	}
	
	// cast
	operator size_t () const
	{
		return this->value;
	}
	
	// sort/compare
	bool operator <  (const Identity<T> & other) { return (this->value <  other.value); }
	bool operator  > (const Identity<T> & other) { return (this->value  > other.value); }
	bool operator <= (const Identity<T> & other) { return (this->value <= other.value); }
	bool operator >= (const Identity<T> & other) { return (this->value >= other.value); }
	bool operator == (const Identity<T> & other) { return (this->value == other.value); }
	bool operator != (const Identity<T> & other) { return (this->value != other.value); }
	
	
	size_t value;
	
	
	// Identity containers
	using Vector = std::vector< Identity<T> >;
	using Set    = std::set< Identity<T> >;
	using Pair   = std::pair< Identity<T>, Identity<T> >;
	
	// hash key
	friend struct std::hash< Identity<T> >;
};


namespace std
{
	template<typename T>
	struct hash< Identity<T> >
	{
		size_t operator () (Identity< T > const & id) const
		{
			return id.value;
		}
	};
}


#endif /* Identity_h */
