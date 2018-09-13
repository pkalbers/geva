//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#include "GenVariant.hpp"


using namespace Gen;


// Variant


// construct

Variant::Variant(const gen_t genotype)
: g(genotype)
, h(genotype_to_haplotypes(genotype))
, p(genotype_is_phased(genotype))
{}

Variant::Variant(const Variant & other)
: g(other.g)
, h(other.h)
, p(other.p)
{}


// return genotype

gen_t Variant::gen() const
{
	return this->g;
}


// return haplotype

hap_t Variant::hap(const ChrType chr) const
{
	if (!this->p)
	{
		throw std::invalid_argument("Genotype is not phased");
	}
	return this->h[ chr ];
}


// check if phased

bool Variant::is_phased() const
{
	return this->p;
}



// Vector


// construct

Variant::Vector::Vector()
: g(0)
, h{{ hap_vector_t(0), hap_vector_t(0) }}
, p(false)
, n(0)
, good(false)
{}

Variant::Vector::Vector(const gen_vector_t & gv)
{
	*this = gv;
}

Variant::Vector::Vector(gen_vector_t && gv)
{
	*this = std::move(gv);
}

Variant::Vector::Vector(const gen_vector_t & gv, const bool phased)
: p(phased)
, good(false)
{
	*this = gv;
}

Variant::Vector::Vector(gen_vector_t && gv, const bool phased)
: g(std::move(gv))
, h{{ hap_vector_t(0), hap_vector_t(0) }}
, p(phased)
, good(false)
{
	this->n = this->g.size();
	
	if (this->n == 0)
	{
		throw std::runtime_error("Empty variant vector was supplied");
	}
}

Variant::Vector::Vector(const Vector & other)
: g(other.g)
, h(other.h)
, p(other.p)
, n(other.n)
, good(other.good)
{}

Variant::Vector::Vector(Vector && other)
: g(std::move(other.g))
, h(std::move(other.h))
, p(other.p)
, n(other.n)
, good(other.good)
{}


// assign

Variant::Vector & Variant::Vector::operator = (const gen_vector_t & gv)
{
	gen_vector_t tmp = gv;
	
	return (*this = std::move(tmp));
}

Variant::Vector & Variant::Vector::operator = (gen_vector_t && gv)
{
	this->guard.lock();
	
	this->g.swap(gv);
	
	this->n = this->g.size();
	
	if (this->n == 0)
	{
		throw std::runtime_error("Variant vector is empty");
	}
	
	bool phased = true;
	
	for (size_t i = 0; i < this->n; ++i)
	{
		if (phased && !genotype_is_phased(this->g.at(i)))
		{
			phased = false;
			break;
		}
	}
	
	this->p = phased;
	
	this->guard.unlock();
	
	return *this;
}


// return variant

Variant Variant::Vector::at(const Marker::Key & i) const
{
	if (this->n == 0)
	{
		throw std::invalid_argument("Invalid variant vector");
	}
	return Variant(this->g.at(i));
}

Variant Variant::Vector::operator [] (const Marker::Key & i) const
{
	return this->at(i);
}


// return genotype

gen_vector_t const & Variant::Vector::gen() const
{
	if (this->n == 0)
	{
		throw std::invalid_argument("Invalid variant vector");
	}
	return this->g;
}

gen_t Variant::Vector::gen(const Marker::Key & i) const
{
	return this->gen().at(i);
}


// return haplotype

hap_vector_t const & Variant::Vector::hap(const ChrType chr)
{
	if (this->n == 0)
	{
		throw std::invalid_argument("Invalid variant vector");
	}
	
	if (!this->p)
	{
		throw std::invalid_argument("Variant vector is not phased");
	}
	
	if (!this->good)
	{
		this->guard.lock();
		
		hap_vector_t & mat = this->h[ MATERNAL ];
		hap_vector_t & pat = this->h[ PATERNAL ];
		
		mat.resize(this->n);
		pat.resize(this->n);
		
		for (size_t i = 0; i < this->n; ++i)
		{
			const hap_pair_t hh = genotype_to_haplotypes(this->g[i]);
			
			mat[i] = hh[ MATERNAL ];
			pat[i] = hh[ PATERNAL ];
		}
		
		this->good = true;
		
		this->guard.unlock();
	}
	
	return this->h[ chr ];
}

hap_t Variant::Vector::hap(const ChrType chr, const Marker::Key & i) const
{
	if (this->good)
	{
		return this->h[chr].at(i);
	}
	
	if (this->n == 0)
	{
		throw std::invalid_argument("Invalid variant vector");
	}
	
	if (!this->p)
	{
		throw std::invalid_argument("Variant vector is not phased");
	}
	
	const hap_pair_t hh = genotype_to_haplotypes(this->g[i]);
	
	return hh[ chr ];
}


// check if phased

bool Variant::Vector::is_phased() const
{
	return this->p;
}


// enforce phasing

void Variant::Vector::is_phased(const bool phased)
{
	this->guard.lock();
	this->p = phased;
	this->guard.unlock();
}


// return size

size_t Variant::Vector::size() const
{
	return this->n;
}

