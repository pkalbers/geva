//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#include "GenSample.hpp"


using namespace Gen;


// constructs

Sample::Sample()
: index(0)
, phase(false)
{}

/*
Sample::Sample(const Sample & other)
: index(other.index)
, label(other.label)
, phase(other.phase)
{}

Sample::Sample(Sample && other)
: index(other.index)
, label(std::move(other.label))
, phase(other.phase)
{}
*/


// print to stream

void Sample::print(std::ostream & stream) const
{
	stream << this->index.value << ' ';
	stream << this->label << ' ';
	stream << ((this->phase) ? 1: 0);
}


// return as string

std::string Sample::str() const
{
	return this->label;
}


// compare

bool Sample::operator == (const Sample & other) const
{
	return (this->index == other.index && this->label == other.label && this->phase == other.phase);
}
bool Sample::operator != (const Sample & other) const
{
	return !(*this == other);
}


// header for printing

const std::string Sample::header = "SampleID Label Phased";

