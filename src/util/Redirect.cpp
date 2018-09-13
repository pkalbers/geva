//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#include "Redirect.hpp"


// constructs

Redirect::Redirect(std::ostream & stream, const std::string & filename)
: source(stream)
, buffer(stream.rdbuf())
, output(filename)
{
	this->source.rdbuf(this->output.rdbuf());
}

Redirect::Redirect(std::ostream & stream, std::streambuf * ptr)
: source(stream)
, buffer(stream.rdbuf())
{
	if (ptr)
		this->source.rdbuf(ptr);
	else
	{
		std::cout << "zzzzzz44342232z" << std::endl;
		this->source.setstate(std::ios_base::failbit);
		std::cout << "zzzzzzz" << std::endl;
	}
}


// destruct

Redirect::~Redirect()
{
	this->source.rdbuf(this->buffer);
	
	if (this->output.is_open())
		this->output.close();
}
