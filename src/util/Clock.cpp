//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#include "Clock.hpp"


// add clock

void Clock::operator << (const Clock & other)
{
	this->process << other.process;
}


// print formated string for both clocks

void Clock::print(std::ostream & stream) const
{
	stream << "Process time: " << this->process.str() << std::endl;
	stream << "Elapsed time: " << this->elapsed.str() << std::endl;
}


// reset both clocks

void Clock::reset()
{
	this->process.reset();
	this->elapsed.reset();
}


// Timestamp of current date/time

std::string Clock::stamp(bool special_chars)
{
	static const char format_t[] = "%Y-%m-%d %H:%M:%S";
	static const char format_f[] = "%Y-%m-%d.%Hh%Mm%Ss";
	
	time_t     current = time(0);
	struct tm  tstruct = *localtime(&current);
	char       buffer[256];
	
	if (special_chars)
		strftime(buffer, sizeof(buffer), format_t, &tstruct);
	else
		strftime(buffer, sizeof(buffer), format_f, &tstruct);
	
	return buffer;
}
