#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <assert.h>
#include <time.h>

#include "stopwatch.h"

using namespace std;

#define debug_class	0
#define debug	0

STOPWATCH::STOPWATCH(){	// default constructor
#if debug_class
	cout << "STOPWATCH: default constructor" << endl;
#endif
	run = false;
	limit = 86400;
   total = 0;
}

STOPWATCH::STOPWATCH( const STOPWATCH &source )
{	// copy constructor
#if debug_class
	cout << "STOPWATCH: copy constructor" << endl;
#endif
	run = source.run;
	limit = source.limit;
   total = source.total;
}

// assignment operator
STOPWATCH& STOPWATCH::operator=( const STOPWATCH& source )
{
#if debug_class
	cout << "STOPWATCH: assignment operator" << endl;
#endif

	if( this != &source )
   {
		run = source.run;
		limit = source.limit;
      total = source.total;
	}

	return *this;
}

// destructor
STOPWATCH::~STOPWATCH()
{
#if debug_class
	cout << "STOPWATCH: destructor" << endl;
#endif
}

void STOPWATCH::start()
{
   assert( run == false );
	begin = time(NULL);
	run = true;
}

void STOPWATCH::stop()
{
   assert( run == true );
	end = time(NULL);
   run = false;
   total += (int)(end - begin);
}

int STOPWATCH::get_time() const
{
   assert( run == true );
	time_t buf = time(NULL);
	return total + (int)(buf - begin);
}

bool STOPWATCH::check_time()
{
   assert( run == true );
	time_t buf = time(NULL);
	bool result = true;
	if ( total + (int)(buf - begin) > limit )
   {
		result = false;
	}else{
		result = true;
	}
	return  result;
}
