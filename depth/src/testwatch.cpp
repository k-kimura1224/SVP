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

#include "testwatch.h"

using namespace std;

#define debug_class	0
#define debug	0

TESTWATCH::TESTWATCH(){	// default constructor
#if debug_class
	cout << "TESTWATCH: default constructor" << endl;
#endif
	run = false;
   total = 0.0;
}

TESTWATCH::TESTWATCH( const TESTWATCH &source )
{	// copy constructor
#if debug_class
	cout << "TESTWATCH: copy constructor" << endl;
#endif
	run = source.run;
   total = source.total;
}

// assignment operator
TESTWATCH& TESTWATCH::operator=( const TESTWATCH& source )
{
#if debug_class
	cout << "TESTWATCH: assignment operator" << endl;
#endif

	if( this != &source )
   {
		run = source.run;
      total = source.total;
	}

	return *this;
}

// destructor
TESTWATCH::~TESTWATCH()
{
#if debug_class
	cout << "TESTWATCH: destructor" << endl;
#endif
}

void TESTWATCH::start()
{
   assert( run == false );
	begin = clock();
	run = true;
}

void TESTWATCH::stop()
{
   assert( run == true );
	end = clock();
   run = false;
   total += (double) (end - begin) / CLOCKS_PER_SEC;
}

double TESTWATCH::get_time() const
{
   assert( run == true );
	clock_t buf = clock();
	return total + (double) (buf - begin) / CLOCKS_PER_SEC;
}

