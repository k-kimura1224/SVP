#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <assert.h>

#include "probdata.h"
#include "vector.h"

using namespace std;

#define debug_class	0
#define debug	0

PROB_DATA::PROB_DATA(){	// default constructor
#if debug_class
	cout << "PROB_DATA: default constructor" << endl;
#endif
	m = -1;
	B = nullptr;
	B_ = nullptr;
	Q = nullptr;
}

PROB_DATA::PROB_DATA( const PROB_DATA &source )
{	// copy constructor
#if debug_class
	cout << "PROB_DATA: copy constructor" << endl;
#endif
	m	=	source.m;

	if( m > 0 )
   {
		assert( source.B != nullptr );
		assert( source.B_ != nullptr );
		assert( source.Q != nullptr );

      int mm = m * m;

		B = new double[mm];
		B_ = new double[mm];
		Q = new double[mm];

		Copy_vec( source.B, B, mm );
		Copy_vec( source.B_, B_, mm );
		Copy_vec( source.Q, Q, mm );
	}
   else
   {
		B = nullptr;
		B_ = nullptr;
		Q = nullptr;
	}
}

// assignment operator
PROB_DATA& PROB_DATA::operator=( const PROB_DATA& source )
{
#if debug_class
	cout << "PROB_DATA: assignment operator" << endl;
#endif

	if( this != &source )
   {
		m = source.m;

		if( m > 0 )
      {
			assert( source.B != nullptr );
			assert( source.B_ != nullptr );
			assert( source.Q != nullptr );

			delete[] B;
			delete[] B_;
			delete[] Q;

         int mm = m * m;

			B = new double[mm];
			B_ = new double[mm];
			Q = new double[mm];

			Copy_vec( source.B, B, mm );
			Copy_vec( source.B_, B_, mm );
			Copy_vec( source.Q, Q, mm );
		}
      else
      {
			B = nullptr;
			B_ = nullptr;
			Q = nullptr;
		}
	}

	return *this;
}

// destructor
PROB_DATA::~PROB_DATA()
{
#if debug_class
	cout << "PROB_DATA: destructor" << endl;
#endif
	delete[] B;
	delete[] B_;
	delete[] Q;
	B = nullptr;
	B_ = nullptr;
	Q = nullptr;
}

void PROB_DATA::set_data(
	const int		s_m,
	const double	*s_B,
	const double	*s_B_,
	const double	*s_Q
	)
{
	assert( m == -1 );
   assert( B == nullptr );
   assert( B_ == nullptr );
   assert( Q == nullptr );

	assert( s_m > 0 );
	assert( s_B != nullptr );
	assert( s_B_ != nullptr );
	assert( s_Q != nullptr );

	m = s_m;

   int mm = m * m;

	B = new double[mm];
	B_ = new double[mm];
	Q = new double[mm];

	Copy_vec( s_B, B, mm );
	Copy_vec( s_B_, B_, mm );
	Copy_vec( s_Q, Q, mm );
}

