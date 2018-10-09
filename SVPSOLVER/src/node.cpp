#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <assert.h>
#include <math.h>

#include "node.h"
#include "vector.h"

using namespace std;

#define debug_class	0
#define debug	0

NODE::NODE(){	// default constructor
#if debug_class
	cout << "NODE: default constructor" << endl;
#endif
	m = -1;
	ub = nullptr;
	lb = nullptr;
	relax_objval = - 1.0e+10;
	relax_solval = nullptr;
	dpt = -1;
	zero = true;
	index = -1;
	sumfixed = nullptr;
}

NODE::NODE( const NODE &source )
{	// copy constructor
#if debug_class
	cout << "NODE: copy constructor" << endl;
#endif
	m = source.m;
	relax_objval = source.relax_objval;
	zero = source.zero;
	dpt = source.dpt;
	index = source.index;

	if( m > 0 ){
		assert( source.ub != nullptr );
		assert( source.lb != nullptr );

		ub = new int[m];
		lb = new int[m];

      for( int i = 0; i < m; ++i )
      {
         ub[i] = source.ub[i];
         lb[i] = source.lb[i];
      }

		if( source.relax_solval != nullptr ){
			relax_solval = new double[m];
			Copy_vec( source.relax_solval, relax_solval, m);
		}else{
			relax_solval = nullptr;
		}
		if( source.sumfixed != nullptr ){
			sumfixed = new double[m];
			Copy_vec( source.sumfixed, sumfixed, m);
		}else{
			sumfixed = nullptr;
		}
	}else{
		ub = nullptr;
		lb = nullptr;
		relax_solval = nullptr;
		sumfixed = nullptr;
	}

}

// assignment operator
NODE& NODE::operator=( const NODE& source )
{
#if debug_class
	cout << "NODE: assignment operator" << endl;
#endif

	if( this != &source ){
		m = source.m;
		relax_objval = source.relax_objval;
		zero = source.zero;
		index = source.index;
		dpt = source.dpt;

		if( m > 0 ){
			assert( source.ub != nullptr );
			assert( source.lb != nullptr );

			delete[] ub;
			delete[] lb;
			ub = new int[m];
			lb = new int[m];

         for( int i = 0; i < m; ++i )
         {
            ub[i] = source.ub[i];
            lb[i] = source.lb[i];
         }

			if( source.relax_solval != nullptr ){
				delete[] relax_solval;
				relax_solval = new double[m];
				Copy_vec( source.relax_solval, relax_solval, m);
			}else{
				relax_solval = nullptr;
			}
			if( source.sumfixed != nullptr ){
				delete[] sumfixed;
				sumfixed = new double[m];
				Copy_vec( source.sumfixed, sumfixed, m);
			}else{
				sumfixed = nullptr;
			}
		}else{
			ub = nullptr;
			lb = nullptr;
			relax_solval = nullptr;
			sumfixed = nullptr;
		}
	}

	return *this;
}

NODE::NODE( NODE &&source ) noexcept
{	// move constructor
#if debug_class
	cout << "NODE: move constructor" << endl;
#endif
   m = source.m;
   ub = source.ub;
   lb = source.lb;
   sumfixed = source.sumfixed;
   relax_objval = source.relax_objval;
   relax_solval = source.relax_solval;
   dpt = source.dpt;
   zero = source.zero;
   index = source.index;
	source.ub = nullptr;
	source.lb = nullptr;
	source.relax_solval = nullptr;
	source.sumfixed = nullptr;
}

// destructor
NODE::~NODE()
{
#if debug_class
	cout << "NODE: destructor" << endl;
#endif
	delete[] ub;
	delete[] lb;
	delete[] relax_solval;
	delete[] sumfixed;
	ub = nullptr;
	lb = nullptr;
	relax_solval = nullptr;
	sumfixed = nullptr;
}

void NODE::set_vals(
	const int		s_m,
	const int	   *s_ub,
	const int	   *s_lb,
	const double	*s_relaxsolval,
	const double	s_relax_objval,
	const int		s_dpt,
	const bool 		s_zero,
	const int		s_index
	)
{
	assert( s_m > 0 );
	assert( ub == nullptr );
	assert( lb == nullptr );
	assert( relax_solval == nullptr );

	m = s_m;
	relax_objval = s_relax_objval;
	dpt = s_dpt;
	zero = s_zero;
	index = s_index;

	ub = new int[m];
	lb = new int[m];
	relax_solval = new double[m];

   for( int i = 0; i < m; ++i )
   {
      ub[i] = s_ub[i];
      lb[i] = s_lb[i];
   }
	Copy_vec( s_relaxsolval, relax_solval, m);

}

void NODE::set_relaxsolval(
	double	*solval
	)
{
	assert( solval != nullptr );
	assert( relax_solval != nullptr );
	assert( m > 0 );

	//Copy_vec( solval, relax_solval, m);
	double ep = 1e-10;
	for(int i=0; i<m ; i++){
		if( ceil(solval[i]) - solval[i] < ep ){
			relax_solval[i] = ceil(solval[i]);
		}else if(solval[i] - floor(solval[i]) <ep ){
			relax_solval[i] = floor(solval[i]);
		}else{
			relax_solval[i] = solval[i];
		}
	}
}

bool NODE::alloc_sumfixed()
{
	assert( m > 0 );

	bool result = true;

	if ( sumfixed == nullptr )
   {
		sumfixed = new double[m];
		result = true;
	}
   else
   {
		result = false;
	}

	return result;
}

void NODE::set_sumfixed(
	double	c,
	double	*s_sumfixed
	)
{
	assert( s_sumfixed != nullptr );
	assert( sumfixed != nullptr );
	assert( m > 0 );
	assert( c != 0.0 );

	Com_scal( s_sumfixed, m, c, sumfixed);
}

void NODE::add_sumfixed(
	double	c,
	double	*s_sumfixed
	)
{
	assert( s_sumfixed != nullptr );
	assert( sumfixed != nullptr );
	assert( m > 0 );
	assert( c != 0.0 );

	for(int i=0; i<m; i++){
		sumfixed[i] += c*s_sumfixed[i];
	}
}

void NODE::NODEdispInformation()
{
   cout << "[index:" << index << "]----------------------------" << endl;
   cout << "m: " << m << endl;
   cout << "dpt: " << dpt << endl;
   cout << "zero: " << zero << endl;
   cout << "bounds&relax_solval: " << endl;
   for ( int i = 0; i < m; i++ )
      cout << i << ": [ " << lb[i] << ", " << ub[i] << "], " << relax_solval[i] << endl;
   cout << "relax_objval: " << relax_objval << endl;

}
