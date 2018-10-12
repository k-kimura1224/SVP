#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <assert.h>

#include "solution.h"
#include "vector.h"

using namespace std;

#define debug_class  0
#define debug  0

SOLUTION::SOLUTION(){   // default constructor
#if debug_class
   cout << "SOLUTION: default constructor" << endl;
#endif
   m = -1;
   val = nullptr;
   objval = -1.0;
}

SOLUTION::SOLUTION( const SOLUTION &source )
{  // copy constructor
#if debug_class
   cout << "SOLUTION: copy constructor" << endl;
#endif
   m = source.m;
   objval = source.objval;

   if( m > 0 ){
      assert( source.val != nullptr );
      val = new double[m];
      Copy_vec( source.val, val, m);
   }else{
      val = nullptr;
   }
}

// assignment operator
SOLUTION& SOLUTION::operator=( const SOLUTION& source )
{
#if debug_class
   cout << "SOLUTION: assignment operator" << endl;
#endif

   if( this != &source ){
      assert( source.val != nullptr );

      m = source.m;
      objval = source.objval;

      if( m > 0 ){
         assert( source.val != nullptr );
         delete[] val;
         val = new double[m];
         Copy_vec( source.val, val, m);
      }else{
         val = nullptr;
      }
   }

   return *this;
}

// destructor
SOLUTION::~SOLUTION()
{
#if debug_class
   cout << "SOLUTION: destructor" << endl;
#endif
   delete[] val;
   val = nullptr;
}

void SOLUTION::set_sol(
   int      s_m,
   double   *s_val,
   double   s_objval
   )
{
   assert( s_m > 0 );
   assert( s_val != nullptr );
   assert( s_objval > 0.0 );
   assert( val == nullptr );

   m = s_m;
   objval = s_objval;

   val = new double[m];

   Copy_vec( s_val, val, m);

}
