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

#define debug_class  0
#define debug  0

NODE::NODE(){  // default constructor
#if debug_class
   cout << "NODE: default constructor" << endl;
#endif
   m = -1;
   relax_objval = - 1.0e+10;
   relax_solval = nullptr;
   dpt = -1;
   index = -1;
   sumfixed = nullptr;
   type = -1;
}

NODE::NODE( const NODE &source )
{  // copy constructor
#if debug_class
   cout << "NODE: copy constructor" << endl;
#endif
   assert( source.m > 0 );
   assert( source.dpt >= 0 );
   assert( source.index >= 0 );
   assert( source.type >= 0 );
   assert( source.relax_objval >= 0.0 );
   assert( source.relax_solval != nullptr );

   m = source.m;
   dpt = source.dpt;
   index = source.index;
   type = source.type;
   relax_objval = source.relax_objval;
   branchinfo = source.branchinfo;

   relax_solval = new double[m];
   Copy_vec( source.relax_solval, relax_solval, m);

   if( source.sumfixed != nullptr )
   {
      sumfixed = new double[m];
      Copy_vec( source.sumfixed, sumfixed, m);
   }
   else
   {
      sumfixed = nullptr;
   }
}

// assignment operator
NODE& NODE::operator=( const NODE& source )
{
#if debug_class
   cout << "NODE: assignment operator" << endl;
#endif

   if( this != &source )
   {
      assert( source.m > 0 );
      assert( source.dpt >= 0 );
      assert( source.index >= 0 );
      assert( source.type >= 0 );
      assert( source.relax_objval >= 0.0 );
      assert( source.relax_solval != nullptr );

      m = source.m;
      dpt = source.dpt;
      index = source.index;
      type = source.type;
      relax_objval = source.relax_objval;
      branchinfo = source.branchinfo;

      delete[] relax_solval;
      relax_solval = new double[m];
      Copy_vec( source.relax_solval, relax_solval, m);

      if( source.sumfixed != nullptr )
      {
         delete[] sumfixed;
         sumfixed = new double[m];
         Copy_vec( source.sumfixed, sumfixed, m);
      }
      else
      {
         sumfixed = nullptr;
      }
   }

   return *this;
}

NODE::NODE( NODE &&source ) noexcept
{  // move constructor
#if debug_class
   cout << "NODE: move constructor" << endl;
#endif
   m = source.m;
   sumfixed = source.sumfixed;
   relax_objval = source.relax_objval;
   relax_solval = source.relax_solval;
   dpt = source.dpt;
   index = source.index;
   branchinfo = source.branchinfo;
   type = source.type;

   source.relax_solval = nullptr;
   source.sumfixed = nullptr;
   source.branchinfo.clear();
   source.branchinfo.shrink_to_fit();
   source.type = -1;
}

NODE& NODE::operator=( NODE&& source ) noexcept
{
#if debug_class
   cout << "NODE: move assignment operator" << endl;
#endif
   if ( this != &source )
   {
      m = source.m;
      dpt = source.dpt;
      index = source.index;
      type = source.type;
      sumfixed = source.sumfixed;
      relax_objval = source.relax_objval;
      relax_solval = source.relax_solval;
      branchinfo = source.branchinfo;

      source.relax_solval = nullptr;
      source.sumfixed = nullptr;
      source.branchinfo.clear();
      source.branchinfo.shrink_to_fit();
      source.type = -1;
   }
   return *this;
}
// destructor
NODE::~NODE()
{
#if debug_class
   cout << "NODE: destructor" << endl;
#endif
   delete[] relax_solval;
   delete[] sumfixed;
   relax_solval = nullptr;
   sumfixed = nullptr;
   branchinfo.clear();
   branchinfo.shrink_to_fit();
}

void NODE::set_vals(
   const int      s_m,
   const double   *s_relaxsolval,
   const double   s_relax_objval,
   const int      s_dpt,
   const int      s_index,
   const int      s_type
   )
{
   assert( s_m > 0 );
   assert( s_relaxsolval != nullptr );
   assert( s_relax_objval >= 0.0 );
   assert( s_dpt >= 0 );
   assert( s_index >= 0 );
   assert( s_type >= 0 );

   m = s_m;
   relax_objval = s_relax_objval;
   dpt = s_dpt;
   index = s_index;
   type = s_type;

   if ( relax_solval == nullptr )
   {
      relax_solval = new double[m];
   }

   Copy_vec( s_relaxsolval, relax_solval, m);
}

void NODE::set_relaxsolval(
   double   *solval
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
   const int      c,
   const double*  s_sumfixed
   )
{
   assert( s_sumfixed != nullptr );
   assert( sumfixed != nullptr );
   assert( m > 0 );
   assert( c != 0.0 );

   for (int i = 0; i < m; ++i )
      sumfixed[i] = c * s_sumfixed[i];
}

void NODE::add_sumfixed(
   const int      c,
   const double*  s_sumfixed
   )
{
   assert( s_sumfixed != nullptr );
   assert( sumfixed != nullptr );
   assert( m > 0 );
   assert( c != 0.0 );

   for (int i = 0; i < m; ++i )
      sumfixed[i] += c * s_sumfixed[i];
}

void NODE::NODEdispInformation()
{
   cout << "[index:" << index << "]----------------------------" << endl;
   cout << "m: " << m << endl;
   cout << "dpt: " << dpt << endl;
   cout << "type: " << type << endl;
   cout << "relax_objval: " << relax_objval << endl;
}
