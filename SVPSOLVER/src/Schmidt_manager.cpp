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

#include "Schmidt_manager.h"
#include "vector.h"

using namespace std;


#define debug_class  0
#define debug  0

SCHMIDT_M::SCHMIDT_M(){ // default constructor
#if debug_class
   cout << "SCHMIDT_M: default constructor" << endl;
#endif
   n = -1;
   min = -1.0;
   z = nullptr;
   B_ = nullptr;
   GS_ = nullptr;
   nrm = nullptr;
}

SCHMIDT_M::SCHMIDT_M( const SCHMIDT_M &source )
{  // copy constructor
#if debug_class
   cout << "SCHMIDT_M: copy constructor" << endl;
#endif
   n = source.n;
   min = source.min;

   if( source.z != nullptr ){
      assert( n > 0 );
      z = new bool[n];
      for(int i=0; i<n; i++) z[i] = source.z[i];
   }else{
      z = nullptr;
   }

   if( source.B_ != nullptr ){
      assert( n > 0 );
      B_ = new double[n*n];
      Copy_vec( source.B_, B_, n*n);
   }else{
      B_ = nullptr;
   }

   if( source.GS_ != nullptr ){
      assert( n > 0 );
      GS_ = new double[n*n];
      Copy_vec( source.GS_, GS_, n*n);
   }else{
      GS_ = nullptr;
   }

   if( source.nrm != nullptr ){
      assert( n > 0 );
      nrm = new double[n];
      Copy_vec( source.nrm, nrm, n);
   }else{
      nrm = nullptr;
   }

}

// assignment operator
SCHMIDT_M& SCHMIDT_M::operator=( const SCHMIDT_M& source )
{
#if debug_class
   cout << "SCHMIDT_M: assignment operator" << endl;
#endif

   if( this != &source ){
      n = source.n;
      min = source.min;

      if( source.z != nullptr ){
         assert( n > 0 );
         z = new bool[n];
         for(int i=0; i<n; i++) z[i] = source.z[i];
      }else{
         z = nullptr;
      }

      if( source.B_ != nullptr ){
         assert( n > 0 );
         B_ = new double[n*n];
         Copy_vec( source.B_, B_, n*n);
      }else{
         B_ = nullptr;
      }

      if( source.GS_ != nullptr ){
         assert( n > 0 );
         GS_ = new double[n*n];
         Copy_vec( source.GS_, GS_, n*n);
      }else{
         GS_ = nullptr;
      }

      if( source.nrm != nullptr ){
         assert( n > 0 );
         nrm = new double[n];
         Copy_vec( source.nrm, nrm, n);
      }else{
         nrm = nullptr;
      }
   }
   return *this;
}

// destructor
SCHMIDT_M::~SCHMIDT_M()
{
#if debug_class
   cout << "SCHMIDT_M: destructor" << endl;
#endif
   delete[] z;
   delete[] B_;
   delete[] GS_;
   delete[] nrm;
   z = nullptr;
   B_ = nullptr;
   GS_ = nullptr;
   nrm = nullptr;
}

void SCHMIDT_M::compute_GS(){

   assert( GS_ != nullptr );
   assert( B_ != nullptr );
   assert( n > 0 );
   assert( z != nullptr );

   min = -1.0;
   for( int i = 0, in = 0; i < n; i++, in += n )
   {
      if( z[i] )
      {
         Copy_vec( &B_[in], &GS_[in], n );

         for ( int j = 0; j < i; j++ )
         {
            if( z[j] )
               Com_linecomb( &GS_[in], &GS_[j*n], n, 1.0, - u(i,j), &GS_[in] );
         }

         nrm[i] = Com_dot( &GS_[in], &GS_[in], n);

         assert( nrm[i] > 0.0 );

         if( min < 0.0 )
            min = nrm[i];
         else if( min > nrm[i] )
            min = nrm[i];
      }
      else
      {
         for ( int j = 0; j < n; j++ )
            GS_[in+j] = 0.0;

         nrm[i] = 0.0;
      }

   }
}

void SCHMIDT_M::setup(
   const int      s_n,
   const double   *s_B_
   )
{
   assert( s_n > 0 );
   assert( s_B_ != nullptr );
   assert( z == nullptr );
   assert( B_ == nullptr );
   assert( GS_ == nullptr );
   assert( nrm == nullptr );

   n = s_n;

   auto nn = n * n;

   B_ = new double[nn];
   GS_ = new double[nn];
   nrm = new double[nn];
   z = new bool[n];

   Copy_vec( s_B_, B_, nn);

   for( int i = 0; i < n; i++ )
   {
      z[i] = true;
      nrm[i] = -1.0;
   }
}

