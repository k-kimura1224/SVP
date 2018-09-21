#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <assert.h>

#include "qpdata.h"
#include "vector.h"

using namespace std;

#define debug_class  0
#define debug  0

QP_DATA::QP_DATA(){  // default constructor
#if debug_class
   cout << "QP_DATA: default constructor" << endl;
#endif
   n = -1;
   Qmat = nullptr;
   lvec = nullptr;
   uvec = nullptr;
   wvec = nullptr;
   pvec = nullptr;
}

QP_DATA::QP_DATA( const QP_DATA &source )
{  // copy constructor
#if debug_class
   cout << "QP_DATA: copy constructor" << endl;
#endif
   n = source.n;

   if( n > 0 )
   {
      assert( source.Qmat != nullptr );
      assert( source.lvec != nullptr );
      assert( source.uvec != nullptr );
      assert( source.wvec != nullptr );
      assert( source.pvec != nullptr );

      auto nn = n * n;

      Qmat = new double[nn];
      lvec = new double[n];
      uvec = new double[n];
      wvec = new double[n];
      pvec = new double[n];

      Copy_vec( source.Qmat, Qmat, nn );
      Copy_vec( source.lvec, lvec, n );
      Copy_vec( source.uvec, uvec, n );
      Copy_vec( source.wvec, wvec, n );
      Copy_vec( source.pvec, pvec, n );
   }
   else
   {
      Qmat = nullptr;
      lvec = nullptr;
      uvec = nullptr;
      wvec = nullptr;
      pvec = nullptr;
   }
}

// assignment operator
QP_DATA& QP_DATA::operator=( const QP_DATA& source )
{
#if debug_class
   cout << "QP_DATA: assignment operator" << endl;
#endif

   if( this != &source )
   {
      n = source.n;

      if( n > 0 )
      {
         assert( source.Qmat != nullptr );
         assert( source.lvec != nullptr );
         assert( source.uvec != nullptr );
         assert( source.wvec != nullptr );
         assert( source.pvec != nullptr );

         delete[] Qmat;
         delete[] lvec;
         delete[] uvec;
         delete[] wvec;
         delete[] pvec;

         auto nn = n * n;

         Qmat = new double[nn];
         lvec = new double[n];
         uvec = new double[n];
         wvec = new double[n];
         pvec = new double[n];

         Copy_vec( source.Qmat, Qmat, nn );
         Copy_vec( source.lvec, lvec, n );
         Copy_vec( source.uvec, uvec, n );
         Copy_vec( source.wvec, wvec, n );
         Copy_vec( source.pvec, pvec, n );
      }
      else
      {
         Qmat = nullptr;
         lvec = nullptr;
         uvec = nullptr;
         wvec = nullptr;
         pvec = nullptr;
      }
   }

   return *this;
}

// destructor
QP_DATA::~QP_DATA()
{
#if debug_class
   cout << "QP_DATA: destructor" << endl;
#endif

   delete[] Qmat;
   delete[] lvec;
   delete[] uvec;
   delete[] wvec;
   delete[] pvec;

   Qmat = nullptr;
   lvec = nullptr;
   uvec = nullptr;
   wvec = nullptr;
   pvec = nullptr;
}

//void QP_DATA::set_data(
//   int      s_m,
//   double   *s_B,
//   double   *s_B_,
//   double   *s_Q
//   )
//{
//   assert( m == -1 );
//   assert( B == nullptr );
//   assert( B_ == nullptr );
//   assert( Q == nullptr );
//
//   assert( s_m > 0 );
//   assert( s_B != nullptr );
//   assert( s_B_ != nullptr );
//   assert( s_Q != nullptr );
//
//   m = s_m;
//
//   int mm = m * m;
//
//   B = new double[mm];
//   B_ = new double[mm];
//   Q = new double[mm];
//
//   Copy_vec( s_B, B, mm );
//   Copy_vec( s_B_, B_, mm );
//   Copy_vec( s_Q, Q, mm );
//}
//
// allocate memories and return result
void QP_DATA::alloc(
      int   s_n
      )
{
   assert( check_allocation() == false );

   n = s_n;

   Qmat = new double[n * n];
   lvec = new double[n];
   uvec = new double[n];
   wvec = new double[n];
   pvec = new double[n];
}

// check whether data is defined and return this result
bool QP_DATA::check_allocation()
{
   auto result = false;
   if( n > 0 )
   {
      assert( Qmat != nullptr );
      assert( lvec != nullptr );
      assert( uvec != nullptr );
      assert( wvec != nullptr );
      assert( pvec != nullptr );

      result = true;
   }
   else
   {
      assert( Qmat == nullptr );
      assert( lvec == nullptr );
      assert( uvec == nullptr );
      assert( wvec == nullptr );
      assert( pvec == nullptr );

      result = false;
   }

   return result;
}

