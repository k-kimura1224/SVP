#ifndef SOLUTION_POOL_H__
#define SOLUTION_POOL_H__

#include "solution.h"

class SOLUTION_POOL{
   int         size;
   int         num;
   SOLUTION    *sols;
   public:

      SOLUTION_POOL();                                // default constructor
      SOLUTION_POOL( const SOLUTION_POOL &source );         // copy constructor
      SOLUTION_POOL& operator=( const SOLUTION_POOL& );     // assignment operator
      ~SOLUTION_POOL();                               // destructor

      void  alloc( int s );
      void  add_solution( SOLUTION sol );
      auto  get_biggestval() {
         assert( num > 0 && size > 0 );
         return sols[num-1].get_objval();
      }
      auto get_num() const { return num; }

      void  disp() const {
         assert( sols != nullptr );
         assert( size > 0 && num > 0 );
         for ( auto i = 0; i < num; ++i )
         {
            printf("(%d) objval:%f solval:", i, sols[i].get_objval() );
            auto solvals = sols[i].get_solval();
            auto m = sols[i].get_dim();
            for ( auto j = 0; j < m; ++j )
               printf("%d ", (int) solvals[j] );
            printf("\n");
         }
      }

      auto& cget_sols() const { return sols; }
};
#endif
