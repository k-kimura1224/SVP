#ifndef CUT_PLANE_H__
#define CUT_PLANE_H__

#include <vector>

using namespace std;

class CUT{
   vector<double>    coef;   //
   double*           lb;     // [1]
   double*           ub;     // [1]

   // lb <= c_1 x_1 + ... c_n x_n <= ub

   public:

      // default constructor
      CUT() { lb = nullptr; ub = nullptr; }

      // constructor
      CUT( const int size, const double* setcoef, const double* setlb, const double* setub ) {
         assert( size > 0 ); assert( setcoef != nullptr );
         coef.resize( size );
         for ( auto i = 0; i < size; ++i ) coef[i] = setcoef[i];
         if ( setlb != nullptr ) {
            lb = new double; *lb = *setlb;
         } else lb = nullptr;
         if ( setub != nullptr ) {
            ub = new double; *ub = *setub;
         } else ub = nullptr;
      }

      // copy constructor
      CUT( const CUT &source ) {
         coef = source.coef;
         if ( source.lb != nullptr ) *lb = *(source.lb);
         else lb = nullptr;
         if ( source.ub != nullptr ) *ub = *(source.ub);
         else ub = nullptr;
      }

      // assignment operator
      CUT& operator=( const CUT& source ) {
         if ( this != &source ) {
            coef = source.coef;
            if ( source.lb != nullptr ) *lb = *(source.lb);
            else lb = nullptr;
            if ( source.ub != nullptr ) *ub = *(source.ub);
            else ub = nullptr;
         }
         return *this;
      }

      // move constructor
      CUT( CUT &&source ) noexcept {
         coef = move( source.coef );
         lb = source.lb;
         ub = source.ub;

         source.coef.clear();
         source.coef.shrink_to_fit();
         source.lb = nullptr;
         source.ub = nullptr;
      }

      // move assignment operator
      CUT& operator=( CUT&& source ) noexcept {
         if ( this != &source ) {
            coef = move( source.coef );
            lb = source.lb;
            ub = source.ub;

            source.coef.clear();
            source.coef.shrink_to_fit();
            source.lb = nullptr;
            source.ub = nullptr;
         }
         return *this;
      }

      // destructor
      ~CUT() { coef.clear(); coef.shrink_to_fit(); lb = nullptr; ub = nullptr; }


//      void set_cut(int s_m, double *s_coef, double s_lb, double s_ub);
//
      const vector<double>&  get_coef() const { return coef; }
      const double  get_lb() const { return *lb; }
      const double  get_ub() const { return *ub; }
};

#endif
