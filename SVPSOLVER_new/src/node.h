#ifndef NODE_H__
#define NODE_H__

#include <vector>

using namespace std;

//enum BRANCH_TYPE {
//   UPPER,   // x_i <= d
//   LOWER,   // d <= x_i
//   EQUAL    // x_i = d
//};

class BRANCH_INFO {
   int     index;
   int     value;
   char    type;

   public:
      // default constructor
      // constructor
      // copy constructor
      // assignment operator
      // move constructor
      // move assignment operator
      // destructor
      BRANCH_INFO() { index = -1; value = 0; type = 'n'; }
      BRANCH_INFO( const int i, const int v, const char t )
      { index = i; value = v; type = t; }
      BRANCH_INFO( const BRANCH_INFO &source )
      { index = source.index; value = source.value; type = source.type; }
      BRANCH_INFO& operator=( const BRANCH_INFO& source )
      {  if ( this != &source )
         { index = source.index; value = source.value; type = source.type; }
         return *this; }
      ~BRANCH_INFO() {}

      int   get_index() const { return index; }
      int   get_value() const { return value; }
      char  get_type() const { return type; }
};


class NODE {
   int            dpt;
   int            index;
   double         relax_objval;
   vector<int>    fixedvalues; // descending order
   vector<double> rsol;

   // after presolving
   //bool           presolved;
   //vector<int> ub;
   //vector<int> lb;
   //vector<int>    sum_fixed;

   //double*     relax_solval;
   //double*  sumfixed;
   //int      m;
   //int      type;
   //vector<BRANCH_INFO>  branchinfo;

   public:

      NODE();                                      // default constructor
      NODE(
      //const bool              s_presolve,
      const int               s_dpt,
      const int               s_index,
      const int               s_relax_objval,
      const vector<int>&      s_fixedavalues
      );
      NODE( const NODE &source );                  // copy constructor
      NODE& operator=( const NODE& );              // assignment operator
      NODE( NODE &&source ) noexcept;              // move constructor
      NODE& operator=( NODE&& source ) noexcept;   // move assignment operator
      ~NODE();                                     // destructor
      bool operator<(const NODE &rhs) const
      { return relax_objval < rhs.relax_objval; }

      //void  set_vals(   const int s_m, const double *s_relaxsolval,
      //                  const double s_relax_objval, const int s_dpt,
      //                  const int s_index, const int s_type );
      //void set_presolved( bool done ) { presolved = done; }
      //void set_ub_i( const int i, const int val )
      //{ assert( !ub.empty() ); assert( i >= 0 && i < (int) ub.size() ); ub[i] = val; }
      //void set_lb_i( const int i, const int val )
      //{ assert( !lb.empty() ); assert( i >= 0 && i < (int) lb.size() ); lb[i] = val; }
      //void set_sumfixed( const vector<int>& s_sum )
      //{ assert( sum_fixed.empty() ); sum_fixed = s_sum; }
      //void add_sumfixed( const int index, const int coef, const int b_i )
      //{ assert( !sum_fixed.empty() ); sum_fixed[index] += coef * b_i; }
      //void alloc_sumfixed( const int size )
      //{ assert( sum_fixed.empty() ); sum_fixed.resize( size ); }

      double   get_lowerbound() const { return relax_objval; }
      //int      get_dimention() const { return m; }
      int      get_index() const { return index; }
      int      get_dpt() const { return dpt; }
      //int      get_type() const { return type; }
      //double*  get_relaxsolval() const { return relax_solval; }
      //double*  get_sumfixed() const { return sumfixed; }
      //const vector<BRANCH_INFO>& get_branchinfo() const { return branchinfo; }
      auto& cget_fixedvalues() const { return fixedvalues; }
      auto& get_fixedvalues()  { return fixedvalues; }
      //auto get_presolved() const { return presolved; }

      //auto& get_sum_fixed() const { return sum_fixed; }
      //auto& geti_sum_fixed() { return sum_fixed; }
      //auto& geti_ub() { return ub; }
      //auto& get_ub() const { return ub; }
      //auto& geti_lb() { return lb; }
      //auto& get_lb() const { return lb; }
      auto& geti_rsol() { return rsol; }
      auto& get_rsol() const { return rsol; }

      void     NODEdispInformation();

      void  set_lowerbound( const double s_ro ) { relax_objval = s_ro; }
      //void  set_relaxsolval( double *solval );
      //bool  alloc_sumfixed();
      //void  set_sumfixed( const int c, const double* s_sumfixed );
      //void  add_sumfixed( const int c, const double* s_sumfixed );

      void  set_dpt( const int s_dpt )
      {
         assert( dpt >= 0 );
         dpt = s_dpt;
      }
      void  set_index( const int s_index )
      {
         assert( index >= 0 );
         index = s_index;
      }
      //void  set_branchinfo( const int index, const int value, const char type ) {
      //   branchinfo.emplace_back( index, value, type );
      //}
      void  inc_dpt() { ++dpt; }
      void  push_back_fixedvalue( const int f ) { fixedvalues.push_back( f ); }
};

#endif
