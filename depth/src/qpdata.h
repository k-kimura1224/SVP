#ifndef QPDATA_H__
#define QPDATA_H__

using namespace std;

class QP_DATA{
   int         n;          // dimention
   double*     Qmat;       // [n*n] Colmajor
   double*     lvec;       // [n]
   double*     uvec;       // [n]
   double*     wvec;       // [n]
   double*     pvec;       // [n]

   public:

      QP_DATA();                                // default constructor
      QP_DATA( const QP_DATA &source );         // copy constructor
      QP_DATA& operator=( const QP_DATA& );     // assignment operator
      ~QP_DATA();                               // destructor

      //void set_data( int s_m, double *s_B, double *s_B_, double *s_Q );
      bool check_allocation();
      void alloc( int s_n );

      // get-function
      int      get_n(){ return n; }
      double*  get_Qmat(){ return Qmat; }
      double*  get_lvec(){ return lvec; }
      double*  get_uvec(){ return uvec; }
      double*  get_wvec(){ return wvec; }
      double*  get_pvec(){ return pvec; }

};

#endif
