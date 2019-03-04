#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <assert.h>

#include "cut.h"
#include "vector.h"

using namespace std;

#define debug_class  0
#define debug  0


//
//void CUT::set_cut(
//   int      s_m,
//   double   *s_coef,
//   double   s_lb,
//   double   s_ub
//   )
//{
//   assert( s_m > 0 );
//   assert( s_coef != nullptr );
//   assert( s_ub > s_lb );
//
//   m = s_m;
//
//   coef = new double[m];
//   lb = new double;
//   ub = new double;
//
//   Copy_vec( s_coef, coef, m);
//   *lb = s_lb;
//   *ub = s_ub;
//
//}
