#ifndef PREPROCESS_H
#define PREPROCESS_H

#include <vector>
#include "solution_pool.h"

using namespace std;

void FindNShorterVec(
      const int                     numofsol,
      const int                     timelimit,
      const vector<vector<double>>& matrix,
      SOLUTION_POOL*                solpool
      );

#endif
