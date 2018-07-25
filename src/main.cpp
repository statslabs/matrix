//
// Created by Yi Pan (Institute of Cancer and Genomic Sciences) on 06/03/2018.
//

#include <cstdio>
#include <iostream>
#include "slab/matrix.h"

using namespace std;

int main() {
  slab::mat x1 = {
      {1,2,3},
      {4,5,6}
  };

  slab::mat x2 = {
      {1,2,3},
      {4,5,6},
      {7,8,9}
  };

  slab::mat x3 = {
      {10,10},
      {10,10}
  };

  cout << slab::join_cols(x1, x2) << endl;
  cout << slab::join_rows(x1, x3) << endl;

  cout << x2(slab::slice{1}, slab::slice{1}) << endl;
  cout << slab::vectorise(x2(slab::slice{1}, slab::slice{1})) << endl;


  return 0;
}