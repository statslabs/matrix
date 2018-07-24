//
// Created by Yi Pan (Institute of Cancer and Genomic Sciences) on 06/03/2018.
//

#include <cstdio>
#include <iostream>
#include "slab/matrix.h"

using namespace std;

int main() {
  slab::vec x1 = {1, 2, 3};
  slab::mat x2 = {
      {1,2,3},
      {4,5,6}
  };
  slab::cube x3 = {
      {{1,2},{3,4}},
      {{5,6},{7,8}},
      {{9,1},{2,3}}
  };

  cout << x1.t() << endl;
  cout << x2.t() << endl;
  //cout << x3.t() << endl;

  return 0;
}