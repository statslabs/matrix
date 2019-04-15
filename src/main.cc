#include "slab/matrix.h"

using namespace std;

int main() {
  slab::mat A = {{1, 2, 3}, {0, 1, 4}, {5, 6, 0}};

  cout << "inverse of A is " << solve(A) << endl;

  slab::mat B;
  slab::inv(B, A);

  cout << "inverse of A is " << B << endl;

  return 0;
}
