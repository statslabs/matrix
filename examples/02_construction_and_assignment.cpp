#include "slab/matrix.h"
using namespace slab;

#include <iostream>
using namespace std;

enum class Piece {none, cross, naught};

int main() {
  Matrix<Piece, 2> board1 {
      {Piece::none, Piece::none, Piece::none},
      {Piece::none, Piece::none, Piece::none},
      {Piece::none, Piece::none, Piece::cross}
  };

  Matrix<Piece, 2> board2(3,3);  // OK

//  Matrix<Piece, 2> board3{3,3};  // error: constructor from initializer_list<int> deleted

  return 0;
}
