#include <iostream>
using namespace std;

#include "slab/matrix.h"
using namespace slab;

enum class Piece { none, cross, naught };

void ConstructFromExtents() {
  cout << "ConstructFromExtents():" << endl;

  Matrix<double, 0> m0;
  Matrix<double, 1> m1(3);
  Matrix<double, 2> m2(3, 4);
  Matrix<double, 3> m3(3, 4, 5);

  cout << "\nm0.desc: " << m0.descriptor() << "\nm0: " << m0
       << "\nm1.desc: " << m1.descriptor() << "\nm1: " << m1
       << "\nm2.desc: " << m2.descriptor() << "\nm2: " << m2
       << "\nm3.desc: " << m3.descriptor() << "\nm3: " << m3 << endl;
}

void ConstructFromMatrixInitializer() {
  cout << "ConstructFromMatrixInitializer():" << endl;

  Matrix<double, 0> m0{1};
  Matrix<double, 1> m1{1, 2};
  Matrix<double, 2> m2{{1, 2}, {3, 4}};
  Matrix<double, 3> m3{{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}};

  cout << "\nm0.desc: " << m0.descriptor() << "\nm0: " << m0
       << "\nm1.desc: " << m1.descriptor() << "\nm1: " << m1
       << "\nm2.desc: " << m2.descriptor() << "\nm2: " << m2
       << "\nm3.desc: " << m3.descriptor() << "\nm3: " << m3 << endl;

  Matrix<Piece, 2> board1{{Piece::none, Piece::none, Piece::none},
                          {Piece::none, Piece::none, Piece::none},
                          {Piece::none, Piece::none, Piece::cross}};

  Matrix<Piece, 2> board2(3, 3);  // OK

  //  Matrix<Piece, 2> board3{3,3};  // error: constructor from
  //  initializer_list<int> deleted
}

void ConstructFromMatrixRef() {
  cout << "ConstructFromMatrixRef():" << endl;

  Matrix<double, 1> m1{1, 2, 3};
  auto mr1 = m1(slice(1));
  Matrix<double, 1> m1sub(mr1);

  Matrix<double, 2> m2{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  auto mr2 = m2(slice(1), slice(1));
  Matrix<double, 2> m2sub(mr2);

  Matrix<double, 3> m3{{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}},
                       {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}},
                       {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}};
  auto mr3 = m3(slice(1), slice(1), slice(1));
  Matrix<double, 3> m3sub(mr3);

  cout << "\nm1sub: " << m1sub.descriptor() << "\nm1sub: " << m1sub
       << "\nm2sub.desc: " << m2sub.descriptor() << "\nm2sub: " << m2sub
       << "\nm3sub.desc: " << m3sub.descriptor() << "\nm3sub: " << m3sub
       << endl;
}

void AssignFromMatrixInitializer() {
  cout << "AssignFromMatrixInitializer():" << endl;

  Matrix<double, 0> m0 = {1};
  Matrix<double, 1> m1 = {1, 2};
  Matrix<double, 2> m2 = {{1, 2}, {3, 4}};
  Matrix<double, 3> m3 = {{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}};

  cout << "\nm0.desc: " << m0.descriptor() << "\nm0: " << m0
       << "\nm1.desc: " << m1.descriptor() << "\nm1: " << m1
       << "\nm2.desc: " << m2.descriptor() << "\nm2: " << m2
       << "\nm3.desc: " << m3.descriptor() << "\nm3: " << m3 << endl;

  Matrix<Piece, 2> board1 = {{Piece::none, Piece::none, Piece::none},
                             {Piece::none, Piece::none, Piece::none},
                             {Piece::none, Piece::none, Piece::cross}};

  Matrix<Piece, 2> board2(3, 3);  // OK

  //  Matrix<Piece, 2> board3 = {3,3};  // error: constructor from
  //  initializer_list<int> deleted
}

void AssignFromMatrixRef() {
  cout << "AssignFromMatrixRef():" << endl;

  Matrix<double, 1> m1{1, 2, 3};
  auto mr1 = m1(slice(1));
  Matrix<double, 1> m1sub = mr1;

  Matrix<double, 2> m2{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  auto mr2 = m2(slice(1), slice(1));
  Matrix<double, 2> m2sub = mr2;

  Matrix<double, 3> m3{{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}},
                       {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}},
                       {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}};
  auto mr3 = m3(slice(1), slice(1), slice(1));
  Matrix<double, 3> m3sub = mr3;

  cout << "\nm1sub: " << m1sub.descriptor() << "\nm1sub: " << m1sub
       << "\nm2sub.desc: " << m2sub.descriptor() << "\nm2sub: " << m2sub
       << "\nm3sub.desc: " << m3sub.descriptor() << "\nm3sub: " << m3sub
       << endl;
}

int main() {
  ConstructFromExtents();
  ConstructFromMatrixInitializer();
  ConstructFromMatrixRef();

  AssignFromMatrixInitializer();
  AssignFromMatrixRef();

  return 0;
}
