add_executable(basic_matrix_uses 01_basic_matrix_uses.cpp)
target_link_libraries(basic_matrix_uses Matrix::matrix)

add_executable(example example.cpp)
target_link_libraries(example Matrix::matrix)

add_executable(dgemm_example dgemm_example.cpp)
target_link_libraries(dgemm_example Matrix::matrix)