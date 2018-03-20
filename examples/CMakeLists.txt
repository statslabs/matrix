add_executable(example example.cpp)
target_link_libraries(example Matrix::matrix)

add_executable(dgemm_example dgemm_example.cpp)
target_link_libraries(dgemm_example Matrix::matrix)