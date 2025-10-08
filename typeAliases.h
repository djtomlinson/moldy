#ifndef TYPEALIASES_H
#define TYPEALIASES_H

//Requires <array> (included in moldy.cpp)
using vec3 = std::array<double,3>;

template <typename T>//, std::size_t size>//, std::size_t Col>
using mat = std::vector<std::vector<T>>;

#endif
