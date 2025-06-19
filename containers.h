#ifndef CONTAINERS_H
#define CONTAINERS_H

#include <cassert>
#include <cstddef>


template <typename T>
class Array
{
private:
    int m_length {};
    T* m_data {};
public:
    //Constructor
    Array(int length)
    {
        assert(length > 0);
        m_data = new T[length] {};
        m_length = length;
    }
    //Destructor
    ~Array()
    {
        delete[] m_data;
    }
    //Indexing array
    T& operator[](int index)
    {
        assert(index >= 0 && index < m_length);
        return m_data[index];
    }

    //Getters
    int getLength() { return m_length; }
};

#endif
