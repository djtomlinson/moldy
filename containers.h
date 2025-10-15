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

    //requires function to increase length by 1
    void extend1()
    {
        m_data[0].idCopied(m_length); //prevents new id being created for copies
        T* data {new T[static_cast<std::size_t>(m_length+1)]};
        
        int elementsToCopy {m_length};
        std::copy_n(m_data, elementsToCopy, data);
        
        delete[] m_data;

        m_data = data;
        m_length = m_length+1;
    }

    //Getters
    int getLength() { return m_length; }
};

#endif
