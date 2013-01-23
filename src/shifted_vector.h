/*  shifted_vector.h

	Produces a virtual array of size X, but where only the block between start and end form valid values;
	ideal when representing a chromosome but not being able to reserve huge blocks of memory, while simultaneously
	having the 'biological' interface of correct chromosomal positions.

	Eric-Wubbo Lameijer, Leiden University Medical Center, 12-08-2011 */

#ifndef SHIFTED_VECTOR_H_
#define SHIFTED_VECTOR_H_

#include <iostream>
#include <vector>
#include <stdlib.h>

template <class T> 
class ShiftedVector {

public:
	ShiftedVector(unsigned int start, unsigned int end, T fillValue );
	T& operator[](unsigned int position);

private:
	ShiftedVector();
	std::vector<T> m_contents;
	unsigned int m_start;
	unsigned int m_end;
};

template< class T>
ShiftedVector<T>::ShiftedVector(unsigned int start, unsigned int end, T fillValue )
{
	unsigned int size = ( end - start ) + 1;
	m_contents.resize( size, fillValue );
	m_start = start;
	m_end = end;
}

template< class T >
T& ShiftedVector<T>::operator[]( unsigned int position )
{
	if ( (position>m_end) || (position<m_start) ) {
		std::cout << "ShiftedVector out-of-range error: position " << position << " falls outside range " << 
				m_start << "-" << m_end << std::endl;
		exit( EXIT_FAILURE );
	} else {
		return m_contents[ position - m_start ];
	}
}

#endif // SHIFTED_VECTOR_H



