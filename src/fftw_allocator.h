#ifndef FFTW_ALLOCATOR_H_
#define FFTW_ALLOCATOR_H_

#include <fftw3.h>

template<typename T>
class fftw_allocator{
public:
	using value_type = T;

	value_type* allocate(std::size_t n){
		return reinterpret_cast<value_type*>(fftw_malloc(n * sizeof(T)));
	}
	void deallocate(value_type *p, std::size_t) noexcept {fftw_free(reinterpret_cast<void*>(p));}
};

template<typename T, typename U>
bool operator==(const fftw_allocator<T>&, const fftw_allocator<U>&) noexcept {return true;}

template<typename T, typename U>
bool operator!=(const fftw_allocator<T>&, const fftw_allocator<U>&) noexcept {return false;}

#endif
