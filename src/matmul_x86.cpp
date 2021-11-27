#include "omp.h"
#include <immintrin.h>
#include "matrix.hpp"

// float
inline void mat_mm256_load_mat1_ps(__m256 &mat1_vec, float * mat1Ptr) {
	mat1_vec = _mm256_broadcast_ss(mat1Ptr);
}

inline void mat_mm256_comb_ps(float * resPtr, __m256 &mat1_vec, float * mat2Ptr) {
	_mm256_store_ps(resPtr, _mm256_fmadd_ps(mat1_vec, _mm256_load_ps(mat2Ptr), _mm256_load_ps(resPtr)));
}

inline void mat_mm256_combu_ps(float * resPtr, __m256 &mat1_vec, float * mat2Ptr) {
	_mm256_store_ps(resPtr, _mm256_fmadd_ps(mat1_vec, _mm256_loadu_ps(mat2Ptr), _mm256_load_ps(resPtr)));
}

// double
inline void mat_mm256_load_mat1_pd(__m256d &mat1_vec, double * mat1Ptr) {
	mat1_vec = _mm256_broadcast_sd(mat1Ptr);
}

inline void mat_mm256_comb_pd(double * resPtr, __m256d &mat1_vec, double * mat2Ptr) {
	_mm256_store_pd(resPtr, _mm256_fmadd_pd(mat1_vec, _mm256_load_pd(mat2Ptr), _mm256_load_pd(resPtr)));
}

inline void mat_mm256_combu_pd(double * resPtr, __m256d &mat1_vec, double * mat2Ptr) {
	_mm256_store_pd(resPtr, _mm256_fmadd_pd(mat1_vec, _mm256_loadu_pd(mat2Ptr), _mm256_load_pd(resPtr)));
}

//int
inline void mat_mm256_load_mat1_epi32(__m256i &mat1_vec, int * mat1Ptr) {
	mat1_vec = _mm256_set1_epi32(*mat1Ptr);
}

inline void mat_mm256_comb_epi32(int * resPtr, __m256i &mat1_vec, int * mat2Ptr) {
	_mm256_store_si256((__m256i *)resPtr, 
		_mm256_add_epi32(
			_mm256_load_si256((__m256i *)resPtr), 
			_mm256_mullo_epi32(mat1_vec, _mm256_load_si256((__m256i *)mat2Ptr) )
		)
	);
}

inline void mat_mm256_combu_epi32(int * resPtr, __m256i &mat1_vec, int * mat2Ptr) {
	_mm256_store_si256((__m256i *)resPtr, 
		_mm256_add_epi32(
			_mm256_load_si256((__m256i *)resPtr), 
			_mm256_mullo_epi32(mat1_vec, _mm256_loadu_si256((__m256i_u *)mat2Ptr) )
		)
	);
}

template <>
const Mat<float> Mat<float>::operator*(const Mat<float> &m) const {
	OP_ASSERT(_cols != m._rows || _channels != m._channels, "Matrix sizes not matched.") 
	Mat<float> res(_rows, m._cols, _channels);	
	if(m.aligned) {
		matmul_simd<__m256>(*this, m, res, mat_mm256_load_mat1_ps, mat_mm256_comb_ps, 8);
	}else{
		matmul_simd<__m256>(*this, m, res, mat_mm256_load_mat1_ps, mat_mm256_combu_ps, 8);
	}
	return res;
}

template <>
const Mat<double> Mat<double>::operator*(const Mat<double> &m) const {
	OP_ASSERT(_cols != m._rows || _channels != m._channels, "Matrix sizes not matched.") 
	Mat<double> res(_rows, m._cols, _channels);	
	if(m.aligned) {
		matmul_simd<__m256d>(*this, m, res, mat_mm256_load_mat1_pd, mat_mm256_comb_pd, 4);
	}else{
		matmul_simd<__m256d>(*this, m, res, mat_mm256_load_mat1_pd, mat_mm256_combu_pd, 4);
	}
	return res;
}

template <>
const Mat<int> Mat<int>::operator*(const Mat<int> &m) const {
	OP_ASSERT(_cols != m._rows || _channels != m._channels, "Matrix sizes not matched.") 
	Mat<int> res(_rows, m._cols, _channels);	
	if(m.aligned) {
		matmul_simd<__m256i>(*this, m, res, mat_mm256_load_mat1_epi32, mat_mm256_comb_epi32, 8);
	}else{
		matmul_simd<__m256i>(*this, m, res, mat_mm256_load_mat1_epi32, mat_mm256_combu_epi32, 8);
	}
	return res;
}