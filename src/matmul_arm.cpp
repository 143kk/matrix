#include <arm_neon.h>
#include "omp.h"


// float
inline void mat_load_mat1_f32(float32x4_t &mat1_vec, float * mat1Ptr) {
	mat1_vec = vdupq_n_f32(*mat1Ptr);
}

inline void mat_comb_f32(float * resPtr, float32x4_t &mat1_vec, float * mat2Ptr) {
	vst1q_f32(resPtr, vmlaq_f32(vld1q_f32(resPtr), mat1_vec, vld1q_f32(mat2Ptr)));
}

// double
inline void mat_load_mat1_f64(float64x2_t &mat1_vec, double * mat1Ptr) {
	mat1_vec = vdupq_n_f64(*mat1Ptr);
}

inline void mat_comb_f64(double * resPtr, float64x2_t &mat1_vec, double * mat2Ptr) {
	vst1q_f64(resPtr, vmlaq_f64(vld1q_f64(resPtr), mat1_vec, vld1q_f64(mat2Ptr)));
}

//int
inline void mat_load_mat1_s32(int32x4_t &mat1_vec, int * mat1Ptr) {
	mat1_vec = vdupq_n_s32(*mat1Ptr);
}

inline void mat_mm256_comb_s32(int * resPtr, int32x4_t &mat1_vec, int * mat2Ptr) {
		vst1q_s32(resPtr, vmlaq_s32(vld1q_s32(resPtr), mat1_vec, vld1q_s32(mat2Ptr)));
}



template <>
const Mat<float> Mat<float>::operator*(const Mat<float> &m) const {
	OP_ASSERT(_cols != m._rows || _channels != m._channels, "Matrix sizes not matched.") 
	Mat<float> res(_rows, m._cols, _channels);	
	matmul_simd<float32x4_t>(*this, m, res, mat_load_mat1_f32, mat_comb_f32, 4);
	return res;
}

template <>
const Mat<double> Mat<double>::operator*(const Mat<double> &m) const {
	OP_ASSERT(_cols != m._rows || _channels != m._channels, "Matrix sizes not matched.") 
	Mat<double> res(_rows, m._cols, _channels);	
	matmul_simd<float64x2_t>(*this, m, res, mat_load_mat1_f64, mat_comb_f64, 2);
	return res;
}

template <>
const Mat<int> Mat<int>::operator*(const Mat<int> &m) const {
	OP_ASSERT(_cols != m._rows || _channels != m._channels, "Matrix sizes not matched.") 
	Mat<int> res(_rows, m._cols, _channels);	
	matmul_simd<int32x4_t>(*this, m, res, mat_load_mat1_s32, mat_comb_s32, 4);
	return res;
}