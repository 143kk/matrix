#pragma once
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cstring>

using namespace std;

#define MAT_ALIGN_BYTES 32
#define OP_ERR(msg) { cerr << "ERROR: " << msg << endl; exit(EXIT_FAILURE); }
#define OP_WARN(msg) { cerr << "WARNING: " << msg << endl; }
#define OP_ASSERT(cond, msg) { if(cond) OP_ERR(msg) }
#define IS_ALIGNED(addr) ((size_t(addr)&(MAT_ALIGN_BYTES-1)) == 0)

struct MatData {
	void * data;
	size_t refcount;
	size_t width;
	size_t _size;

	MatData(size_t size, size_t width): refcount(1), width(width), _size(size){
		data = aligned_alloc(MAT_ALIGN_BYTES, size);
		memset(data, 0, size);
	}

	~MatData(){
		if(this->data) {
			free(data);
			data = NULL;
		}
		refcount = 0;
		width = 0;
	}

	void addRef() {
		refcount++;
	}

	void delRef() {
		refcount--;
		if(!refcount) {
			free(data);
			data = NULL;
		}
	}
};

template <typename T>
class Mat {
	size_t _cols, _rows, _channels;
	//T * data;
	T ** data;
	MatData * md;
	bool aligned;

	T * getAddr(size_t row, size_t col, size_t ch = 0) const{
		return data[ch] + (row-1)*(md->width) + (col-1);
	}

	public:
		Mat(size_t rows, size_t cols, size_t ch = 1);

		Mat(size_t rows, size_t cols, size_t ch, T ** d, MatData * matd, bool aligned);
		
		Mat(const Mat &m);

		~Mat() {
			if(data) release();
		}

		void release() {
			md->delRef();
			delete [] data;
			if(!md->refcount) free(md);
			data = NULL;
			md = NULL;
		}

		#define SAME_TYPE (m._rows == this->_rows && m._cols == this->_cols && m._channels == this->_channels)
		
		Mat& operator=(const Mat &m) {
			if(m == *this) return *this;
			release();
			_rows = m._rows;
			_cols = m._cols;
			_channels = m._channels;
			md = m.md;
			data = m.data;
			aligned = m.aligned;
			md->addRef();
			return *this;
		}

		const Mat operator+(const Mat &m) const {
			OP_ASSERT(!SAME_TYPE, "Matrix types not matched")

			Mat ret(_rows, _cols, _channels);
			
			for(size_t i = 1; i <= m._rows; ++i) {
				for(size_t j = 1; j <= m._cols; ++j) {
					for(size_t k = 0; k < m._channels; ++k) {
						*(ret.getAddr(i, j, k)) = *getAddr(i, j, k) + *(m.getAddr(i, j, k));
					}
				}
			}

			return ret;
		}

		const Mat operator-(const Mat &m) const {
			OP_ASSERT(!SAME_TYPE, "Matrix types not matched")

			Mat ret(_rows, _cols, _channels);
			for(size_t i = 1; i <= m._rows; ++i) {
				for(size_t j = 1; j <= m._cols; ++j) {
					for(size_t k = 0; k < m._channels; ++k) {
						*(ret.getAddr(i, j, k)) = *getAddr(i, j, k) - *(m.getAddr(i, j, k));
					}
				}
			}
			
			return ret;
		}

		void printRef() {
			cout << md->refcount << endl;
		}

		const Mat operator*(const Mat &m) const;

		const Mat mul_bf(const Mat & m) const{
			cout << "brute force!" << endl;
			OP_ASSERT(_cols != m._rows || _channels != m._channels, "Matrix sizes not matched.") 
			Mat ret(_rows, m._cols, _channels);
			for(size_t c = 0; c < _channels; ++c)
				for(size_t i = 1; i <= _rows; ++i) {
					for(size_t k = 1; k <= _cols; ++k) {
						for(size_t j = 1; j <= m._cols; ++j) {
							*(ret.getAddr(i, j, c)) += *(getAddr(i, k, c)) * *(m.getAddr(k, j, c));
						}
					}
				}
			return ret;
		}

		const Mat operator-() const {
			Mat ret(_rows, _cols, _channels);
	
			for(size_t i = 1; i <= _rows; ++i) {
				for(size_t j = 1; j <= _cols; ++j) {
					for(size_t k = 0; k < _channels; ++k) 
						*(ret.getAddr(i, j, k)) = -(*getAddr(i, j, k));
				}
			}

			return ret;
		}

		Mat& operator+= (const Mat &m) {
			OP_ASSERT(!SAME_TYPE, "Matrix types not matched")
			// If data is only used for the current matrix, then no need to request another memory.
			if(md->refcount > 1) {
				*this = m + *this;
			}else{
				for(size_t i = 1; i <= _rows; ++i) {
					for(size_t j = 1; j <= _cols; ++j) {
						for(size_t k = 0; j < _channels; ++k) {
								*getAddr(i, j, k) += *m.getAddr(i, j, k);
						}
					}
				}
			}
			return *this;
		}

		Mat& operator-= (const Mat &m) {
			OP_ASSERT(!SAME_TYPE, "Matrix types not matched")
			// If data is only used for the current matrix, then no need to request another memory.
			if(md->refcount > 1) {
				*this = *this - m;
			}else{
				for(size_t i = 1; i <= _rows; ++i) {
					for(size_t j = 1; j <= _cols; ++j) {
						for(size_t k = 0; j < _channels; ++k)
							*getAddr(i, j, k) -= *m.getAddr(i, j, k);
					}
				}
			}
			return *this;
		}

		// Mat& operator*= (const Mat &m);

		Mat operator*(T x) const {
			Mat ret(_rows, _cols, _channels);
			for(size_t i = 1; i <= _rows; ++i) {
				for(size_t j = 1; j <= _cols; ++j) {
					for(size_t k = 0; k < _channels; ++k) {
						*(ret.getAddr(i, j, k)) = (*getAddr(i, j, k)) * x;
					}
				}
			}

			return ret;
		}	

		Mat& operator*=(T x) {
			// If data is only used for the current matrix, then no need to request another memory.
			if(md->refcount > 1) {
				*this = *this * x;
			}else{
				for(size_t i = 1; i <= _rows; ++i) {
					for(size_t j = 1; j <= _cols; ++j) {
						for(size_t k = 0; k < _channels; ++k)
							(*getAddr(i, j, k)) *= x;
					}
				}
			}
			
			return *this;	
		}

		T& operator()(size_t i, size_t j, size_t c = 0) const{
			OP_ASSERT(i > _rows || j > _cols || i < 1 || j < 1 || c < 0 || c >= _channels, 
				"Argument(s) out of range")
			return *getAddr(i, j, c);
		}
		bool operator==(const Mat &m) const {
			if((_cols != m._cols) || (_rows != m._rows) || (_channels != m._channels)) return false;
			for(size_t i = 0; i < _channels; ++i) {
				if(data[i] != m.data[i]) return false;
			}
			return true;
		}

		template <typename Q>
		friend ostream& operator<<(ostream& os, const Mat<Q>& m);


		Mat select(size_t rbegin, size_t rend, size_t cstart, size_t cend) const;

		// Mat clone() const;
		
		bool sizeEqual(const Mat &m) const { return m._cols == _cols && m._rows == _rows; }

		size_t rows() const { return _rows; }
		size_t cols() const { return _cols; }
		size_t total() const { return _rows * _cols * _channels; }
		bool is_aligned() const { return aligned; }

		static Mat eye(size_t size) {
			Mat ret(size, size);
			memset(ret.md->data, 0, ret.md->_size);
			for(size_t i = 1; i <= size; ++i) {
				ret(i, i) = 1;
			}
			return ret;
		}
		static Mat ones(size_t size);

		static Mat zeros(size_t size) {
			Mat ret(size, size);
			memset(ret.md->data, 0, ret.md->_size);
			return ret;
		}

		static Mat random(size_t rows, size_t cols, size_t channels = 1);

		friend const Mat<float> matmul(const Mat<float> &mat1, const Mat<float> &mat2);

		template<typename Q, typename TT, typename F_FILL, typename F_OP>
		friend void matmul_simd(const Mat<TT> &mat1, const Mat<TT> &mat2, Mat<TT> &res,
			F_FILL simd_fill,
			F_OP simd_op,
			size_t step
		);
};

template <typename T>
Mat<T> operator*(const T x, const Mat<T> &m) {
	return m*x;
}

template <typename Q>
ostream& operator<<(ostream& os, const Mat<Q>& m) {
	if(m._rows < 1 || m._cols < 1) return os;
	if(m._channels == 1) {
		cout << "[";
		for(size_t i = 1; i <= m._rows; ++i) {
			for(size_t j = 1; j < m._cols; ++j) {
				if(i == 1 && j == 1) os << setw(7);
				else os << setw(8);
				os << m(i, j) << ",";
			}
			if(i == m._rows) {
				os << setw(8) << m(i, m._cols) << "]";
			}else{
				os << setw(8) << m(i, m._cols) << ";" << endl;
			}
		}	
	}else{
		for(size_t i = 1; i <= m._rows; ++i) {
			for(size_t j = 1; j <= m._cols; ++j) {
				os << "(";
				for(size_t k = 0; k < m._channels-1; ++k) {
					os << m(i, j, k) << ",";
				}
				os << m(i, j, m._channels-1) << ") ";
			}
			os << endl;
		}
	}

	return os;
}

template <typename T>
Mat<T>::Mat(size_t rows, size_t cols, size_t ch) {
	OP_ASSERT(rows < 1 || cols < 1 || ch < 1, "Invalid matrix size.")
	_rows = rows;
	_cols = cols;
	_channels = ch;

	size_t r_cols = cols;
	
	#ifdef MAT_ALIGN_BYTES
	if(sizeof(T) * cols % MAT_ALIGN_BYTES != 0) {
		size_t e_count = MAT_ALIGN_BYTES / sizeof(T);
		r_cols = (cols / e_count + 1) * e_count;
	}
	#endif

	md = new MatData(sizeof(T) * rows * r_cols * ch, r_cols);

	// init channel
	data = new T* [ch];
	for(size_t i = 0; i < ch; ++i) {
		data[i] = (T*)md->data + (rows * r_cols * i);
	}

	aligned = true;
}

template <typename T>
Mat<T>::Mat(size_t rows, size_t cols, size_t ch, T ** d, MatData * matd, bool aligned):
	_cols(cols),_rows(rows), _channels(ch),
	md(matd),
	aligned(aligned)
{
	data = new T* [ch];
	for(int i = 0; i < ch; ++i) {
		data[i] = d[i];
	}
	matd->addRef();
}

template <typename T>
Mat<T>::Mat(const Mat &m) : _rows(m._rows), _cols(m._cols), md(m.md), aligned(m.aligned), _channels(m._channels) {
	data = new T* [m._channels];
	for(int i = 0; i < m._channels; ++i) {
		data[i] = m.data[i];
	}
	m.md->addRef();
}

template <typename T>
Mat<T> Mat<T>::select(size_t rbegin, size_t rend, size_t cstart, size_t cend) const {
	OP_ASSERT(rbegin < 1 || rend > _rows || cstart < 1 || cend > _cols || rbegin > rend || cstart > cend, "Argument(s) out of range")

	T ** ch_ptrs = new T*[_channels];
	for(size_t i = 0; i < _channels; ++i)
		ch_ptrs[i] = getAddr(rbegin, cstart, i); 
	
	Mat<T> ret = Mat<T>(rend-rbegin+1, cend-cstart+1, _channels, ch_ptrs, this->md, IS_ALIGNED(ch_ptrs[0]));
	delete [] ch_ptrs;
	
	return ret;
}

template <typename T>
Mat<T> Mat<T>::random(size_t rows, size_t cols, size_t channels) {
	Mat<T> ret(rows, cols, channels);
	srand((unsigned int)time(NULL));
	for(int i = 0; i < channels; ++i) {
		for(int j = 1; j <= rows; ++j) {
			for(int k = 1; k <= cols; ++k) {
				*ret.getAddr(j, k, i) = rand()%10 + 1;
			}
		}
	}
	return ret;
}

template <typename T>
Mat<T> Mat<T>::ones(size_t size) {
	Mat<T> ret(size, size);
	size_t total = ret.total();
	for(size_t i = 1; i <= size; ++i)
		for(size_t j = 1; j <= size; ++j)
			*(ret.getAddr(i, j)) = 1;
	return ret;
}

template<typename Q, typename TT, typename F_FILL, typename F_OP>
void matmul_simd(const Mat<TT> &mat1, const Mat<TT> &mat2, Mat<TT> &res,
	F_FILL simd_fill,
	F_OP simd_op,
	const size_t step
) {
	register size_t mat1_rows = mat1._rows;
	register size_t mat1_cols = mat1._cols;
	register size_t mat2_cols = mat2._cols;

	Q mat1_vec;
	TT * resPtr, * mat2Ptr, * mat1Ptr, * resPtr_init;
	for(size_t c = 0; c < mat1._channels; ++c)
	#pragma omp parallel for schedule(static) \
	private(mat1_vec, resPtr, mat1Ptr, mat2Ptr, resPtr_init) \
	if(mat1_rows > 64)
		for(size_t i = 1; i <= mat1_rows; i++) {
			mat1Ptr = mat1.getAddr(i, 1, c);
			resPtr_init = res.getAddr(i, 1, c);
			for(size_t k = 1; k <= mat1_cols; ++k) {
					mat2Ptr = mat2.getAddr(k, 1, c);
					resPtr = resPtr_init;
					simd_fill(mat1_vec, mat1Ptr);
					size_t j = 0;
					for(; j+step < mat2_cols; j += step) {
							simd_op(resPtr, mat1_vec, mat2Ptr);
							resPtr += step;
							mat2Ptr += step;
					}
					for(; j < mat2_cols; ++j) {
						*resPtr += (*mat1Ptr) * (*mat2Ptr);
						resPtr++;
						mat2Ptr++;
					}

					mat1Ptr++;
			}
		}
}

template <typename T>
const Mat<T> Mat<T>::operator*(const Mat<T> &m) const {
	return this->mul_bf(m);
}

template <>
const Mat<double> Mat<double>::operator*(const Mat<double> &m) const;

template <>
const Mat<float> Mat<float>::operator*(const Mat<float> &m) const;

template <>
const Mat<int> Mat<int>::operator*(const Mat<int> &m) const;
