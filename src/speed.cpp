#include <iostream>
#include <chrono>
#include "matrix.hpp"

using namespace std::chrono;
#define TIME_START t1 = steady_clock::now();
#define TIME_END t2 = steady_clock::now();\
	cout << duration_cast<milliseconds>(t2-t1).count() << "ms.\n";\


int main() {
	//Mat<float> m2 = Mat<float>::eye(99);
	//Mat<float> m1 = Mat<float>::ones(2544);
	//auto m3 = matmul(m1, m1);
	//cout << m3 << endl;

	auto t1 = steady_clock::now();
	auto t2 = steady_clock::now();
	double elasped;	

	auto m1 = Mat<float>::random(2048, 2048);
	auto m2 = Mat<float>::random(2048, 2048);


	TIME_START
	auto m9 = m1*m2;
	TIME_END

	TIME_START
	auto m10 = m1.mul_bf(m2);
	TIME_END
	// // //cout << "m3=" << endl << m3 << endl;
	// // //cout << "m4=" << endl << m4 << endl;
	// // cout << "m1=" << endl << m1 << endl;
	// // cout << "m2=" << endl << m2 << endl;
	
	// // cout << "m1*m2 = " << endl;
	// // cout << m1.mul_bf(m2) << endl;
	// // cout << "-----comparsion-----" << endl;
	// // cout << m1*m2 << endl;
}