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

	auto m3 = Mat<float>::random(1020, 2000, 2);
	auto m4 = Mat<float>::random(2000, 20, 2);
	auto m1 = m3.select(2, 20, 2, 5);
	auto m2 = m4.select(2, 5, 2, 20);
	//auto m1 = Mat<short>::ones(20);
	//auto m2 = Mat<short>::eye(20);

	// auto m5 = Mat<float>::random(20, 20, 1);
	// auto m6 = Mat<float>::random(20, 20, 1);

	// cout << "----- Some information -----" << endl;
	// cout << "[Alignment] m1: " << m1.is_aligned() << " m2: " << m2.is_aligned() << endl; 

	// cout << "m1=" << endl << m1 << endl;
	// cout << "m2=" << endl << m2 << endl;
	// cout << "m1*m2 = " << endl;
	// cout << m1.mul_bf(m2) << endl;
	// cout << "-----comparsion-----" << endl;
	// cout << m1*m2 << endl;

	// cout << "m5*m6 = " << endl;
	// cout << m5.mul_bf(m6) << endl;
	// cout << "-----comparsion-----" << endl;
	// cout << m5*m6 << endl;

	// TIME_START
	// auto m9 = m3*m4;
	// TIME_END
	//cout << "m3=" << endl << m3 << endl;
	//cout << "m4=" << endl << m4 << endl;
	cout << "m1=" << endl << m1 << endl;
	cout << "m2=" << endl << m2 << endl;
	
	cout << "m1*m2 = " << endl;
	cout << m1.mul_bf(m2) << endl;
	cout << "-----comparsion-----" << endl;
	cout << m1*m2 << endl;
}