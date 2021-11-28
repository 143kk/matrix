#include "matrix.hpp"
#include <iostream>

using namespace std;

int main() {
	srand((unsigned int)time(NULL));
	auto m1 = Mat<int>::ones(5);
	for(int i = 1; i <= 5; ++i) {
		m1(i, i) = 100000;
	}
	cout << "m1=" << endl;
	cout << m1 << endl;

	Mat<int> m2 = m1.select(2, 3, 2, 5);
	cout << "m2=" << endl;
	cout << m2 << endl;
	
	

	cout << "-m2=" << endl << (-m2) << endl;

	auto m3 = m1.select(1, 2, 2, 3);
	cout << "m3=" << endl << m3 << endl;

	auto ones = Mat<int>::ones(2);
	cout << "m3+ones=" << endl << (m3 + ones) << endl;
	cout << "m3-ones=" << endl << (m3 - ones) << endl;

	cout << "m2*9=" << endl << (m2*9) << endl;
	cout << "9*m2=" << endl << (9*m2) << endl;

	auto f1 = Mat<float>::ones(4);
	for(int i = 1; i <= 4; ++i) {
		f1(i, i) = 10.f;
	} 
	cout << "f1=" << endl << f1 << endl;

	cout << "f1+=ones" << endl << (f1 += Mat<float>::ones(4)) << endl;
	cout << "f1-=ones" << endl << (f1 -= Mat<float>::ones(4)) << endl;
	cout << "f1*=0.5" << endl << (f1 *= 0.5) << endl;
	cout << "f1*=f1" << endl << (f1 *= f1) << endl;


	cout << "identity*m2=" << endl << Mat<int>::eye(2) * m2 << endl;
	cout << "ones*m2=" << endl << Mat<int>::ones(2) * m2  << endl;
}