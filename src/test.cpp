#include "matrix.hpp"
#include <iostream>

using namespace std;

int main() {
	Mat<int> m1 = Mat<int>::ones(5);
	for(int i = 1; i <= 3; ++i) {
		m1(i, i) = 100000;
	}
	cout << "m1=" << endl;
	cout << m1 << endl;

	Mat<int> m2 = m1.select(2, 3, 2, 5);
	cout << "m2=" << endl;
	cout << m2 << endl;
	
	Mat<int> m3 = m1.select(1, 2, 2, 3);
	cout << "m3=" << endl;
	cout << m3 << endl;

	Mat<int> id = Mat<int>::ones(2);
	cout << "m3+one=" << endl;
	cout << (id + m3) << endl;

	cout << "identity + m3 - identity=" << endl;
	cout << (id + m3 - id) << endl;

	cout << "-m2=" << endl;
	cout << (-m2) << endl;

	cout << "------------------------------------" << endl;

	cout << "m3-=ones" << endl;
	m3 -= Mat<int>::ones(2);

	cout << m3 << endl;
	cout << m2 << endl;
	cout << "m2*9" << endl;
	cout << (m2*9) << endl;
	cout << (9*m2) << endl;
	cout << "m3*=3" << endl;
	m3 *= 3;

	cout << m3 << endl;
	cout << "---------------------------------------" << endl;
}