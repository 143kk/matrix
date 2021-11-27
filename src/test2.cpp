#include <iostream>
#include "matrix.hpp"

using namespace std;

int main() {
	Mat<short> m1(4, 3, 3);
	srand(time(0));
	for(int i = 1; i <= 4; ++i) {
		for(int j = 1; j <= 3; ++j) {
			for(int k = 0; k < 3; ++k) {
				m1(i, j, k) = rand()%256;
			}
		}
	}
	m1(2,2,1) = 0;
	Mat<short> m2 = m1;
	cout << m1 << endl;
	Mat<short> m3 = m1 + m2;
	cout << "m3 = m1*2" << endl;
	cout << m3 << endl;
	cout << (m1*2) << endl;

	cout << "m1*=2" << endl;
	m1*=2;
	cout << m1 << endl;

	Mat<short> m4 = m1.select(2, 3, 2, 3);
	cout << m4 << endl;
	auto m5 = m1.select(1, 2, 1, 2);
	cout << m5 << endl;
	cout << "m4+m5" << endl;
	cout << m4+m5 << endl;
	m4 += m5;
	cout << m4 << endl;
	cout << m1 << endl;
}