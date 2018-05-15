#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <iomanip>
#include <algorithm>    // std::min_element, std::max_element

using namespace std;

//global
const int S = 134862;

void create_matr(ifstream& file, double *matr, double size) {
	for (int i = 0; i < size; ++i) {
		file >> setprecision(16) >> matr[i];
	}
}

void main() {
	double *matr_A = new double[S * 7];

	ifstream A;
	A.open("result.dat");
	create_matr(A, matr_A, S * 7);

	//cout << "Max element: " << max_element(matr_A, matr_A + S - 1) << endl;
	//cout << "Min element: " << min_element(matr_A, matr_A + S - 1) << endl;

	auto min_and_max = minmax_element(matr_A, matr_A + S - 1);

	cout << "Min: " << setprecision(16) << *min_and_max.first << endl << "Max: " << *min_and_max.second << endl;

	A.close();
	system("pause");
}