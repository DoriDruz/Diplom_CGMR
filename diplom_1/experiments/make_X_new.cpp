#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <iomanip>

using namespace std;

//global
const int S = 134862;

void create_X_new(ofstream& file, double size) {
	file.open("X_new.dat");
	double tmp = 0;
	double range = 0.00001;
	for (int i = 0; i < size; ++i) {
		cout << i + 1 << " / " << endl;
		tmp += range;
		file << setprecision(16) << tmp << endl;
	}
	file.close();
	cout << tmp << endl;
}

void main() {
	ofstream X_new;
	create_X_new(X_new, S);
	system("pause");
}