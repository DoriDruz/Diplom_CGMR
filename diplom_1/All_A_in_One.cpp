#include <iostream>
#include <fstream>
#include <string>
#include <ctime>

using namespace std;

void main() {
	string header;
	string tmp;
	int size = 134862;

	ifstream A1;
	ifstream A2;
	ifstream A3;
	ifstream A4;
	ifstream A5;
	fstream result;

	clock_t begin = clock();

	A1.open("A1.dat");
	A2.open("A2.dat");
	A3.open("A3.dat");
	A4.open("A4.dat");
	A5.open("A5.dat");

	getline(A1, header);
	getline(A2, header);
	getline(A3, header);
	getline(A4, header);
	getline(A5, header);

	result.open("result.dat");

	//REWORK FOR CORRECT NULLS?

	for (int i = 1; i <= size; ++i) {
		cout << i <<  " / " << size << endl;
		if (i == 1) {
			A3 >> tmp;
			result << tmp << " ";
			A4 >> tmp;
			result << tmp << " ";
			A5 >> tmp;
			result << tmp << " ";
			result << "1 0 0 0" << endl; 
		}
		else if (i >= 2 && i <= 247) {
			A2 >> tmp;
			result << tmp << " ";
			A3 >> tmp;
			result << tmp << " ";
			A4 >> tmp;
			result << tmp << " ";
			A5 >> tmp;
			result << tmp << " ";
			result << "1 0 0" << endl;
		}
		else if (i >= 248 && i <= 494) {
			A1 >> tmp;
			result << tmp << " ";
			A2 >> tmp;
			result << tmp << " ";
			A3 >> tmp;
			result << tmp << " ";
			A4 >> tmp;
			result << tmp << " ";
			A5 >> tmp;
			result << tmp << " ";
			result << "1 0" << endl;
		}
		else if (i >= 495 && i <= 134368) {
			result << "1" << " ";
			A1 >> tmp;
			result << tmp << " ";
			A2 >> tmp;
			result << tmp << " ";
			A3 >> tmp;
			result << tmp << " ";
			A4 >> tmp;
			result << tmp << " ";
			A5 >> tmp;
			result << tmp << " ";
			result << "1" << endl;
		}
		else if (i >= 134369 && i <= 134615) {
			result << "0 1" << " ";
			A1 >> tmp;
			result << tmp << " ";
			A2 >> tmp;
			result << tmp << " ";
			A3 >> tmp;
			result << tmp << " ";
			A4 >> tmp;
			result << tmp << " ";
			A5 >> tmp;
			result << tmp << endl;
		}
		else if (i >= 134616 && i <= 134861) {
			result << "0 0 1" << " ";
			A1 >> tmp;
			result << tmp << " ";
			A2 >> tmp;
			result << tmp << " ";
			A3 >> tmp;
			result << tmp << " ";
			A4 >> tmp;
			result << tmp << endl;
		}
		else if (i == 134862) {
			result << "0 0 0 1" << " ";
			A1 >> tmp;
			result << tmp << " ";
			A2 >> tmp;
			result << tmp << " ";
			A3 >> tmp;
			result << tmp << endl;
		}
	}

	result << "End";
	result.close();

	A1.close();
	A2.close();
	A3.close();
	A4.close();
	A5.close();

	clock_t end = clock();
	double time = double(end - begin) / CLOCKS_PER_SEC;
	cout << "Time of calc: " << time << endl;

	system("pause");
}