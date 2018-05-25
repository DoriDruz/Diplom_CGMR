#include <iostream>
#include <fstream>
#include <string>
#include <ctime>

using namespace std;

void main() {
	clock_t begin = clock();
	string header;
	string tmp = "2";
	int size = 134862;

	ofstream result;

	result.open("F_with_two.dat", ofstream::trunc);
	if(result.is_open()){
		for (int i = 0; i < size; ++i) {
			result << tmp << endl;
		}
	}
	else {
		cout << "Error" << endl;
	}

	/*for (int i = 1; i <= size; ++i) {
		cout << i << " / " << size << endl;
		if (i == 1) {
			result << tmp << " ";
			
			result << tmp << " ";
			
			result << tmp << " ";
			result << "1 0 0 0" << endl;
		}
		else if (i >= 2 && i <= 247) {
			
			result << tmp << " ";
			
			result << tmp << " ";
			
			result << tmp << " ";
			
			result << tmp << " ";
			result << "1 0 0" << endl;
		}
		else if (i >= 248 && i <= 494) {
			
			result << tmp << " ";
			
			result << tmp << " ";
			
			result << tmp << " ";
			
			result << tmp << " ";
			
			result << tmp << " ";
			result << "1 0" << endl;
		}
		else if (i >= 495 && i <= 134368) {
			result << "1" << " ";
			
			result << tmp << " ";
			
			result << tmp << " ";
			
			result << tmp << " ";
			
			result << tmp << " ";
			
			result << tmp << " ";
			result << "1" << endl;
		}
		else if (i >= 134369 && i <= 134615) {
			result << "0 1" << " ";
			
			result << tmp << " ";
			
			result << tmp << " ";
			
			result << tmp << " ";
			
			result << tmp << " ";
			
			result << tmp << endl;
		}
		else if (i >= 134616 && i <= 134861) {
			result << "0 0 1" << " ";
			
			result << tmp << " ";
			
			result << tmp << " ";
			
			result << tmp << " ";
			
			result << tmp << endl;
		}
		else if (i == 134862) {
			result << "0 0 0 1" << " ";
			
			result << tmp << " ";
			
			result << tmp << " ";
			
			result << tmp << endl;
		}
	}*/

	result.close();
	
	clock_t end = clock();
	double time = double(end - begin) / CLOCKS_PER_SEC;
	cout << "Time of calc: " << time << endl;

	system("pause");
}