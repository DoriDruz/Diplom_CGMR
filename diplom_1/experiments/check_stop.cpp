//Проверка конечной невязки по result_example.txt
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <iomanip>

using namespace std;

//global
const int S = 134862;

void create_matr(ifstream& file, double *matr, double size) {
	for (int i = 0; i < size; ++i) {
		file >> setprecision(16) >> matr[i];
	}
}

void matr_on_vec(double *matr, double *vec, double *res_vec) {
	int second_count = 0;
	int third_count = 0;
	int fourth_count = 0;
	int fifth_count = 0;
	int sixth_count = 0;

	for (int i = 0; i < S; ++i) {
		if (i == 0) {
			res_vec[i] = \
				matr[0] * vec[0] \
				+ matr[1] * vec[1] \
				+ matr[2] * vec[247] \
				+ matr[3] * vec[494];
		}
		else if (i > 0 && i < 247) {
			res_vec[i] = \
				matr[7 * i] * vec[second_count] \
				+ matr[7 * i + 1] * vec[second_count + 1] \
				+ matr[7 * i + 2] * vec[second_count + 2] \
				+ matr[7 * i + 3] * vec[second_count + 248] \
				+ matr[7 * i + 4] * vec[second_count + 495];
			second_count++;
		}
		else if (i > 246 && i < 494) {
			res_vec[i] = \
				matr[7 * i] * vec[third_count] \
				+ matr[7 * i + 1] * vec[third_count + 246] \
				+ matr[7 * i + 2] * vec[third_count + 247] \
				+ matr[7 * i + 3] * vec[third_count + 248] \
				+ matr[7 * i + 4] * vec[third_count + 494] \
				+ matr[7 * i + 5] * vec[third_count + 741];
			third_count++;
		}
		else if (i > 493 && i < 134368) {
			res_vec[i] = \
				matr[7 * i] * vec[fourth_count] \
				+ matr[7 * i + 1] * vec[fourth_count + 247] \
				+ matr[7 * i + 2] * vec[fourth_count + 493] \
				+ matr[7 * i + 3] * vec[fourth_count + 494] \
				+ matr[7 * i + 4] * vec[fourth_count + 495] \
				+ matr[7 * i + 5] * vec[fourth_count + 741] \
				+ matr[7 * i + 6] * vec[fourth_count + 988];
			fourth_count++;
		}
		else if (i > 134367 && i < 134615) {
			res_vec[i] = \
				matr[7 * i + 1] * vec[fifth_count + 133874] \
				+ matr[7 * i + 2] * vec[fifth_count + 134121] \
				+ matr[7 * i + 3] * vec[fifth_count + 134367] \
				+ matr[7 * i + 4] * vec[fifth_count + 134368] \
				+ matr[7 * i + 5] * vec[fifth_count + 134369] \
				+ matr[7 * i + 6] * vec[fifth_count + 134615];
			fifth_count++;
		}
		else if (i > 134614 && i < 134861) {
			res_vec[i] = \
				matr[7 * i + 2] * vec[sixth_count + 134121] \
				+ matr[7 * i + 3] * vec[sixth_count + 134368] \
				+ matr[7 * i + 4] * vec[sixth_count + 134614] \
				+ matr[7 * i + 5] * vec[sixth_count + 134615] \
				+ matr[7 * i + 6] * vec[sixth_count + 134616];
			sixth_count++;
		}
		else if (i == 134861) {
			res_vec[i] = \
				matr[7 * i + 3] * vec[134367] \
				+ matr[7 * i + 4] * vec[134614] \
				+ matr[7 * i + 5] * vec[134860] \
				+ matr[7 * i + 6] * vec[134861];
		}
	}
}

void dif_vec(double *vec_one, double *vec_two, double *res_vec) {
	for (int i = 0; i < S; ++i) {
		res_vec[i] = vec_one[i] - vec_two[i];
	}
}

double norm_vec(double *vec) {
	double res = 0;
	for (int i = 0; i < S; ++i) {
		res += pow(vec[i], 2);
	}
	return sqrt(res);
}

void nev(double *A, double *F, double *X) {
	
	double *tmp = new double[S];
	double *stop_tmp = new double[S];
	double stop_num = 1;
	double stop_up = 0;
	double stop_down = 0;


	//Условие останова
	// ||A*x_k - b|| / ||b|| < Eps
	//A*x_k
	matr_on_vec(A, X, tmp);
	//A*x_k - b
	dif_vec(tmp, F, stop_tmp);
	//||A*x_k - b||
	stop_up = norm_vec(stop_tmp);
	//||b||
	stop_down = norm_vec(F);
	//||..|| / ||..||
	stop_num = 1;
	stop_num = stop_up / stop_down;

	//Показываем невязку
	cout << "Nevyazka: " << stop_num << endl;

	delete(tmp);
	delete(stop_tmp);
}

void main() {
	double *matr_A = new double[S * 7];
	double *matr_F = new double[S];
	double *matr_X = new double[S];
	double *matr_X2 = new double[S];
	
	ifstream A;
	ifstream F;
	ifstream X;
	ifstream X2;
	
	A.open("A_with_01.dat");
	F.open("F.dat");
	X.open("X_1.dat");
	X2.open("X_ars.dat");
	
	create_matr(A, matr_A, S * 7);
	create_matr(F, matr_F, S);
	create_matr(X, matr_X, S);
	create_matr(X2, matr_X2, S);
	
	if (A.is_open() && F.is_open() && X.is_open()) {
		nev(matr_A, matr_F, matr_X);
		nev(matr_A, matr_F, matr_X2);
	}
	else {
		cout << "Error reading files" << endl;
	}

	A.close();
	F.close();
	X.close();

	delete(matr_A);
	delete(matr_F);
	delete(matr_X);
	delete(matr_X2);

	system("pause");
}