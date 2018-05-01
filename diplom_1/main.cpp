// ������. ��� � �������������� �� ����������, ������������ ���������� � ����������
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>

using namespace std;

//global
const int S = 134862;

void create_matr(ifstream& file, double *matr, double size) {
	for (int i = 0; i < size; ++i) {
		file >> matr[i];
	}
}

//bad idea
void create_big_matr(ifstream& file, double *matr) {
	int mid_count = 0;
	int low_count = 0;
	int lower_count = 0;
	for (int i = 0; i < S; ++i) {
				if (i == 0) {
					file >> matr[0];
					file >> matr[1];
					file >> matr[247];
					file >> matr[494];
					file >> matr[495];
					file >> matr[496];
					file >> matr[497];
				}
				else if (i > 0 && i < 247) {
					file >> matr[i*S + i - 1];
					file >> matr[i*S + i];
					file >> matr[i*S + i + 1];
					file >> matr[i*S + i + 247];
					file >> matr[i*S + i + 494];
					file >> matr[i*S + i + 495];
					file >> matr[i*S + i + 496];
				}
				else if (i > 246 && i < 494) {
					file >> matr[i*S + mid_count];
					file >> matr[i*S + mid_count + 246];
					file >> matr[i*S + mid_count + 247];
					file >> matr[i*S + mid_count + 248];
					file >> matr[i*S + mid_count + 494];
					file >> matr[i*S + mid_count + 741];
					file >> matr[i*S + mid_count + 742];
					mid_count++;
				}
				else if (i > 493 && i < 134368) {
					file >> matr[i*S + low_count];
					file >> matr[i*S + low_count + 247];
					file >> matr[i*S + low_count + 493];
					file >> matr[i*S + low_count + 494];
					file >> matr[i*S + low_count + 495];
					file >> matr[i*S + low_count + 741];
					file >> matr[i*S + low_count + 988];
					low_count++;
				}
				else if (i > 134367 && i < 134615) {
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

					file >> matr[i*S + lower_count];
				}
				else if (i > 134614 && i < 134861) {
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
				else if (i == 134861) {
					result << "0 0 0 1" << " ";
					A1 >> tmp;
					result << tmp << " ";
					A2 >> tmp;
					result << tmp << " ";
					A3 >> tmp;
					result << tmp << endl;
				}
	}
}


//size = sizeof(matr) ?(!warn about sizeof!)?
void show_matr(double *matr, double size) {
	for (int i = 0; i < size; ++i) {
		cout << matr[i] << "; ";
	}
	cout << endl;
}

//add res_matr?
void mult_matr_on_num(double *matr, double num, double *res_matr, double size) {
	for (int i = 0; i < size; ++i) {
		res_matr[i] = matr[i] * num;

		//debug
		//cout << "mult matr on num debug: " << res_matr[i] << " = " << matr[i] << " * " << num << endl;
	}
}

void mult_A_on_alpha_E(double *matr, double alpha, double *res_matr) {
	for (int i = 0; i < S; ++i) {
		if (i == 0) {
			matr[i] = matr[i] + alpha;
		}
		else if (i > 0 && i < 247) {
			matr[i * 7 + 1] = matr[i] + alpha;
		}
		else if (i > 246 && i < 494) {
			matr[i * 7 + 2] = matr[i] + alpha;
		}
		else if (i > 493 && i < 134368) {
			matr[i * 7 + 3] = matr[i] + alpha;
		}
		else if (i > 134367 && i < 134615) {
			matr[i * 7 + 4] = matr[i] + alpha;
		}
		else if (i > 134614 && i < 134861) {
			matr[i * 7 + 5] = matr[i] + alpha;
		}
		else if (i == 134861) { //944034
			matr[i * 7 + 6] = matr[i] + alpha;
		}
	}
}

void sum_matr(double *matr_one, double *matr_two, double *res_matr, double size) {
	for (int i = 0; i < size; ++i) {
		res_matr[i] = matr_one[i] + matr_two[i];
	}
}


//warn about size! (now for vector's size)
void dif_matr(double *matr_one, double *matr_two, double *res_matr) {
	for (int i = 0; i < 3; ++i) {
		res_matr[i] = matr_one[i] - matr_two[i];
		
		//debug
		//cout << "debug: " << res_matr[i] << " = " << matr_one[i] << " - " << matr_two[i] << endl;
	}
}

void matr_on_vec(double *matr, double *vec, double *res_vec) {
	for (int i = 0; i < S; ++i) {
		for (int j = 0; j < 7; ++j) {
			res_vec[i] += matr[7*i + j] * vec[j];

			//debug
			//cout << "i: " << i << " j: " << j  << " res_vec: " << res_vec[i] << " matr: " << matr[3 * i + j] << " vec: " << vec[j] <<  endl;
		}
	}
}

double vec_on_vec(double *vec_one, double *vec_two) {
	double res = 0;
	for (int i = 0; i < S; ++i) {
		res += vec_one[i] * vec_two[i];
	}
	//debug
	//cout << res << endl;
	//system("pause");
	
	return res;
}

double norm_vec(double *vec) {
	double res = 0;
	for (int i = 0; i < S; ++i) {
		res += pow(vec[i], 2);
	}
	return sqrt(res);
}

//size?
void copy_vec(double *matr_one, double *matr_two) {
	for (int i = 0; i < S; ++i) {
		matr_one[i] = matr_two[i];
	}
}

void nullify(double *vec) {
	for (int i = 0; i < S; ++i) {
		vec[i] = 0;
	}
}
//stop deal (rework!)
double stop_deal(double *rA, double *x_k1, double *F) {
	double *stop = NULL;
	double *tmp = NULL;
	double stop_norm = 0;
	double stop_norm2 = 0;
	//rA*x_k
	matr_on_vec(rA, x_k1, stop);
	//rA*x_k - b
	dif_matr(stop, F, tmp);
	//norm(rA*x_k - b)
	stop_norm = norm_vec(tmp);
	//norm(b)
	stop_norm2 = norm_vec(F);
	//stop_deal = stop_norm / stop_norm2
	return (stop_norm / stop_norm2);
}

void CGMR(double *A, double *F) {
	const double Eps = 0.1;

	//const double l_size = 9;
	double *rA = new double[S * 7];
	double *x = new double[S];
	double *r = new double[S];
	double *p = new double[S];

	//coef in cycle
	double ak = 0;
	double bk = 0;
	
	//const coef a
	double alpha = 0;

	//tmp
	double *matr_tmp = new double[S * 7];
	double *tmp = new double[S];
	double *r_k1 = new double[S];
	double *x_k1 = new double[S];
	double *p_k1 = new double[S];
	double *stop = new double[S];
	double stop_norm = 0;
	double stop_norm2 = 0;
	double norm = 0;
	double norm2 = 0;
	double vov = 0;
	double stop_eps = 1;

	for (int i = 0; i < S; ++i) {
		x[i] = r[i] = p[i] = 0;
		stop[i] = tmp[i] = r_k1[i] = x_k1[i] = p_k1[i] = 0;
	}

	//making const rA
	//rA = A + alpha*E
	mult_A_on_alpha_E(A, alpha, rA);

	//r = b-rA*x0	|rA*x0 = 0| |r = b|
	copy_vec(r, F);
	//p = r;
	copy_vec(p, r);

	//stop deal before cycle
	//rA*x_k						| = 0 |
	//matr_on_vec(rA, x_k1, stop);
	//rA*x_k - b					| -b |
	//dif_matr(stop, F, tmp);
	//norm(rA*x_k - b)				| norm(-b)|
	//stop_norm = norm_vec(F);
	//norm(b)
	//stop_norm2 = norm_vec(F);
	//stop_deal = stop_norm / stop_norm2	|always 1|
	//stop_eps = (stop_norm / stop_norm2);
	//stop_eps = 1;


	//mb norm like scalar not norm (|| ||) itself - yep
	while( !(stop_eps < Eps) ) {
		cout << "Start of iter" << endl;
		
		//ak = norm(r_k)^2 / (rA*p_k, p_k)
		//norm(r_k)^2
			//norm = pow(norm_vec(r), 2);
		norm = vec_on_vec(r, r);
		//rA*p_k
		matr_on_vec(rA, p, tmp);
		//(rA*p_k, p_k)
		vov = vec_on_vec(tmp, p);
		//ak = ...
		ak = (isnan(norm / vov) ? 0 : (norm / vov));
		
		nullify(tmp);
		//r_k+1 = r_k - ak*rA*p_k
		//rA*p_k
		matr_on_vec(rA, p, tmp);
			//cout << "rA*pk: " << endl;
			//show_matr(tmp, 3);
		//ak*(rA*p_k)
		mult_matr_on_num(tmp, ak, tmp, 3);
			//cout << "ak*rA*pk: " << endl;
			//show_matr(tmp, 3);
		//r_k+1 = r_k - ...
		dif_matr(r, tmp, r_k1);

		//cout << "r_k: " << endl;
		//show_matr(r, 3);
		//cout << "r_k+1:" << endl;
		//show_matr(r_k1, 3);

		nullify(tmp);
		//x_k+1 = x_k + ak*p_k
			//tmp = p;
			//copy_vec(tmp, p);
		//ak*p_k
		mult_matr_on_num(p, ak, tmp, 3);
		//x_k+1 = x_k + ...
		sum_matr(x, tmp, x_k1, 3);

		//bk = norm(r_k+1)^2 / norm(r_k)^2
		//norm(r_k+1)^2
			//norm = pow(norm_vec(r_k1), 2);
		norm = vec_on_vec(r_k1, r_k1);
		//norm(r_k)^2
			//norm2 = pow(norm_vec(r), 2);
		norm2 = vec_on_vec(r, r);
		//bk = ...
		bk = (isnan(norm / norm2) ? 0 : norm / norm2);

		nullify(tmp);
		//p_k+1 = r_k+1 + bk*p_k
		//bk*p_k
			//tmp = p;
			//copy_vec(tmp, p);
		mult_matr_on_num(p, bk, tmp, 3);
		//p_k+1 = r_k+1 + ...
		sum_matr(r_k1, tmp, p_k1, 3);

		//cout << "p_k: " << endl;
		//show_matr(p, 3);
		//cout << "p_k+1: " << endl;
		//show_matr(p_k1, 3);
		
		//show results
		cout << "coefs:" << endl << "alpha: " << alpha << "; ak: " << ak << "; bk: " << bk << endl;

		cout << "x_k1: " << endl;
		show_matr(x_k1, 3);

		nullify(tmp);
		//stop deal
		//rA*x_k
		matr_on_vec(rA, x_k1, stop);
		//rA*x_k - b
		dif_matr(stop, F, tmp);
		//norm(rA*x_k - b)
		stop_norm = norm_vec(tmp);
		//norm(b)
		stop_norm2 = norm_vec(F);
		//stop_deal = stop_norm / stop_norm2
		stop_eps = (stop_norm / stop_norm2);

		cout << "Stop deal: " << stop_eps << " " << ((stop_eps < Eps) ? "<" : ">") << " " << Eps << endl;
		//boolalpha << 

		//copy for next iter
		//p = p_k1;
		copy_vec(p, p_k1);
		//x = x_k1;
		copy_vec(x, x_k1);
		//r = r_k1;
		copy_vec(r, r_k1);
		//clear stop
		nullify(stop);
		nullify(p_k1);
		nullify(x_k1);
		nullify(r_k1);
		cout << "End of iter" << endl << endl;
	}

}

void create_big_matrix() {
	
}

void main() {

	//const double matr_size = 9;
	//double *matr_A = new double[matr_size];
	//double *matr_E = new double[matr_size];
	//double *matr_F = new double[matr_size / 3];

	//MAKE BIG MATRIX WITH POSITIONS BUT ONLY WITH NECCESERY VALUES

	double *matr_A = new double[S * 7];
	double *matr_F = new double[S];

	ifstream A;
	ifstream F;
	ifstream E;

	clock_t begin = clock();

	A.open("result.dat");	// 134862 * 3
	F.open("F.dat");		// 134862

	create_matr(A, matr_A, S * 7);
	create_matr(F, matr_F, S);

	A.close();
	F.close();

	CGMR(matr_A, matr_F);

	clock_t end = clock();
	double time = double(end - begin) / CLOCKS_PER_SEC;
	cout << "Time of calc: " << time << endl;

	system("pause");
}