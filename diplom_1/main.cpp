// Итерац. МСГ с регуляризацией на процессоре[, многоядерном процессоре и видеокарте]
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <iomanip>

using namespace std;

//global
const int S = 134862;

void add_in_file(double nev) {
	ofstream nevyazka;
	nevyazka.open("nevyazka.dat", ios_base::app);
	nevyazka << nev << endl;
	nevyazka.close();
}

void write_in_file(double *vec, double size, string name) {
	clock_t begin_write = clock();

	cout << "Start write in file" << endl;
	ofstream result_file;
	result_file.open(name, ofstream :: trunc); //x_res_from_main.txt

	for (int i = 0; i < size; ++i) {
		result_file << setprecision(16) << vec[i] << endl;
	}

	result_file.close();

	clock_t end_write = clock();
	double write_time = double(end_write - begin_write) / CLOCKS_PER_SEC;
	cout << "Stop write in file: " << write_time << endl;
}

void create_matr(ifstream& file, double *matr, double size) {
	for (int i = 0; i < size; ++i) {
		file >> setprecision(16) >> matr[i];
	}
}

/*void create_big_matr(ifstream& file, double *matr) {
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
}*/

void show_vec(double *matr, double size) {
	for (int i = 0; i < size; ++i) {
		cout << setprecision(16) << matr[i] << "; ";
	}
	cout << endl;
}

void mult_vec_on_num(double *vec, double num, double *res_vec) {
	for (int i = 0; i < S; ++i) {
		res_vec[i] = vec[i] * num;
	}
}

/*useless
void mult_matr_on_num(double *matr, double num, double *res_matr, double size) {
	for (int i = 0; i < size; ++i) {
		res_matr[i] = matr[i] * num;

		//debug
		//cout << "mult matr on num debug: " << res_matr[i] << " = " << matr[i] << " * " << num << endl;
	}
}*/

//function changing original matrix A
void mult_A_on_alpha_E(double *matr, double alpha) {
	for (int i = 0; i < S; ++i) {
		if (i == 0) {
			matr[i] = matr[i] + alpha;
		}
		else if (i > 0 && i < 247) {
			matr[i * 7 + 1] = matr[i * 7 + 1] + alpha;
		}
		else if (i > 246 && i < 494) {
			matr[i * 7 + 2] = matr[i * 7 + 2] + alpha;
		}
		else if (i > 493 && i < 134368) {
			matr[i * 7 + 3] = matr[i * 7 + 3] + alpha;
		}
		else if (i > 134367 && i < 134615) {
			matr[i * 7 + 4] = matr[i * 7 + 4] + alpha;
		}
		else if (i > 134614 && i < 134861) {
			matr[i * 7 + 5] = matr[i * 7 + 5] + alpha;
		}
		else if (i == 134861) { //944034
			matr[i * 7 + 6] = matr[i * 7 + 6] + alpha;
		}
	}
}

void sum_vec(double *vec_one, double *vec_two, double *res_vec) {
	for (int i = 0; i < S; ++i) {
		res_vec[i] = vec_one[i] + vec_two[i];
	}
}

/*useless
void sum_matr(double *matr_one, double *matr_two, double *res_matr, double size) {
	for (int i = 0; i < size; ++i) {
		res_matr[i] = matr_one[i] + matr_two[i];
	}
}*/

void dif_vec(double *vec_one, double *vec_two, double *res_vec) {
	for (int i = 0; i < S; ++i) {
		res_vec[i] = vec_one[i] - vec_two[i];
	}
}

/*warn about size! (now for vector's size)
void dif_matr(double *matr_one, double *matr_two, double *res_matr) {
	for (int i = 0; i < 3; ++i) {
		res_matr[i] = matr_one[i] - matr_two[i];
		
		//debug
		//cout << "debug: " << res_matr[i] << " = " << matr_one[i] << " - " << matr_two[i] << endl;
	}
}*/

void matr_on_vec(double *matr, double *vec, double *res_vec) {
	int second_count = 0;
	int third_count = 0;
	int fourth_count = 0;
	int fifth_count = 0;
	int sixth_count = 0;

	//debug
	//ofstream debug_f;
	//debug_f.open("debug.dat", ios_base::app);

	for (int i = 0; i < S; ++i) {
		if (i == 0) {
			res_vec[i] = \
				matr[0] * vec[0] \
				+ matr[1] * vec[1] \
				+ matr[2] * vec[247] \
				+ matr[3] * vec[494];

			/*debug
			debug_f << i << ": " << res_vec[i] << " = " \
				<< matr[0] << " * " << vec[0] << " + " \
				<< matr[1] << " * " << vec[1] << " + " \
				<< matr[2] << " * " << vec[247] << " + " \
				<< matr[3] << " * " << vec[494] \
				<< endl;
				*/
		}
		else if (i > 0 && i < 247) {
			res_vec[i] = \
				matr[7 * i] * vec[second_count] \
				+ matr[7 * i + 1] * vec[second_count + 1] \
				+ matr[7 * i + 2] * vec[second_count + 2] \
				+ matr[7 * i + 3] * vec[second_count + 248] \
				+ matr[7 * i + 4] * vec[second_count + 495];

			/*debug
			if(i == 1) {
			debug_f << i << ": " << res_vec[i] << " = " \
				<< matr[7 * i] << " * " << vec[second_count] << " + " \
				<< matr[7 * i + 1] << " * " << vec[second_count + 1] << " + " \
				<< matr[7 * i + 2] << " * " << vec[second_count + 2] << " + " \
				<< matr[7 * i + 3] << " * " << vec[second_count + 248] << " + " \
				<< matr[7 * i + 4] << " * " << vec[second_count + 495] \
				<< endl;
			}*/

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

			/*debug
			if (i == 247) {
				debug_f << i << ": " << res_vec[i] << " = " \
					<< matr[7 * i] << " * " << vec[third_count] << " + " \
					<< matr[7 * i + 1] << " * " << vec[third_count + 246] << " + " \
					<< matr[7 * i + 2] << " * " << vec[third_count + 247] << " + " \
					<< matr[7 * i + 3] << " * " << vec[third_count + 248] << " + " \
					<< matr[7 * i + 4] << " * " << vec[third_count + 494] << " + " \
					<< matr[7 * i + 5] << " * " << vec[third_count + 741] \
					<< endl;
			}*/

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

			/*debug
			if (i == 494) {
				debug_f << i << ": " << res_vec[i] << " = " \
					<< matr[7 * i] << " * " << vec[fourth_count] << " + " \
					<< matr[7 * i + 1] << " * " << vec[fourth_count + 247] << " + " \
					<< matr[7 * i + 2] << " * " << vec[fourth_count + 493] << " + " \
					<< matr[7 * i + 3] << " * " << vec[fourth_count + 494] << " + " \
					<< matr[7 * i + 4] << " * " << vec[fourth_count + 495] << " + " \
					<< matr[7 * i + 5] << " * " << vec[fourth_count + 741] << " + " \
					<< matr[7 * i + 6] << " * " << vec[fourth_count + 988] \
					<< endl;
			}*/

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

			/*debug
			if (i == 134368) {
				debug_f << i << ": " << res_vec[i] << " = " \
					<< matr[7 * i + 1] << " * " << vec[fifth_count + 133874] << " + " \
					<< matr[7 * i + 2] << " * " << vec[fifth_count + 134121] << " + " \
					<< matr[7 * i + 3] << " * " << vec[fifth_count + 134367] << " + " \
					<< matr[7 * i + 4] << " * " << vec[fifth_count + 134368] << " + " \
					<< matr[7 * i + 5] << " * " << vec[fifth_count + 134369] << " + " \
					<< matr[7 * i + 6] << " * " << vec[fifth_count + 134615] \
					<< endl;
			}*/

			fifth_count++;
		}
		else if (i > 134614 && i < 134861) {
			res_vec[i] = \
				matr[7 * i + 2] * vec[sixth_count + 134121] \
				+ matr[7 * i + 3] * vec[sixth_count + 134368] \
				+ matr[7 * i + 4] * vec[sixth_count + 134614] \
				+ matr[7 * i + 5] * vec[sixth_count + 134615] \
				+ matr[7 * i + 6] * vec[sixth_count + 134616];

			/*debug
			if (i == 134615) {
				debug_f << i << ": " << res_vec[i] << " = " \
					<< matr[7 * i + 2] << " * " << vec[sixth_count + 134121] << " + " \
					<< matr[7 * i + 3] << " * " << vec[sixth_count + 134368] << " + " \
					<< matr[7 * i + 4] << " * " << vec[sixth_count + 134614] << " + " \
					<< matr[7 * i + 5] << " * " << vec[sixth_count + 134615] << " + " \
					<< matr[7 * i + 6] << " * " << vec[sixth_count + 134616] \
					<< endl;
			}*/

			sixth_count++;
		}
		else if (i == 134861) { //last_element_position = 944034
			res_vec[i] = \
				matr[7 * i + 3] * vec[134367] \
				+ matr[7 * i + 4] * vec[134614] \
				+ matr[7 * i + 5] * vec[134860] \
				+ matr[7 * i + 6] * vec[134861];

			/*debug
			debug_f << i << ": " << res_vec[i] << " = " \
				<< matr[7 * i + 3] << " * " << vec[134367] << " + " \
				<< matr[7 * i + 4] << " * " << vec[134614] << " + " \
				<< matr[7 * i + 5] << " * " << vec[134860] << " + " \
				<< matr[7 * i + 6] << " * " << vec[134861] \
				<< endl;
		}*/

		//debug
		//cout << "i: " << i << " j: " << j  << " res_vec: " << res_vec[i] << " matr: " << matr[3 * i + j] << " vec: " << vec[j] <<  endl;

		}

		//debug
		//debug_f.close();
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
/*stop deal (rework!)
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
}*/

void CGMR(double *A, double *F, clock_t begin_algo) {
	const double Eps = 0.1;

	//const double l_size = 9;
	//double *rA = new double[S * 7];
	double *x = new double[S];
	double *r = new double[S];
	double *p = new double[S];

	//coef in cycle
	double ak = 0;
	double bk = 0;
	
	//const coef a
	double alpha = 0;

	//tmp
	double *tmp = new double[S];
	double *rApk = new double[S];
	double *stop = new double[S];
	double *r_k1 = new double[S];
	double *x_k1 = new double[S];
	double *p_k1 = new double[S];
	double stop_norm = 0;
	double stop_norm2 = 0;
	double norm = 0;
	double norm2 = 0;
	double vov = 0;
	double stop_eps = 1;

	int count_tmp = 0;

	//clear file with nevyazka
	ofstream nevyazka_file;
	nevyazka_file.open("nevyazka.dat", ofstream::trunc);
	nevyazka_file.close();

	//clear file with debug
	ofstream debug_f;
	debug_f.open("debug.dat", ofstream::trunc);
	debug_f.close();

	//fill all arrays with zeroes
	clock_t begin_zero_filling = clock();
	for (int i = 0; i < S; ++i) {
		x[i] = r[i] = p[i] = 0;
		stop[i] = tmp[i] = r_k1[i] = x_k1[i] = p_k1[i] = 0;
	}
	clock_t end_zero_filling = clock();

	double z_time = double(end_zero_filling - begin_zero_filling) / CLOCKS_PER_SEC;
	cout << "Zero_filling_preperation: " << z_time << endl;

	clock_t begin_make_prep = clock();
	
	//preps before algo

	//making const rA
	//rA = A + alpha*E
	mult_A_on_alpha_E(A, alpha);

	//r = b-rA*x0	|rA*x0 = 0| |r = b|
	copy_vec(r, F);
	//p = r;
	copy_vec(p, r);
	
	clock_t end_make_prep = clock();

	double p_time = double(end_make_prep - begin_make_prep) / CLOCKS_PER_SEC;
	cout << "Preperation rA, r, p: " << p_time << endl << endl;

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

	//WHY NEVYAZKA INCREASING ? - Different data in arrays? Make zero iter out of cycle?

	while( !(stop_eps < Eps) && count_tmp < S / 16) {
		//system("cls");

		cout << "Start of iter" << endl;
		clock_t begin_CGMR = clock();

		//ak = norm(r_k)^2 / (rA*p_k, p_k)
		//norm(r_k)^2
			//norm = pow(norm_vec(r), 2);
		norm = vec_on_vec(r, r);
		//rA*p_k
		matr_on_vec(A, p, rApk);
		//(rA*p_k, p_k)
		vov = vec_on_vec(rApk, p);
		//ak = ...
		ak = (isnan(norm / vov) ? 0 : (norm / vov));
		
		//r_k+1 = r_k - ak*rA*p_k
		//rA*p_k
		//matr_on_vec(A, p, tmp);
			//cout << "rA*pk: " << endl;
			//show_matr(tmp, 3);
		//ak*(rA*p_k)
		mult_vec_on_num(rApk, ak, tmp);
			//cout << "ak*rA*pk: " << endl;
			//show_matr(tmp, 3);
		//r_k+1 = r_k - ...
		dif_vec(r, tmp, r_k1);

		//cout << "r_k: " << endl;
		//show_matr(r, 3);
		//cout << "r_k+1:" << endl;
		//show_matr(r_k1, 3);

		nullify(tmp);
		//x_k+1 = x_k + ak*p_k
			//tmp = p;
			//copy_vec(tmp, p);
		//ak*p_k
		mult_vec_on_num(p, ak, tmp);
		//x_k+1 = x_k + ...
		sum_vec(x, tmp, x_k1);

		//bk = norm(r_k+1)^2 / norm(r_k)^2
		//norm(r_k+1)^2
			//norm = pow(norm_vec(r_k1), 2);
		norm = vec_on_vec(r_k1, r_k1);
		//norm(r_k)^2
			//norm2 = pow(norm_vec(r), 2);
		norm2 = vec_on_vec(r, r);
		//bk = ...
		bk = (isnan(norm / norm2) ? 0 : (norm / norm2));

		nullify(tmp);
		//p_k+1 = r_k+1 + bk*p_k
		//bk*p_k
			//tmp = p;
			//copy_vec(tmp, p);
		mult_vec_on_num(p, bk, tmp);
		//p_k+1 = r_k+1 + ...
		sum_vec(r_k1, tmp, p_k1);

		//cout << "p_k: " << endl;
		//show_matr(p, 3);
		//cout << "p_k+1: " << endl;
		//show_matr(p_k1, 3);
		
		//show results
		cout << "alpha: " << alpha << "; ak: " << ak << "; bk: " << bk << endl;

		cout << "x_k1[0-9]: " << endl;
		show_vec(x_k1, 10);

		//write_in_file(x_k1, 1000);

		nullify(tmp);
		//stop deal
		//rA*x_k
		matr_on_vec(A, x_k1, stop);
		//rA*x_k - b
		dif_vec(stop, F, tmp);
		//norm(rA*x_k - b)
		stop_norm = norm_vec(tmp);
			//stop_norm = sqrt(vec_on_vec(tmp, tmp));
		//norm(b) - можно вынести как константу
		stop_norm2 = norm_vec(F);
			//stop_norm2 = sqrt(vec_on_vec(F, F));
		//stop_deal = stop_norm / stop_norm2
		stop_eps = 0;
		stop_eps = (isnan(stop_norm / stop_norm2) ? 0 : (stop_norm / stop_norm2));

		add_in_file(stop_eps);

		cout << "Stop deal: " << stop_eps << " " << ((stop_eps < Eps) ? "<" : ">") << " " << Eps << endl;
		//boolalpha << 

		//copy for next iter
		//p = p_k1;
		copy_vec(p, p_k1);
		//x = x_k1;
		copy_vec(x, x_k1);
		//r = r_k1;
		copy_vec(r, r_k1);
		//clear vectors
		nullify(stop);
		nullify(p_k1);
		//nullify(x_k1);
		nullify(r_k1);
		nullify(tmp);
		nullify(rApk);
		cout << "End of iter" << endl;
				
		clock_t end_CGMR = clock();

		double CGMR_time = double(end_CGMR - begin_CGMR) / CLOCKS_PER_SEC;
		cout << "Runtime of iter: " << CGMR_time << endl << endl;

		clock_t end_algo = clock();
		double algo_time = double(end_algo - begin_algo) / CLOCKS_PER_SEC;
		cout << "Runtime of algo: " << algo_time << endl;

		count_tmp++;
		cout << "Number of iter: " << count_tmp << " / " << S << endl << endl;
		
		//tmp - delete last
		if (count_tmp == S / 32) {
			write_in_file(x_k1, 1000, "debug.dat");
		}

		if (count_tmp == S / 16) {
			write_in_file(x_k1, S, "x_res_from_main.txt");
		}
		else {
			nullify(x_k1);
		}
	}
	delete(x);
	delete(r);
	delete(p);
	delete(x_k1);
	delete(r_k1);
	delete(p_k1);
	delete(tmp);
	delete(rApk);
	delete(stop);

}

//tmp
/*void tmp_func() {
	double *matr_1 = new double[S];
	double *matr_2 = new double[S];
	double *matr_tmp = new double[S];
	double num_tmp = 0;

	for (int i = 0; i < S; ++i) {
		matr_1[i] = 1;
		matr_2[i] = 1;
	}

	num_tmp = vec_on_vec(matr_1, matr_2);
	write_in_file(matr_tmp, S);
	cout << num_tmp << endl;
}*/

void main() {
	clock_t begin_algo = clock();

	double *matr_A = new double[S * 7];
	double *matr_F = new double[S];
	
	/*double *matr_tmp = new double[S];
	double *matr_tmp_res = new double[S];

	for (int i = 0; i < S; ++i) {
		matr_tmp[i] = 1;
	}*/

	ifstream A;
	ifstream F;
	//ifstream E;

	clock_t begin = clock();

	A.open("result.dat");				// 134862 * 7	| result.dat
	F.open("F.dat");		// 134862		| F.dat

	create_matr(A, matr_A, S * 7);
	create_matr(F, matr_F, S);

	A.close();
	F.close();
	
	CGMR(matr_A, matr_F, begin_algo);

	clock_t end = clock();
	double time = double(end - begin) / CLOCKS_PER_SEC;
	cout << "Runtime: " << time << endl;

	delete(matr_A);
	delete(matr_F);

	system("pause");
}