#include "CBLAB_Taegun.h"

int main() {
	
	
	std::string snp_path = "X.txt";  //X.txt
	std::string pheno_path = "Y.txt"; //Y.txt
	

//	inputMS(snp_path, pheno_path);	
	inputMS2(snp_path, pheno_path);
	

	//std::ifstream input_x("X_1.txt");
	//std::ifstream input_y("Y.txt");
	//std::ifstream input_k("K_1.txt");
	//std::ifstream input_posterior("posterior.txt");
	//std::ofstream NICE("NICE.txt");


	////필요 변수들 파악해둘것.
	//Eigen::MatrixXd X = read_mat(input_x, count_matrix_row(input_x), count_matrix_col(input_x)); // vector 형태
	//Eigen::MatrixXd K = read_mat(input_k, count_matrix_row(input_k), count_matrix_col(input_k));
	//Eigen::MatrixXd Y = read_mat(input_y, count_matrix_row(input_y), count_matrix_col(input_y));


	//NICE.precision(5);
	//emma(X, Y, K, NICE);


	//input_x.close();
	//input_y.close();
	//input_k.close();
	//input_posterior.close();
	
	
	
	NICE.close();
	return 0;
}