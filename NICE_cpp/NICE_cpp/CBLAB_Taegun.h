#pragma once

#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/LU"
#include "Eigen/Eigen/Eigenvalues"
#include "Eigen/unsupported/Eigen/MatrixFunctions"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <cmath>
#include <ctime>
#include <random>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/distributions/inverse_gaussian.hpp>
#include <boost/math/special_functions/beta.hpp>

static std::mt19937 ran;
static std::ofstream NICE("NICE.txt");

class REMLE {
public:
    double REML; double delta; double ve; double vg;
    REMLE() {
        REML = 0.0; delta = 0.0; ve = 0.0; vg = 0.0;
    }
};

class input_MS {
public:
    double cal_estimate(input_MS* pheno);
    double cal_RSE(input_MS* pheno);

    std::vector<double> element_vec;
    std::vector<double> residual_vec;

    double mean;
    double sum;
    double squared_residual_sum;
    // y = ax + b
    double slope; // a
    double intercept; // b
    int count;

    input_MS();
    input_MS(std::string arr, int check);
    input_MS(Eigen::MatrixXd& Y, int row_idx, int col);

};
class Metasoft {

};
class MetaSnp {
public:
    std::vector<double>* betas_; //later delete plz
    std::vector<double>* standardErrors_; //later delete plz

    MetaSnp();
    MetaSnp(int snp_count);
};


void emma_REMLE(Eigen::ArrayXd& y, Eigen::MatrixXd& x, Eigen::MatrixXd& K, Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>& eig_R, REMLE& REM, int indi);

int count_matrix_col(std::ifstream& matrix);
int count_matrix_row(std::ifstream& matrix);
void emma_eigen_L_wo_Z(Eigen::MatrixXd& Kinship, Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>& eig_L);
void emma_eigen_R_wo_Z(Eigen::MatrixXd& Kinship, Eigen::MatrixXd& X, Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>& eig_R, int indi);
double emma_delta_REML_dLL_wo_Z(double logdelta, Eigen::Map<Eigen::ArrayXd>& eig_value, Eigen::Map<Eigen::ArrayXd>& etas, Eigen::ArrayXd& etasq);
double emma_delta_REML_LL_wo_Z(double logdelta, Eigen::Map<Eigen::ArrayXd>& eig_value, Eigen::Map<Eigen::ArrayXd>& etas);
double uniroot_emma_delta_REML_dLL_wo_Z(double a, double b, Eigen::Map<Eigen::ArrayXd>& eig_value, Eigen::Map<Eigen::ArrayXd>& etas, Eigen::ArrayXd& etasq);
long double p_value(boost::math::students_t_distribution<double> t_dist, long double x);
void emma(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, Eigen::MatrixXd& K, std::ofstream& out);

double zscore(double x, double t_value);
long double p_value(boost::math::students_t_distribution<double> t_dist, long double x);
long double t_value(boost::math::students_t_distribution<double> t_dist, long double x);

void inputMS(std::string snp, std::string pheno);
void inputMS2(std::string snp, std::string pheno);

void computeMvaluesMCMC(std::vector<double>& betas, std::vector<double>& std_, int sample, int pheno_num, std::string X, Eigen::MatrixXd Y, int seed = 0);

double observationLogLikelihood(std::vector<double>& betas, std::vector<double>& std_, std::vector<int>& H1, int numH1);
double makeRandomDouble();
int makeRandomInteger(int number);

Eigen::MatrixXd read_mat(std::ifstream& input_file, int row, int col);
Eigen::MatrixXd read_mat(std::string X, int col);
Eigen::MatrixXd cov(Eigen::MatrixXd mat);
Eigen::MatrixXd normMe(Eigen::MatrixXd mat);
Eigen::MatrixXd cbind(Eigen::MatrixXd& a, Eigen::MatrixXd& b);

