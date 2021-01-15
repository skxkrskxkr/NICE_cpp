#include "CBLAB_Taegun.h"

input_MS::input_MS() { //do not use
    mean = 0.0; sum = 0.0; count = 0; squared_residual_sum = 0.0; slope = 0.0; intercept = 0.0;
}

input_MS::input_MS(std::string arr, int check) {  //check 0 = snp, check 1 = pheno
    mean = 0.0; sum = 0.0; count = 0; squared_residual_sum = 0.0; slope = 0.0; intercept = 0.0;
    std::string token;
    std::stringstream stream;
    stream.str(arr);
    while (stream >> token) {
        element_vec.push_back(std::stod(token));
        sum += element_vec[count];
        count++;
    }
    stream.clear();
    mean = sum / count;

    if (check == 0) {
        for (int i = 0; i < count; i++) {
            residual_vec.push_back(element_vec[i] - mean);
            squared_residual_sum += residual_vec[i] * residual_vec[i];
        }
    }
    else if (check == 1) {
        for (int i = 0; i < count; i++) {
            residual_vec.push_back(element_vec[i] - mean);
        }
    }
    else {
        std::cout << "Check your data is snp or pheno" << std::endl;
        exit(1);
    }
}

input_MS::input_MS(Eigen::MatrixXd& Y, int row_idx, int col) {  //check 0 = snp, check 1 = pheno
    mean = 0.0; sum = 0.0; count = 0; squared_residual_sum = 0.0; slope = 0.0; intercept = 0.0;
    for(int i = 0; i < col; i++) {
        element_vec.push_back(Y(row_idx,i));
        sum += element_vec[count];
        count++;
    }
    mean = sum / count;
    int check = 1;
    if (check == 0) {
        for (int i = 0; i < count; i++) {
            residual_vec.push_back(element_vec[i] - mean);
            squared_residual_sum += residual_vec[i] * residual_vec[i];
        }
    }
    else if (check == 1) {
        for (int i = 0; i < count; i++) {
            residual_vec.push_back(element_vec[i] - mean);
        }
    }
    else {
        std::cout << "Check your data is snp or pheno" << std::endl;
        exit(1);
    }
}

input_MS::input_MS(double** Y, int row_idx, int col) {  //check 0 = snp, check 1 = pheno
    mean = 0.0; sum = 0.0; count = 0; squared_residual_sum = 0.0; slope = 0.0; intercept = 0.0;
    for (int i = 0; i < col; i++) {
        element_vec.push_back(Y[row_idx][i]);
        sum += element_vec[count];
        count++;
    }
    mean = sum / count;
    int check = 1;
    if (check == 0) {
        for (int i = 0; i < count; i++) {
            residual_vec.push_back(element_vec[i] - mean);
            squared_residual_sum += residual_vec[i] * residual_vec[i];
        }
    }
    else if (check == 1) {
        for (int i = 0; i < count; i++) {
            residual_vec.push_back(element_vec[i] - mean);
        }
    }
    else {
        std::cout << "Check your data is snp or pheno" << std::endl;
        exit(1);
    }
}

double input_MS::cal_estimate(input_MS* pheno) { //only use snp instance
    double sr_pr_multiply_sum = 0.0;

    if (this->count != pheno->count) {
        std::cout << "snp is " << this->count << std::endl;
        std::cout << "pheno is " << pheno->count << std::endl;
        std::cout << "snp, pheno count is not same pleas check" << std::endl;
        exit(1);
    }
    for (int i = 0; i < this->count; i++) {
        sr_pr_multiply_sum += this->residual_vec[i] * pheno->residual_vec[i];
    }
    this->slope = sr_pr_multiply_sum / this->squared_residual_sum;
    this->intercept = pheno->mean - this->slope * this->mean;
    //std::cout << "intercept is " << this->intercept << std::endl;
    return this->slope;
}

double input_MS::cal_RSE(input_MS* pheno) {
    double y_residual_sum = 0.0;

    for (int i = 0; i < this->count; i++) {
        y_residual_sum += std::pow(pheno->element_vec[i] - ((this->slope * this->element_vec[i]) + this->intercept), (double)2);
    }
    return std::sqrt((y_residual_sum / (this->count - 2)) / this->squared_residual_sum);
}

long double t_value(boost::math::students_t_distribution<double> t_dist, long double x) {
    double probability = 0.0;

    probability = std::abs(boost::math::quantile(t_dist, x));
    return probability;
}



double zscore(double x, double t_value) {
    if (t_value >= 0) {
        return std::abs(boost::math::quantile(boost::math::normal_distribution<double>(0, 1), x));
    }
    else {
        return -1 * std::abs(boost::math::quantile(boost::math::normal_distribution<double>(0, 1), x)); //�����Ҷ��� ��������ִ°� ������ �ð����غ���.
    }
}

void inputMS(std::string snp_p, std::string pheno_p) {
    std::ifstream snp(snp_p); // row : SNP   col : indi || x
    std::ifstream pheno(pheno_p); // row : pheno or expression col : indi || y

    int line_num = 0; double beta = 0.0; double zsc = 0.0;
    int pheno_num = 0;
    std::string temp = "";
    std::ofstream inputMS("inputMS_cpp.txt"); inputMS.precision(15);
    std::ofstream p_test("p_test_cpp.txt"); p_test.precision(15);

    int snp_check = 0;
    long double p_val = 0.0;
    if (!snp.is_open()) {
        std::cout << "snp file is not exist" << std::endl;
    }
    if (!pheno.is_open()) {
        std::cout << "pheno file is not exist" << std::endl;
    }

    while (!snp.eof()) { //snp�̶� pheno�� �Ųٷ� �ϸ� �ذ�!
        line_num++;
        std::vector<double> betas;
        std::vector<double> std_;
        getline(snp, temp);
        if (temp != "") {
            if (line_num > 1) {
                p_test << "\n";
                inputMS << "\n";
            }
            input_MS* snp_c = new input_MS(temp, 0);

            inputMS << line_num;

            pheno_num = 0;
            while (!pheno.eof()) {
                pheno_num++;
                getline(pheno, temp);
                if (temp != "") {
                    inputMS << " ";
                    if (pheno_num != 1) {
                        p_test << " ";
                    }
                    input_MS* pheno_c = new input_MS(temp, 1);
                    beta = snp_c->cal_estimate(pheno_c);
                    //std::cout << p_value(abc, snp_c->cal_estimate(pheno_c) / snp_c->cal_RSE(pheno_c)) << std::endl; //  p-value
                    //std::cout << std::endl;
                    //std::cout << "beta =" << beta <<std::endl;
                    //std::cout << "p value =" << p_value(boost::math::students_t_distribution<double>(snp_c->count - 2), snp_c->cal_estimate(pheno_c) / snp_c->cal_RSE(pheno_c)) * 2 << std::endl;
                    p_val = p_value(boost::math::students_t_distribution<double>(snp_c->count - 2), snp_c->cal_estimate(pheno_c) / snp_c->cal_RSE(pheno_c)) * 2;
                    zsc = zscore(p_val / 2, snp_c->cal_estimate(pheno_c) / snp_c->cal_RSE(pheno_c));

                    //std::cout.precision(15);
                    //std::cout << "Estimate =" << snp_c->cal_estimate(pheno_c) << std::endl; // coeff
                    //std::cout << "Std. Error =" << snp_c->cal_RSE(pheno_c) << std::endl; // Std. Error
                    //std::cout << "t value =" << snp_c->cal_estimate(pheno_c) / snp_c->cal_RSE(pheno_c) << std::endl; //t-value
                    //std::cout << "qnorm = " << std::abs(boost::math::quantile(boost::math::normal_distribution<double>(0, 1), p_value(boost::math::students_t_distribution<double>(snp_c->count - 2), snp_c->cal_estimate(pheno_c) / snp_c->cal_RSE(pheno_c)))) << std::endl;
                    //std::cout << "beta = " << beta << std::endl;
                    //std::cout << "zscore =" << zsc << std::endl;
                    p_test << p_val;
                    inputMS << beta << " " << beta / zsc;
                    betas.push_back(beta); std_.push_back(1.0 / std::pow((beta / zsc), 2));
                    delete pheno_c;
                }
            }
            //computeMvaluesMCMC(betas, std_, 1000000, pheno_num, Y); //1000000 is samplenumber;
            pheno.clear();
            pheno.seekg(0, std::ios::beg);

            delete snp_c;
        }
    }
    p_test.close();
    inputMS.close();
    snp.close();
    pheno.close();
}

void inputMS2(std::string snp_p, std::string pheno_p) {
    std::ifstream snp(snp_p); // row : SNP   col : indi || x
    std::ifstream pheno(pheno_p); // row : pheno or expression col : indi || y

    clock_t start, end;
    double result;

    int line_num = 0; double beta = 0.0; double zsc = 0.0;
    int pheno_num = 0;
    std::string temp = "";

    int snp_check = 0;
    long double p_val = 0.0;
    if (!snp.is_open()) {
        std::cout << "snp file is not exist" << std::endl;
    }
    if (!pheno.is_open()) {
        std::cout << "pheno file is not exist" << std::endl;
    }
    
    int y_row = count_matrix_row(pheno); int y_col = count_matrix_col(pheno);
    Eigen::MatrixXd Y = read_mat(pheno, y_row, y_col);
    Eigen::ArrayXd P_val;
    pheno.close();


    while (!snp.eof()) { //snp�̶� pheno�� �Ųٷ� �ϸ� �ذ�!
        P_val = Eigen::ArrayXd(y_row);
        std::cout << line_num << std::endl;
        line_num++;
        std::vector<double> betas;
        std::vector<double> std_;
        getline(snp, temp);
        if (temp != "") {
            input_MS* snp_c = new input_MS(temp, 0);


            pheno_num = 0;
            for(int i = 0; i < y_row; i++){
                pheno_num++;
                   
                    input_MS* pheno_c = new input_MS(Y, i, y_col);
                    beta = snp_c->cal_estimate(pheno_c);
                    p_val = p_value(boost::math::students_t_distribution<double>(snp_c->count - 2), snp_c->cal_estimate(pheno_c) / snp_c->cal_RSE(pheno_c)) * 2;
                    zsc = zscore(p_val / 2, snp_c->cal_estimate(pheno_c) / snp_c->cal_RSE(pheno_c));

                    P_val(i) = p_val;
                    betas.push_back(beta); std_.push_back(1.0 / std::pow((beta / zsc), 2));
                    delete pheno_c;
                
            }


            computeMvaluesMCMC(betas, std_, 1000000, pheno_num, temp, Y, P_val); //1000000 is samplenumber;



            delete snp_c;
        }
    }
    snp.close();
}

void inputMS3(std::string snp_p, std::string pheno_p) {
    std::ifstream snp(snp_p); // row : SNP   col : indi || x
    std::ifstream pheno(pheno_p); // row : pheno or expression col : indi || y

    int line_num = 0; double beta = 0.0; double zsc = 0.0;
    int pheno_num = 0;
    std::string temp = "";

    int snp_check = 0;
    long double p_val = 0.0;
    if (!snp.is_open()) {
        std::cout << "snp file is not exist" << std::endl;
    }
    if (!pheno.is_open()) {
        std::cout << "pheno file is not exist" << std::endl;
    }

    int y_row = count_matrix_row(pheno); int y_col = count_matrix_col(pheno);
    Eigen::MatrixXd Y = read_mat(pheno, y_row, y_col);
    Eigen::ArrayXd P_val;
    pheno.close();

    while (!snp.eof()) { //snp�̶� pheno�� �Ųٷ� �ϸ� �ذ�!
        P_val = Eigen::ArrayXd(y_row);
        std::cout << line_num << std::endl;
        line_num++;
        double* betas = new double[y_row];
        double* std_ = new double[y_row];
        double* std_tm = new double[y_row];
        double* std_tmm = new double[y_row];
        double* std_logt = new double[y_row];
        double* logProbNullPoints_ = new double[y_row];
        getline(snp, temp);
        if (temp != "") {
            input_MS* snp_c = new input_MS(temp, 0);
            for (int i = 0; i < y_row; i++) {
                input_MS* pheno_c = new input_MS(Y, i, y_col);
                beta = snp_c->cal_estimate(pheno_c);
                p_val = p_value(boost::math::students_t_distribution<double>(snp_c->count - 2), snp_c->cal_estimate(pheno_c) / snp_c->cal_RSE(pheno_c)) * 2;
                zsc = zscore(p_val / 2, snp_c->cal_estimate(pheno_c) / snp_c->cal_RSE(pheno_c));
                P_val(i) = p_val;
                betas[i] = beta; std_[i] = 1.0 / std::pow((beta / zsc), 2);
                std_tm[i] = std_[i] * beta;
                std_tmm[i] = std_tm[i] * beta;
                std_logt[i] = std::log(std_[i]);
                logProbNullPoints_[i] = std_logt[i] - std_[i] * beta * beta / 2;
                delete pheno_c;
            }
            //std::log(std_[i]) - std_[i] * betas[i] * betas[i] / 2;
            //,double* std_tm, double* std_tmm, double* std_logt,

            computeMvaluesMCMC2(betas, std_,std_tm, std_tmm, std_logt, logProbNullPoints_, 1000000, y_row, temp, Y, P_val); //1000000 is samplenumber;


            delete snp_c;
        }
        delete[] betas;
        delete[] std_;
        delete[] std_tm;
        delete[] std_tmm;
        delete[] std_logt;
        delete[] logProbNullPoints_;

    }
    snp.close();
}


void computeMvaluesMCMC(std::vector<double>& betas, std::vector<double>& std_, int sample, int pheno_num, std::string X_line, Eigen::MatrixXd& Y,
                                                                                                                                            Eigen::ArrayXd& P_val, int seed) {
    int maxNumFlip = 1;

    maxNumFlip = (std::floor(0.1 * pheno_num) > 1) ? std::floor(0.1 * pheno_num) : 1;
    //std -> ts �� ��ȯ�ϴ� �۾��� �ʿ�. �Ʒ� �ڵ�� �Ǿ��ִٰ� ���� ������ ��ġ�ų� ���⼭ ����
    std::vector<int> H1; //true false;
    std::vector<int> tmp;
    std::vector<int> shuffleBuffer;

    std::vector<double> logPriorConfig;
    std::vector<int> accumCntH0; std::vector<int> accumCntH1;

    ran.seed(seed);

    int numH1 = 0;
    long int burninCount = 1000;
    long int chainCount = 0;
    for (int i = 0; i < pheno_num; i++) {
        logPriorConfig.push_back(boost::math::beta(i + 1.0, pheno_num - (double)i + 1.0) - boost::math::beta(1.0, 1.0));
        accumCntH0.push_back(0); accumCntH1.push_back(0);

        tmp.push_back(0);
        H1.push_back(makeRandomInteger(1)); numH1 += H1[i];
        shuffleBuffer.push_back(i);
    }

    while (chainCount < sample) {
        double currentLogProb = observationLogLikelihood(betas, std_, H1, numH1) + logPriorConfig[numH1];
        if (makeRandomDouble() > 0.01) {
            // Usual jump
            int numFlip = (makeRandomInteger(maxNumFlip) + 1 > pheno_num) ? pheno_num : makeRandomInteger(maxNumFlip) + 1;

            for (int i = 0; i < numFlip; i++) {
                int pick = makeRandomInteger(pheno_num - i - 1);
                //   std::cout << "test = " << i + pick << std::endl;
                int t = shuffleBuffer[i];
                shuffleBuffer[i] = shuffleBuffer[i + pick];
                shuffleBuffer[i + pick] = t;

            }
            for (int i = 0; i < numFlip; i++) {
                int j = shuffleBuffer[i];
                H1[j] = !H1[j];
                numH1 += H1[j] ? 1 : -1;
            }

            double nextLogProb = observationLogLikelihood(betas, std_, H1, numH1) + logPriorConfig[numH1];

            if (nextLogProb > currentLogProb || makeRandomDouble() < std::exp(nextLogProb - currentLogProb)) {
                // Move
                currentLogProb = nextLogProb;
            }
            else {
                // Stay ... revert back
                for (int i = 0; i < numFlip; i++) {
                    int j = shuffleBuffer[i];
                    H1[j] = !H1[j];
                    numH1 += H1[j] ? 1 : -1;
                }
            }
        }
        else {

            // Randomization move
            int tmpNumH1 = 0;
            for (int i = 0; i < pheno_num; i++) {
                tmp[i] = makeRandomInteger(1);
                if (tmp[i]) tmpNumH1++;
            }
            double nextLogProb = observationLogLikelihood(betas, std_, tmp, tmpNumH1) + logPriorConfig[tmpNumH1];
            if (nextLogProb > currentLogProb || makeRandomDouble() < std::exp(nextLogProb - currentLogProb)) {
                // Move 
                for (int i = 0; i < pheno_num; i++) {
                    H1[i] = tmp[i];
                }
                numH1 = tmpNumH1;
                currentLogProb = nextLogProb;
            }
            else {
                // Stay ...
            }
        }

        // Are we still in Burn-in?
        if (burninCount > 0) {
            burninCount--;
        }
        else {
            for (int i = 0; i < pheno_num; i++) {
                if (H1[i]) {
                    accumCntH1[i]++;
                }
                else {
                    accumCntH0[i]++;
                }
            }
            chainCount++;
        }
    }


    int k_size = 0;
    Eigen::MatrixXd K = Eigen::MatrixXd(pheno_num, Y.cols());
    for (int i = 0; i < pheno_num; i++) {
        double mvalue = (double)accumCntH1[i] / (accumCntH0[i] + accumCntH1[i]);
        //mvalues_.add(mvalue);
        if (mvalue < 0.5) {
            K.row(k_size) = Y.row(i);
            k_size++;
        }
    }

    if (k_size >= 10) {
        K.resize(k_size, Y.cols());
        K = normMe(K);
        Eigen::MatrixXd cov_mat = cov(K);
        Eigen::MatrixXd X = read_mat(X_line, Y.cols());
        emma(X, Y, cov_mat, NICE);
    }
    else {
        for (int i = 0; i < Y.rows(); NICE << " ", i++) {
            NICE << P_val(i);
        }
        NICE << "\n";
    }
}

void computeMvaluesMCMC2(double* betas, double* std_,double* std_tm, double* std_tmm, double* std_logt, double* logProbNullPoints_, int sample, int pheno_num, std::string X_line, Eigen::MatrixXd& Y,
    Eigen::ArrayXd& P_val, int seed) {

    clock_t start, end;
    double result;


    int maxNumFlip = 1;
    maxNumFlip = (std::floor(0.1 * pheno_num) > 1) ? std::floor(0.1 * pheno_num) : 1;
    //std -> ts �� ��ȯ�ϴ� �۾��� �ʿ�. �Ʒ� �ڵ�� �Ǿ��ִٰ� ���� ������ ��ġ�ų� ���⼭ ����
    int* H1 = new int[pheno_num]; //true false;
    int* tmp = new int[pheno_num];
    int* shuffleBuffer = new int[pheno_num];

    double* logPriorConfig = new double[pheno_num];;
    int* accumCntH0 = new int[pheno_num]; int* accumCntH1 = new int[pheno_num];

    ran.seed(seed);

    int numH1 = 0;
    long int burninCount = 1000;
    long int chainCount = 0;
    for (int i = 0; i < pheno_num; i++) {
        logPriorConfig[i] = boost::math::beta(i + 1.0, pheno_num - (double)i + 1.0) - boost::math::beta(1.0, 1.0);
        accumCntH0[i]= 0; accumCntH1[i] = 0;

        tmp[i] = 0;
        H1[i] = makeRandomInteger(1); numH1 += H1[i];
        shuffleBuffer[i] = i;
    }


    start = clock();


    while (chainCount < sample) {
        double currentLogProb = observationLogLikelihood(betas, std_,std_tm, std_tmm, std_logt, logProbNullPoints_, H1, numH1, pheno_num) + logPriorConfig[numH1];
        if (makeRandomDouble() > 0.01) {
            // Usual jump
            int numFlip = (makeRandomInteger(maxNumFlip) + 1 > pheno_num) ? pheno_num : makeRandomInteger(maxNumFlip) + 1;

            for (int i = 0; i < numFlip; i++) {
                int pick = makeRandomInteger(pheno_num - i - 1);
                //   std::cout << "test = " << i + pick << std::endl;
                int t = shuffleBuffer[i];
                shuffleBuffer[i] = shuffleBuffer[i + pick];
                shuffleBuffer[i + pick] = t;

                int j = shuffleBuffer[i];
                H1[j] = !H1[j];
                numH1 += H1[j] ? 1 : -1;
            }
            //for (int i = 0; i < numFlip; i++) {
            //    int j = shuffleBuffer[i];
            //    H1[j] = !H1[j];
            //    numH1 += H1[j] ? 1 : -1;
            //}

            double nextLogProb = observationLogLikelihood(betas, std_, std_tm, std_tmm, std_logt, logProbNullPoints_, H1, numH1, pheno_num) + logPriorConfig[numH1];

            if (nextLogProb > currentLogProb || makeRandomDouble() < std::exp(nextLogProb - currentLogProb)) {
                // Move
                currentLogProb = nextLogProb;
            }
            else {
                // Stay ... revert back
                for (int i = 0; i < numFlip; i++) {
                    int j = shuffleBuffer[i];
                    H1[j] = !H1[j];
                    numH1 += H1[j] ? 1 : -1;
                }
            }
        }
        else {
            // Randomization move
            int tmpNumH1 = 0;
            for (int i = 0; i < pheno_num; i++) {
                tmp[i] = makeRandomInteger(1);
                if (tmp[i]) tmpNumH1++;
            }
            double nextLogProb = observationLogLikelihood(betas, std_, std_tm, std_tmm, std_logt, logProbNullPoints_, tmp, tmpNumH1, pheno_num) + logPriorConfig[tmpNumH1];
            if (nextLogProb > currentLogProb || makeRandomDouble() < std::exp(nextLogProb - currentLogProb)) {
                // Move 
                for (int i = 0; i < pheno_num; i++) {
                    H1[i] = tmp[i];
                }
                numH1 = tmpNumH1;
                currentLogProb = nextLogProb;
            }
            else {
                // Stay ...
            }
        }

        // Are we still in Burn-in?
        if (burninCount > 0) {
            burninCount--;
        }
        else {
            for (int i = 0; i < pheno_num; i++) {
                if (H1[i]) {
                    accumCntH1[i]++;
                }
                else {
                    accumCntH0[i]++;
                }
            }
            chainCount++;
        }

    }

    end = clock();

    result = (double)(end - start);

    std::cout << "working time: " << ((result) / CLOCKS_PER_SEC) << " seconds" << std::endl;


    int k_size = 0;
    int y_col = Y.cols();
    std::vector<long double*> K;
    // Eigen::MatrixXd K = Eigen::MatrixXd(pheno_num, Y.cols());
    for (int i = 0; i < pheno_num; i++) {
        double mvalue = (double)accumCntH1[i] / (accumCntH0[i] + accumCntH1[i]);
        //mvalues_.add(mvalue);
        if (mvalue < 0.5) {
            long double sum = 0.0;
            long double mean = 0.0;
            long double deviation = 0.0;
            long double* temp = new long double[y_col];
            for (int j = 0; j < y_col; j++) {
                sum += Y(i,j);
            }
            mean = sum / y_col;

            for (int j = 0; j < y_col; j++) {
                double temp = Y(i, j) - mean;
                deviation += temp * temp;
            }
            for (int j = 0; j < y_col; j++) {
                temp[j] = (Y(i,j) - mean) / std::sqrt(deviation / (y_col - 1));
            }
            K.push_back(temp);
        }
    }


    if (K.size() >= 10) {
        //Eigen::MatrixXd cov_mat = cov(K);
        //emma(X, Y, cov_mat, NICE);
        Eigen::MatrixXd Ki = Eigen::MatrixXd(y_col, y_col);
        for (int i = 0; i < y_col; i++) {
            for (int j = i; j < y_col; j++) {
                double cov_ij = 0.0;
                double sum_i = 0.0; double sum_j = 0.0;
                double mean_i = 0.0; double mean_j = 0.0;
                for (int k = 0; k < K.size(); k++) {
                    sum_i += K[k][i]; sum_j = K[k][j];
                }
                mean_i = sum_i / pheno_num; mean_j = sum_j / pheno_num;
                for (int k = 0; k < K.size(); k++) {
                    cov_ij += (K[k][i]-mean_i) * (K[k][j]-mean_j);
                }
                Ki(i, j) = cov_ij/(pheno_num-1);
                Ki(j, i) = Ki(i, j);
            }
        }
        Eigen::MatrixXd X = read_mat(X_line, Y.cols());
  

        emma(X, Y, Ki, NICE);
    }
    else {
        for (int i = 0; i < Y.rows(); NICE << " ", i++) {
            NICE << P_val(i);
        }
        NICE << "\n";
    }

    delete[] H1;
    delete[] tmp;
    delete[] shuffleBuffer;
    delete[] logPriorConfig;
    delete[] accumCntH0;
    delete[] accumCntH1;

}

double observationLogLikelihood(std::vector<double>& betas, std::vector<double>& std_, std::vector<int>& H1, int numH1) {
    int n = betas.size();
    double logProbNullPoints = 0;
    double logProbAltPoints = 0;

    double log_sqrt2pi = 0.5 * std::log(2 * std::atan(1) * 4);

    for (int i = 0; i < n; i++) {
        if (!H1[i]) logProbNullPoints += std::log(std_[i]) - std_[i] * betas[i] * betas[i] / 2;
    }
    if (numH1 > 0) {
        double sum_t = 0.0;
        double sum_tm = 0.0;
        double sum_tmm = 0.0;
        double sum_logt = 0.0;
        for (int i = 0; i < n; i++) {
            if (H1[i]) {
                sum_t += std_[i];
                sum_tm += std_[i] * betas[i];
                sum_tmm += std_[i] * betas[i] * betas[i];
                sum_logt += std::log(std_[i]);
            }
        }
        double betaJoint = sum_tm / sum_t;
        double tJoint = sum_t;
        double tFinal = 1.0 / ((1.0 / tJoint) + 0.04);
        double logScaleFactor =
            -(numH1 - 1) * log_sqrt2pi + 0.5 * sum_logt - 0.5 * std::log(sum_t)
            - (sum_tmm - sum_tm * sum_tm / sum_t) / 2;
        double logJointPDF =
            0.5 * std::log(tFinal) - log_sqrt2pi - tFinal * betaJoint * betaJoint / 2;
        logProbAltPoints = logJointPDF + logScaleFactor;
    }

    return logProbNullPoints + logProbAltPoints;
}

double observationLogLikelihood(double* betas, double* std_,double* std_tm, double* std_tmm, double* std_logt, double* logProbNullPoints_, int* H1, int numH1, int n_size) {
    int n = n_size;
    double logProbNullPoints = 0;
    double logProbAltPoints = 0;

    for (int i = 0; i < n; i++) {
        if (!H1[i]) logProbNullPoints += logProbNullPoints_[i];
            //std::log(std_[i]) - std_[i] * betas[i] * betas[i] / 2;
    }
    if (numH1 > 0) {
        double sum_t = 0.0;
        double sum_tm = 0.0;
        double sum_tmm = 0.0;
        double sum_logt = 0.0;
        for (int i = 0; i < n; i++) {
            if (H1[i]) {
                sum_t += std_[i];
                sum_tm += std_tm[i];
                sum_tmm += std_tmm[i];
                sum_logt += std_logt[i];
            }
        }
        double betaJoint = sum_tm / sum_t;
        double tJoint = sum_t;
        double tFinal = 1.0 / ((1.0 / tJoint) + 0.04);
        double logScaleFactor =
            -(numH1 - 1) * log_sqrt2pi + 0.5 * sum_logt - 0.5 * std::log(sum_t)
            - (sum_tmm - sum_tm * sum_tm / sum_t) / 2;
        double logJointPDF =
            0.5 * std::log(tFinal) - log_sqrt2pi - tFinal * betaJoint * betaJoint / 2;
        logProbAltPoints = logJointPDF + logScaleFactor;
    }

    return logProbNullPoints + logProbAltPoints;
}


int makeRandomInteger(int number) {
    std::uniform_int_distribution<int> dis(0, number);
    return dis(ran);
}

double makeRandomDouble() {
    std::uniform_real_distribution<double> dis(0, 1);
    return dis(ran);
}

Eigen::MatrixXd read_mat(std::ifstream& input_file, int row, int col) {

    std::string read_buffer;
    std::string token;
    std::stringstream stream;
    Eigen::MatrixXd ret_mat = Eigen::MatrixXd(row, col);
    for (int i = 0; i < row; i++) { //row
        std::getline(input_file, read_buffer);
        stream.str(read_buffer);
        for (int j = 0; j < col; j++) { //col
            stream >> token;
            ret_mat(i, j) = std::stold(token);
        }
        stream.clear();  //for coursor reset
    }
    return ret_mat;
}

double** read_mat_darray(std::ifstream& input_file, int row, int col) {

    std::string read_buffer;
    std::string token;
    std::stringstream stream;
    double** ret_mat = new double*[row];
    for (int i = 0; i < row; i++) { //row
        ret_mat[i] = new double[col];
        std::getline(input_file, read_buffer);
        stream.str(read_buffer);
        for (int j = 0; j < col; j++) { //col
            stream >> token;
            ret_mat[i][j] = std::stold(token);
        }
        stream.clear();  //for coursor reset
    }
    return ret_mat;
}


Eigen::MatrixXd read_mat(std::string X, int col) {

    std::string token;
    std::stringstream stream;
    Eigen::MatrixXd ret_mat = Eigen::MatrixXd(1, col);
    stream.str(X);
    for (int j = 0; j < col; j++) { //col
        stream >> token;
        ret_mat(0, j) = std::stold(token);
    }
    stream.clear();  //for coursor reset

    return ret_mat;
}


int count_matrix_col(std::ifstream& matrix) {
    std::string token;
    std::stringstream stream; int count = 0;
    std::string read_line;

    std::getline(matrix, read_line);
    stream.str(read_line);
    while (stream >> token) {
        count++;
    }
    matrix.clear(); //coursor reset
    matrix.seekg(0, std::ios_base::beg);
    return count;
}

int count_matrix_row(std::ifstream& matrix) {
    std::string read_line; int count = 0;
    while (matrix.peek() != EOF) {
        std::getline(matrix, read_line);
        count++;
    }
    matrix.clear();//coursor reset
    matrix.seekg(0, std::ios_base::beg);
    return count;
}

void emma_eigen_L_wo_Z(Eigen::MatrixXd& Kinship, Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>& eig_L) {
    eig_L.compute(Kinship);
}

void emma_eigen_R_wo_Z(Eigen::MatrixXd& Kinship, Eigen::MatrixXd& X, Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>& eig_R, int indi) { //X�� 2�÷� ��Ʈ������ ù�÷��� ������ 1
    Eigen::MatrixXd unit_mat; unit_mat.setIdentity(indi, indi);
    Eigen::MatrixXd S = unit_mat - X * (X.transpose() * X).inverse() * X.transpose();
    Eigen::MatrixXd eig = S * (Kinship + unit_mat) * S;
    eig_R.compute(eig);
}

double emma_delta_REML_LL_wo_Z(double logdelta, Eigen::Map<Eigen::ArrayXd>& eig_value, Eigen::Map<Eigen::ArrayXd>& etas) {
    int et_size = etas.size(); double delta = std::exp(logdelta);

    return 0.5 * (et_size * (std::log((double)et_size / (2 * std::atan(1) * 4)) - 1 - std::log((etas.square() / (eig_value + delta)).sum())) - (eig_value + delta).log().sum());
}

double emma_delta_REML_dLL_wo_Z(double logdelta, Eigen::Map<Eigen::ArrayXd>& eig_value, Eigen::Map<Eigen::ArrayXd>& etas, Eigen::ArrayXd& etasq) { //etasq = etas * etas
    int et_size = etas.size(); double delta = std::exp(logdelta);
    Eigen::ArrayXd ldelta = eig_value + delta;
    return 0.5 * (et_size * (etasq / (ldelta * ldelta)).sum() / (etasq / ldelta).sum() - (1 / ldelta).sum());
}

double uniroot_emma_delta_REML_dLL_wo_Z(double a, double b, Eigen::Map<Eigen::ArrayXd>& eig_value, Eigen::Map<Eigen::ArrayXd>& etas, Eigen::ArrayXd& etasq) {
    double c;
    int m = 0;
    while (std::abs((b - a)) > 0.0001220703125) {
        //�ݺй����� �ϸ� ��갪�� �ٸ�.. �ٸ� ��� ���� �ʿ�
        m++;
        c = (a + b) / 2;
        if (emma_delta_REML_dLL_wo_Z(a, eig_value, etas, etasq) * emma_delta_REML_dLL_wo_Z(c, eig_value, etas, etasq) < 0 && m < 1000) {
            b = c;
        }
        else if (emma_delta_REML_dLL_wo_Z(c, eig_value, etas, etasq) * emma_delta_REML_dLL_wo_Z(b, eig_value, etas, etasq) < 0 && m < 1000) {
            a = c;
        }
        else {
            break;
        }
    }
    return (a + b) / 2;
}

void emma_REMLE(Eigen::ArrayXd& y, Eigen::MatrixXd& x, Eigen::MatrixXd& K, Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>& eig_R, REMLE& REM, int indi) {
    //y is y[i,] cols x is cbind x0
    Eigen::ArrayXd logdelta(101);
    logdelta << -10.0, -9.8, -9.6, -9.4, -9.2, -9.0, -8.8, -8.6, -8.4, -8.2, -8.0, -7.8, -7.6, -7.4, -7.2, -7.0, -6.8, -6.6, -6.4,
        -6.2, -6.0, -5.8, -5.6, -5.4, -5.2, -5.0, -4.8, -4.6, -4.4, -4.2, -4.0, -3.8, -3.6, -3.4, -3.2, -3.0, -2.8, -2.6,
        -2.4, -2.2, -2.0, -1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2,
        1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0,
        5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6, 8.8,
        9.0, 9.2, 9.4, 9.6, 9.8, 10.0;
    Eigen::ArrayXd exp_logdelta(101);
    exp_logdelta << 4.539992976248485417315e-05, 5.545159943217694512724e-05, 6.772873649085389821320e-05, 8.272406555663222842839e-05, 1.010394018370934184937e-04, 1.234098040866795612143e-04, 1.507330750954765036589e-04, 1.841057936675791944017e-04,
        2.248673241788481946563e-04, 2.746535699721425367870e-04, 3.354626279025118532236e-04, 4.097349789797868077970e-04, 5.004514334406108327577e-04, 6.112527611295723002985e-04, 7.465858083766798531863e-04, 9.118819655545162437479e-04,
        1.113775147844803248276e-03, 1.360368037547893867861e-03, 1.661557273173933915386e-03, 2.029430636295734020208e-03, 2.478752176666358490731e-03, 3.027554745375815332259e-03, 3.697863716482932029683e-03, 4.516580942612670283853e-03,
        5.516564420760771553232e-03, 6.737946999085467000845e-03, 8.229747049020030152944e-03, 1.005183574463358597839e-02, 1.227733990306844810703e-02, 1.499557682047770318379e-02, 1.831563888873417866865e-02, 2.237077185616560132120e-02,
        2.732372244729256924312e-02, 3.337326996032610043619e-02, 4.076220397836624598220e-02, 4.978706836786394446248e-02, 6.081006262521792410380e-02, 7.427357821433390461241e-02, 9.071795328941247016363e-02, 1.108031583623339672018e-01,
        1.353352832366127023178e-01, 1.652988882215864208103e-01, 2.018965179946554666657e-01, 2.465969639416063785564e-01, 3.011942119122023031608e-01, 3.678794411714423340243e-01, 4.493289641172220627574e-01, 5.488116360940256122092e-01,
        6.703200460356391054972e-01, 8.187307530779824871203e-01, 1.000000000000000000000e+00, 1.221402758160168966484e+00, 1.491824697641270791593e+00, 1.822118800390511550447e+00, 2.225540928492469205935e+00, 2.718281828459045090796e+00,
        3.320116922736551234863e+00, 4.055199966844669212662e+00, 4.953032424395113153537e+00, 6.049647464412939434908e+00, 7.389056098930650406942e+00, 9.025013499434114905284e+00, 1.102317638064160476574e+01, 1.346373803500168619962e+01,
        1.644464677109706229885e+01, 2.008553692318766792368e+01, 2.453253019710937365971e+01, 2.996410004739702515053e+01, 3.659823444367803801924e+01, 4.470118449330077226023e+01, 5.459815003314423620395e+01, 6.668633104092509711336e+01,
        8.145086866496799871129e+01, 9.948431564193377596439e+01, 1.215104175187349682119e+02, 1.484131591025765999348e+02, 1.812722418751510531365e+02, 2.214064162041871668407e+02, 2.704264074261529913201e+02, 3.302995599096489058866e+02,
        4.034287934927351102488e+02, 4.927490410932576310188e+02, 6.018450378720812068423e+02, 7.350951892419712976334e+02, 8.978472916504183558573e+02, 1.096633158428458500566e+03, 1.339430764394416883079e+03, 1.635984429995924301693e+03,
        1.998195895104120836550e+03, 2.440601977624500705133e+03, 2.980957987041728301847e+03, 3.640950307332352167577e+03, 4.447066747699865118193e+03, 5.431659591362988066976e+03, 6.634244006277865992161e+03, 8.103083927575384223019e+03,
        9.897129058743908899487e+03, 1.208838073021696800424e+04, 1.476478156557729380438e+04, 1.803374492782852394157e+04, 2.202646579480671789497e+04;

    if ((x.transpose() * x).determinant() == 0) {
        return;
    }
    int m = 101;
    Eigen::ArrayXd eig_v = eig_R.eigenvalues().reverse().array() - 1;
    Eigen::Map<Eigen::ArrayXd> eigenvalues(eig_v.data(), indi - 2);
    Eigen::ArrayXd etas_100 = (eig_R.eigenvectors().transpose() * y.matrix()).reverse();
    Eigen::Map<Eigen::ArrayXd> etas_98(etas_100.data(), indi - 2);
    Eigen::MatrixXd Lambdas = Eigen::MatrixXd(indi - 2, m);
    for (int i = 0; i < m; i++) {
        Lambdas.col(i) = eigenvalues + exp_logdelta(i);
    }
    Eigen::ArrayXd Etasq = etas_98.square(); //��Ʈ���� �ʿ䰡����.
    Eigen::MatrixXd Eta_Lambdas = Eigen::MatrixXd(indi - 2, m);
    Eigen::MatrixXd Eta_Lambdas_sqr = Eigen::MatrixXd(indi - 2, m);
    for (int i = 0; i < m; i++) {
        Eta_Lambdas.col(i) = Etasq / Lambdas.col(i).array();
        Eta_Lambdas_sqr.col(i) = Eta_Lambdas.col(i).array() / Lambdas.col(i).array();
    }

    Eigen::ArrayXd LL = 0.5 * ((indi - 2) * (std::log((double)(indi - 2) / (2 * std::atan(1) * 4)) - 1 - (Eta_Lambdas.colwise().sum()).array().log()) - Lambdas.array().log().colwise().sum());
    Eigen::ArrayXd dLL = 0.5 * exp_logdelta.transpose() * ((indi - 2) * Eta_Lambdas_sqr.colwise().sum().array() / Eta_Lambdas.colwise().sum().array() - Eigen::MatrixXd(1 / Lambdas.array()).colwise().sum().array());

    std::vector<double> optlogdelta;
    std::vector<double> optLL;

    //std::cout << eigenvalues << std::endl;
    //std::cout << emma_delta_REML_LL_wo_Z(-10, eigenvalues, etas_98) << std::endl;
    if (dLL(0) < 1e-10) {
        optlogdelta.push_back(-10.0);
        optLL.push_back(emma_delta_REML_LL_wo_Z(-10, eigenvalues, etas_98));
    }
    if (dLL(99) > 0 - (1e-10)) {
        optlogdelta.push_back(10.0);
        optLL.push_back(emma_delta_REML_LL_wo_Z(10, eigenvalues, etas_98));
    }
    for (int i = 0; i < 100; i++) {
        if ((dLL(i) * dLL(i + 1) < 0 - std::pow(1e-10, 2)) && (dLL(i) > 0) && (dLL(i + 1) < 0)) {
            double root = uniroot_emma_delta_REML_dLL_wo_Z(logdelta[i], logdelta[i + 1], eigenvalues, etas_98, Etasq);
            optlogdelta.push_back(root);
            optLL.push_back(emma_delta_REML_LL_wo_Z(root, eigenvalues, etas_98));

            //std::cout <<"root is "<< root << " i is" << i << std::endl;
        }
    }
    if (optLL.size() == 0) {
        REM.delta = 0; REM.REML = 0; REM.ve = 0; REM.vg = 0;
        return;
    }
    int optLL_maxidx = 0;
    for (int i = 0; i < optLL.size() - 1; i++) {
        if (optLL[i] < optLL[i + 1]) {
            optLL_maxidx = i + 1;
        }
    }

    REM.delta = std::exp(optlogdelta[optLL_maxidx]);
    REM.REML = optLL[optLL_maxidx];
    REM.vg = (Etasq / (eigenvalues + REM.delta)).sum() / (indi - 2);
    REM.ve = REM.vg * REM.delta;
}

Eigen::MatrixXd cbind(Eigen::MatrixXd& a, Eigen::MatrixXd& b) {
    int a_rows = a.rows();
    int b_rows = b.rows();
    int b_cols = b.cols();
    if (a_rows != b_rows) {
        std::cout << "can't bind different row matrix" << std::endl;
        exit(1);
    }
    Eigen::MatrixXd ret_mat = Eigen::MatrixXd(a_rows, b_cols + 1);
    ret_mat.col(0) = a.col(0);
    for (int i = 1; i < b_cols + 1; i++) {
        ret_mat.col(i) = b.col(i - 1);
    }
    return ret_mat;
}

long double p_value(boost::math::students_t_distribution<double> t_dist, long double x) {

    double x2 = x * x;
    double probability;
    if (t_dist.degrees_of_freedom() > 2 * x2)
    {
        double z = x2 / (t_dist.degrees_of_freedom() + x2);
        probability = boost::math::ibetac(static_cast<double>(0.5), t_dist.degrees_of_freedom() / 2, z) / 2;
    }
    else
    {
        double z = t_dist.degrees_of_freedom() / (t_dist.degrees_of_freedom() + x2);
        probability = boost::math::ibeta(t_dist.degrees_of_freedom() / 2, static_cast<double>(0.5), z) / 2;
    }

    return (probability);
    //return (1 - (long double)boost::math::cdf(t_dist, x));
    //return boost::math::pvalue(t_dist, x);
}

void emma(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, Eigen::MatrixXd& K, std::ofstream& out) {
    REMLE rem = REMLE();

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig_L;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig_R1;

    Eigen::MatrixXd X_0 = Eigen::MatrixXd(100, 1); X_0.setOnes();
    X.transposeInPlace();
    X_0 = cbind(X_0, X);

    emma_eigen_L_wo_Z(K, eig_L);
    emma_eigen_R_wo_Z(K, X_0, eig_R1, 100);
    //std::cout << (X_0.transpose() * X_0).determinant() << std::endl;

    /////start for 

    for (int i = 0; i < Y.rows(); out << " ", i++) {
        Eigen::ArrayXd Y_row = Y.row(i);
        emma_REMLE(Y_row, X_0, K, eig_R1, rem, 100);

        Eigen::ArrayXd U_sqrt_egv = (1 / (eig_L.eigenvalues().reverse().array() + rem.delta)).sqrt().array();
        Eigen::MatrixXd U = Eigen::MatrixXd(100, 100);// ���߿� �����ϳ� ����
        for (int j = 0; j < 100; j++) {
            U.col(j) = eig_L.eigenvectors().col(99 - j) * U_sqrt_egv(j);
        }
        Eigen::ArrayXd yt = U.transpose() * Y_row.matrix();
        Eigen::MatrixXd xt = U.transpose() * X_0;
        Eigen::MatrixXd	iXX = (xt.transpose() * xt).inverse();
        Eigen::ArrayXd beta = iXX * xt.transpose() * yt.matrix();

        double stats = (beta(1) / std::sqrt(iXX(1, 1) * rem.vg));

        out << 2 * p_value(boost::math::students_t_distribution<double>(98), stats);
    }
    out << "\n";
}



Eigen::MatrixXd cov(Eigen::MatrixXd& mat) {
    Eigen::MatrixXd centered = mat.rowwise() - mat.colwise().mean();
    return (centered.adjoint() * centered) / double(mat.rows() - 1);
}

Eigen::MatrixXd normMe(Eigen::MatrixXd& mat) {
    Eigen::ArrayXd std = (((mat.colwise() - mat.rowwise().mean()).rowwise().squaredNorm()) / (mat.cols() - 1)).cwiseSqrt();
    mat = (mat.colwise() - mat.rowwise().mean());
    for (int i = 0; i < std.size(); i++) {
        mat.row(i) = mat.row(i) / std(i);
    }
    return mat;
}
