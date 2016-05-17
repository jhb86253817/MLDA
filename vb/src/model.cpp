#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "model.h"
#include "dataset.h"
#include "utils.h"
#include "constants.h"
#include <boost/math/special_functions/digamma.hpp> // digamma
#include <boost/math/special_functions/gamma.hpp> // lgamma
#include <math.h>

using namespace std;

void model::set_default_values() {
    wordmapfile = "wordmap.txt";
    wordmapfile2 = "wordmap2.txt";
    theta_suffix = ".theta";
    /////////////////////////////////////////////////////////////////////////////////////
    beta_suffix = ".phi"; // to be more convenient to do inference using Gibbs
    others_suffix = ".others";
    twords_suffix = ".twords";

    dir = "./";
    dir2 = "./";
    dfile = "";
    dfile2 = "";
    train_method = 1;
    model_name = "model-final-1";
    model_name2 = "model-final-2";
    ptrndata = NULL;
    ptrndata2 = NULL;
    M = 0;
    M2 = 0;
    V = 0;
    V2 = 0;
    K = 50;
    alpha = 50.0 / K;
    eta = 0.1;
    niters = 40;
    niters2 = 10;
    twords = 20;
    nested_iters = 2;

    train_time_lang1 = 0.0;
    train_time_lang2 = 0.0;

    phi = NULL;
    phi2 = NULL;
    gamma = NULL;
    gamma_sum = NULL;
    lambda = NULL;
    lambda2 = NULL;
    lambda_sum = NULL;
    lambda_sum2 = NULL;
    dig_gamma = NULL;
    dig_lambda = NULL;
    dig_lambda2 = NULL;
    dig_lambda_sum = NULL;
    dig_lambda_sum2 = NULL;

    theta = NULL;
    beta = NULL;
    beta2 = NULL;
}

int model::parse_args(int argc, char ** argv) {
    return utils::parse_args(argc, argv, this);
}

int model::init(int argc, char ** argv) {

    srandom(1); // set random seed

    if (parse_args(argc, argv)) {
        return 1;
    }
    if (init_est()) {
        return 1;
    }

    return 0;
}

int model::init_est() {
    int m, n, k, v;
    int m2, n2, k2, v2;

    // read training data
    ptrndata = new dataset;
    if (ptrndata->read_trndata(dir + dfile, dir + wordmapfile)) {
        printf("Failed to read training data!\n");
        return 1;
    }

    // allocate memory and assign values for variables
    M = ptrndata->M;
    V = ptrndata->V;

    // read training data 2
    ptrndata2 = new dataset;
    if (ptrndata2->read_trndata(dir2 + dfile2, dir2 + wordmapfile2)) {
        printf("Failed to read training data2!\n");
        return 1;
    }

    // allocate memory and assign values for variables
    M2 = ptrndata2->M;
    V2 = ptrndata2->V;

    ////////////////////////////////////////////////////////////
    // use corpus to initialize lambda and lambda2
    lambda = new double*[K];
    for (k = 0; k < K; k++) {
        lambda[k] = new double[V];
        for (v = 0; v < V; v++) {
            lambda[k][v] = eta;
        }
    }
    lambda2 = new double*[K];
    for (k2 = 0; k2 < K; k2++) {
        lambda2[k2] = new double[V2];
        for (v2 = 0; v2 < V2; v2++) {
            lambda2[k2][v2] = eta;
        }
    }
    for (k = 0; k < K; k++) {
        // just fix the seed doc number as one, which is better than two to five
        for (int i = 0; i < 4; i++) {
            int d = (int)(((double)random() / RAND_MAX) * M);
            int N = ptrndata->docs[d]->length;
            for (n = 0; n < N; n++) {
                int word = ptrndata->docs[d]->words[n];
                lambda[k][word] += 1;
            }
            int N2 = ptrndata2->docs[d]->length;
            for (n2 = 0; n2 < N2; n2++) {
                int word2 = ptrndata2->docs[d]->words[n2];
                lambda2[k][word2] += 1;
            }
        }
    }
    lambda_sum = new double[K];
    for (k = 0; k < K; k++) {
        lambda_sum[k] = 0.0;
        for (v = 0; v < V; v++) {
            lambda_sum[k] += lambda[k][v];
        }
    }
    lambda_sum2 = new double[K];
    for (k2 = 0; k2 < K; k2++) {
        lambda_sum2[k2] = 0.0;
        for (v2 = 0; v2 < V2; v2++) {
            lambda_sum2[k2] += lambda2[k2][v2];
        }
    }
    ////////////////////////////////////////////////////////////

    phi = new double*[ptrndata->N_max];
    for (n = 0; n < ptrndata->N_max; n++) {
        phi[n] = new double[K];
        for (k = 0; k < K; k++) {
            phi[n][k] = 1.0 / K;
        }
    }

    // initialize gamma 
    gamma = new double*[M];
    for (m = 0; m < M; m++) {
        gamma[m] = new double[K];
        int N = ptrndata->docs[m]->length;
        for (k = 0; k < K; k++) {
            gamma[m][k] = alpha + 1.0 * N / K;
        }
    }
    dig_gamma = new double*[M];
    for (m = 0; m < M; m++) {
        dig_gamma[m] = new double[K];
        for (k = 0; k < K; k++) {
            dig_gamma[m][k] = 0.0;
        }
    }
    gamma_sum = new double[M];
    for (m = 0; m < M; m++) {
        gamma_sum[m] = 0.0;
    }

    dig_lambda = new double*[K];
    for (k = 0; k < K; k++) {
        dig_lambda[k] = new double[V];
        for (v = 0; v < V; v++) {
            dig_lambda[k][v] = 0.0;
        }
    }
    dig_lambda_sum = new double[K];
    for (k = 0; k < K; k++) {
        dig_lambda_sum[k] = 0.0;
    }

    theta = new double*[M];
    for (m = 0; m < M; m++) {
        theta[m] = new double[K];
        for (k = 0; k < K; k++) {
            theta[m][k] = 0.0;
        }
    }

    beta = new double*[K];
    for (k = 0; k < K; k++) {
        beta[k] = new double[V];
        for (v = 0; v < V; v++) {
            beta[k][v] = 0.0;
        }
    }

    phi2 = new double*[ptrndata2->N_max];
    for (n2 = 0; n2 < ptrndata2->N_max; n2++) {
        phi2[n2] = new double[K];
        for (k2 = 0; k2 < K; k2++) {
            phi2[n2][k2] = 1.0 / K;
        }
    }

    dig_lambda2 = new double*[K];
    for (k2 = 0; k2 < K; k2++) {
        dig_lambda2[k2] = new double[V2];
        for (v2 = 0; v2 < V2; v2++) {
            dig_lambda2[k2][v2] = 0.0;
        }
    }
    dig_lambda_sum2 = new double[K];
    for (k2 = 0; k2 < K; k2++) {
        dig_lambda_sum2[k2] = 0.0;
    }

    beta2 = new double*[K];
    for (k2 = 0; k2 < K; k2++) {
        beta2[k2] = new double[V2];
        for (v2 = 0; v2 < V2; v2++) {
            beta2[k2][v2] = 0.0;
        }
    }

    return 0;
}

//int model::init_est_greedy() {
//    int m2, n2, k2;
//
//    // compute theta using gamma
//    for (m2 = 0; m2 < M2; m2++) {
//        for (k2 = 0; k2 < K; k2++) {
//            theta[m2][k2] = gamma[m2][k2] / gamma_sum[m2];
//        }
//    }
//
//    for (m2 = 0; m2 < M2; m2++) {
//        int N2 = ptrndata2->docs[m2]->length;
//        for (n2 = 0; n2 < N2; n2++) {
//            for (k2 = 0; k2 < K; k2++) {
//                phi2[m2][n2][k2] = theta[m2][k2];
//            }
//        }
//    }
//
//    return 0;
//}

void model::estimate() {
    int m, n, k, v;
    int m2, n2, k2, v2;

    // read wordmap file and map id to word
    if (twords > 0) {
        dataset::read_wordmap(dir + wordmapfile, &id2word);
        dataset::read_wordmap(dir2 + wordmapfile2, &id2word2);
    }

    if (train_method == 1) {
        printf("Use training method 1!\n");
        printf("Do %d iterations for 2 languages!\n", niters);

        // recording training time
        clock_t tStart = clock();

        // coordinate ascent for niters iterations
        for (int iter = 1; iter <= niters; iter++) {
            printf("Iteration %d...\n", iter);

            // compute digamma function of lambda and lambda_sum
            for (k = 0; k < K; k++) {
                for (v = 0; v < V; v++) {
                    dig_lambda[k][v] = exp( boost::math::digamma(lambda[k][v]) );
                }
                dig_lambda_sum[k] = exp( boost::math::digamma(lambda_sum[k]) );
            }
            // compute digamma function of lambda2 and lambda_sum2
            for (k2 = 0; k2 < K; k2++) {
                for (v2 = 0; v2 < V2; v2++) {
                    dig_lambda2[k2][v2] = exp( boost::math::digamma(lambda2[k2][v2]) );
                }
                dig_lambda_sum2[k2] = exp( boost::math::digamma(lambda_sum2[k2]) );
            }
            // compute digamma function of gamma 
            for (m = 0; m < M; m++) {
                for (k = 0; k < K; k++) {
                    dig_gamma[m][k] = exp( boost::math::digamma(gamma[m][k]) );
                }
            }

            // initialize lambda with eta
            for (k = 0; k < K; k++) {
                for (v = 0; v < V; v++) {
                    lambda[k][v] = eta;
                }
            }
            // initialize lambda2 with eta
            for (k2 = 0; k2 < K; k2++) {
                for (v2 = 0; v2 < V2; v2++) {
                    lambda2[k2][v2] = eta;
                }
            }
            // initialize gamma with alpha
            for (m = 0; m < M; m++) {
                for (k = 0; k < K; k++) {
                    gamma[m][k] = alpha;
                }
            }

            // update phi and add new phi to lambda and gamma statistics
            for (m = 0; m < M; m++) {
                int N = ptrndata->docs[m]->length;
                // update phi of a doc
                for (n = 0; n < N; n++) {
                    int word = ptrndata->docs[m]->words[n];
                    double phi_sum_temp = 0.0;
                    for (k = 0; k < K; k++) {
                        phi[n][k] =  dig_gamma[m][k] * dig_lambda[k][word] / dig_lambda_sum[k];
                        phi_sum_temp += phi[n][k];
                    }
                    //phi_sum[n] = phi_sum_temp;
                    // normalize phi
                    for (k = 0; k < K; k++) {
                        phi[n][k] /= phi_sum_temp;
                        // update lambda
                        lambda[k][word] += phi[n][k];
                        // update gamma
                        gamma[m][k] += phi[n][k];
                    }
                }
            }
            // update phi2 and add new phi2 to lambda2 and gamma statistics
            for (m2 = 0; m2 < M2; m2++) {
                int N2 = ptrndata2->docs[m2]->length;
                // update phi2 of a doc
                for (n2 = 0; n2 < N2; n2++) {
                    int word2 = ptrndata2->docs[m2]->words[n2];
                    double phi_sum_temp2 = 0.0;
                    for (k2 = 0; k2 < K; k2++) {
                        phi2[n2][k2] =  dig_gamma[m2][k2] * dig_lambda2[k2][word2] / dig_lambda_sum2[k2];
                        phi_sum_temp2 += phi2[n2][k2];
                    }
                    //phi_sum2[m2][n2] = phi_sum_temp2;
                    // normalize phi2
                    for (k2 = 0; k2 < K; k2++) {
                        phi2[n2][k2] /= phi_sum_temp2;
                        // update lambda2
                        lambda2[k2][word2] += phi2[n2][k2];
                        // update gamma
                        gamma[m2][k2] += phi2[n2][k2];
                    }
                }
            }

            // update lambda_sum
            for (k = 0; k < K; k++) {
                double lambda_sum_temp = 0.0;
                for (v = 0; v < V; v++) {
                    lambda_sum_temp += lambda[k][v];
                }
                lambda_sum[k] = lambda_sum_temp;
            }
            // update lambda_sum2
            for (k2 = 0; k2 < K; k2++) {
                double lambda_sum_temp2 = 0.0;
                for (v2 = 0; v2 < V2; v2++) {
                    lambda_sum_temp2 += lambda2[k2][v2];
                }
                lambda_sum2[k2] = lambda_sum_temp2;
            }

            // update gamma_sum
            for (m = 0; m < M; m++) {
                double gamma_sum_temp = 0.0;
                for ( k = 0; k < K; k++) {
                    gamma_sum_temp += gamma[m][k];
                }
                gamma_sum[m] = gamma_sum_temp;
            }

            //compute_parameter();
            //printf("elbo: %f\n", compute_elbo());
        } 
        train_time_lang1 = ((double)(clock()-tStart) / CLOCKS_PER_SEC) / 2.0;
        train_time_lang2 = train_time_lang1;

    } else if (train_method == 2) {
        printf("Use training method 2!\n");

        printf("Do %d iterations for language 1!\n", niters);

        // recording training time
        clock_t tStart = clock();

        for (int iter = 1; iter <= niters; iter++) {
            printf("Iteration %d...\n", iter);

            // compute digamma function of lambda and lambda_sum
            for (k = 0; k < K; k++) {
                for (v = 0; v < V; v++) {
                    dig_lambda[k][v] = exp( boost::math::digamma(lambda[k][v]) );
                }
                dig_lambda_sum[k] = exp( boost::math::digamma(lambda_sum[k]) );
            }
            // compute digamma function of gamma 
            for (m = 0; m < M; m++) {
                for (k = 0; k < K; k++) {
                    dig_gamma[m][k] = exp( boost::math::digamma(gamma[m][k]) );
                }
            }

            // initialize lambda with eta
            for (k = 0; k < K; k++) {
                for (v = 0; v < V; v++) {
                    lambda[k][v] = eta;
                }
            }
            // initialize gamma with alpha
            for (m = 0; m < M; m++) {
                for (k = 0; k < K; k++) {
                    gamma[m][k] = alpha;
                }
            }

            // update phi
            for (m = 0; m < M; m++) {
                int N = ptrndata->docs[m]->length;
                // update phi
                for (n = 0; n < N; n++) {
                    int word = ptrndata->docs[m]->words[n];
                    double phi_sum_temp = 0.0;
                    for (k = 0; k < K; k++) {
                        phi[n][k] =  dig_gamma[m][k] * dig_lambda[k][word] / dig_lambda_sum[k];
                        phi_sum_temp += phi[n][k];
                    }
                    //phi_sum[m][n] = phi_sum_temp;
                    // normalize phi
                    for (k = 0; k < K; k++) {
                        phi[n][k] /= phi_sum_temp;
                        // update lambda
                        lambda[k][word] += phi[n][k];
                        // update gamma
                        gamma[m][k] += phi[n][k];
                    }
                }
            }

            // update lambda_sum
            for (k = 0; k < K; k++) {
                double lambda_sum_temp = 0.0;
                for (v = 0; v < V; v++) {
                    lambda_sum_temp += lambda[k][v];
                }
                lambda_sum[k] = lambda_sum_temp;
            }
            
            // update gamma_sum
            for (m = 0; m < M; m++) {
                double gamma_sum_temp = 0.0;
                for ( k = 0; k < K; k++) {
                    gamma_sum_temp += gamma[m][k];
                }
                gamma_sum[m] = gamma_sum_temp;
            }
        }
        train_time_lang1 = ((double)(clock()-tStart) / CLOCKS_PER_SEC);

        //if (greedy_init) {
        //    printf("Use greedy initialization!\n");
        //    init_est_greedy();
        //}

        printf("Do %d iterations for language 2!\n", niters2);

        // recording training time
        tStart = clock();

        for (int iter = 1; iter <= niters2; iter++) {
            printf("Iteration %d...\n", iter);

            // compute digamma function of lambda2 and lambda_sum2
            for (k2 = 0; k2 < K; k2++) {
                for (v2 = 0; v2 < V2; v2++) {
                    dig_lambda2[k2][v2] = exp( boost::math::digamma(lambda2[k2][v2]) );
                }
                dig_lambda_sum2[k2] = exp( boost::math::digamma(lambda_sum2[k2]) );
            }

            // initialize lambda2 with eta
            for (k2 = 0; k2 < K; k2++) {
                for (v2 = 0; v2 < V2; v2++) {
                    lambda2[k2][v2] = eta;
                }
            }

            // update phi2
            for (m2 = 0; m2 < M2; m2++) {
                int N2 = ptrndata2->docs[m2]->length;
                // update phi2
                for (n2 = 0; n2 < N2; n2++) {
                    int word2 = ptrndata2->docs[m2]->words[n2];
                    double phi_sum_temp2 = 0.0;
                    for (k2 = 0; k2 < K; k2++) {
                        phi2[n2][k2] =  dig_gamma[m2][k2] * dig_lambda2[k2][word2] / dig_lambda_sum2[k2];
                        phi_sum_temp2 += phi2[n2][k2];
                    }
                    //phi_sum2[m2][n2] = phi_sum_temp2;
                    // normalize phi2
                    for (k2 = 0; k2 < K; k2++) {
                        phi2[n2][k2] /= phi_sum_temp2;
                        // update lambda2
                        lambda2[k2][word2] += phi2[n2][k2];
                    }
                }
            }

            // update lambda_sum2
            for (k2 = 0; k2 < K; k2++) {
                double lambda_sum_temp2 = 0.0;
                for (v2 = 0; v2 < V2; v2++) {
                    lambda_sum_temp2 += lambda2[k2][v2];
                }
                lambda_sum2[k2] = lambda_sum_temp2;
            }
        }
        train_time_lang2 = ((double)(clock()-tStart) / CLOCKS_PER_SEC);

    }

    printf("Training time of language 1: %f\n", train_time_lang1);
    printf("Training time of language 2: %f\n", train_time_lang2);

    compute_parameter();
    printf("Saving the parameters...\n");
    save_model();
}

//double model::compute_elbo() {
//    int m, n, k, v;
//    int m2, n2, k2, v2;
//    double elbo = 0.0;
//    // for language 1
//    for (k = 0; k < K; k++) {
//        elbo += lgamma(V * eta) - lgamma(lambda_sum[k]);
//        double temp = log(dig_lambda_sum[k]);
//        for (v = 0; v < V; v++) {
//            elbo += lgamma(lambda[k][v]) - lgamma(eta);
//            elbo += (eta - lambda[k][v]) * (log(dig_lambda[k][v]) - temp);
//        }
//    }
//    // for language 2
//    for (k2 = 0; k2 < K; k2++) {
//        elbo += lgamma(V2 * eta) - lgamma(lambda_sum2[k2]);
//        double temp = log(dig_lambda_sum2[k2]);
//        for (v2 = 0; v2 < V2; v2++) {
//            elbo += lgamma(lambda2[k2][v2]) - lgamma(eta);
//            elbo += (eta - lambda2[k2][v2]) * (log(dig_lambda2[k2][v2]) - temp);
//        }
//    }
//    // for both languages
//    for (m = 0; m < M; m++) {
//        elbo += lgamma(K * alpha) - lgamma(gamma_sum[m]);
//        double dig_gamma_sum = boost::math::digamma(gamma_sum[m]);
//        for (k = 0; k < K; k++) {
//            elbo += lgamma(gamma[m][k]) - lgamma(alpha);
//            elbo += (alpha - gamma[m][k]) * (log(dig_gamma[m][k]) - dig_gamma_sum);
//        }
//    }
//    // for language 1
//    for (m = 0; m < M; m++) {
//        int N = ptrndata->docs[m]->length;
//        for (n = 0; n < N; n++) {
//            for (k = 0; k < K; k++) {
//                elbo += phi[m][n][k] * (log(theta[m][k]) - log(phi[m][n][k]));
//            }
//        }
//    }
//    // for language 2
//    for (m2 = 0; m2 < M2; m2++) {
//        int N2 = ptrndata2->docs[m2]->length;
//        for (n2 = 0; n2 < N2; n2++) {
//            for (k2 = 0; k2 < K; k2++) {
//                elbo += phi2[m2][n2][k2] * (log(theta[m2][k2]) - log(phi2[m2][n2][k2]));
//            }
//        }
//    }
//    // for language 1
//    for (m = 0; m < M; m++) {
//        int N = ptrndata->docs[m]->length;
//        for (n = 0; n < N; n++) {
//            int word = ptrndata->docs[m]->words[n];
//            for (k = 0; k < K; k++) {
//                elbo += phi[m][n][k] * (log(dig_lambda[k][word]) - log(dig_lambda_sum[k]));
//            }
//        }
//    }
//    // for language 2
//    for (m2 = 0; m2 < M2; m2++) {
//        int N2 = ptrndata2->docs[m2]->length;
//        for (n2 = 0; n2 < N2; n2++) {
//            int word2 = ptrndata2->docs[m2]->words[n2];
//            for (k2 = 0; k2 < K; k2++) {
//                elbo += phi2[m2][n2][k2] * (log(dig_lambda2[k2][word2]) - log(dig_lambda_sum2[k2]));
//            }
//        }
//    }
//
//    return elbo;
//}

void model::compute_parameter() {
    int m, n, k, v;
    int k2, v2;

    // compute theta using gamma
    for (m = 0; m < M; m++) {
        for (k = 0; k < K; k++) {
            theta[m][k] = gamma[m][k] / gamma_sum[m];
        }
    }
    
    // compute beta using lambda and lambda_sum
    for (k = 0; k < K; k++) {
        for (v = 0; v < V; v++) {
            beta[k][v] = lambda[k][v] / lambda_sum[k];
        }
    }
    // compute beta2 using lambda2 and lambda_sum2
    for (k2 = 0; k2 < K; k2++) {
        for (v2 = 0; v2 < V2; v2++) {
            beta2[k2][v2] = lambda2[k2][v2] / lambda_sum2[k2];
        }
    }
}

int model::save_model() {
    if (save_model_others(dir + model_name + others_suffix, 1)) {
        return 1;
    }
    if (save_model_others(dir2 + model_name2 + others_suffix, 2)) {
        return 1;
    }
    if (save_model_theta(dir + model_name + theta_suffix)) {
        return 1;
    }
    if (save_model_beta(dir + model_name + beta_suffix, 1)) {
        return 1;
    }
    if (save_model_beta(dir2 + model_name2 + beta_suffix, 2)) {
        return 1;
    }
    if (twords > 0) {
        if (save_model_twords(dir + model_name + twords_suffix, 1)) {
            return 1;
        }
        if (save_model_twords(dir2 + model_name2 + twords_suffix, 2)) {
            return 1;
        }
    }
}

int model::save_model_theta(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
        printf("Cannot open file %s to save!\n", filename.c_str());
        return 1;
    }
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < K; j++) {
            fprintf(fout, "%f ", theta[i][j]);
        }
        fprintf(fout, "\n");
    }
    fclose(fout);
    return 0;
}

int model::save_model_beta(string filename, int language) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
        printf("Cannot open file %s to save!\n", filename.c_str());
        return 1;
    }
    if (language == 1) {
        for (int i = 0; i < K; i++) {
            for (int j = 0; j < V; j++) {
                fprintf(fout, "%f ", beta[i][j]);
            }
            fprintf(fout, "\n");
        }
    } else if (language == 2) {
        for (int i = 0; i < K; i++) {
            for (int j = 0; j < V2; j++) {
                fprintf(fout, "%f ", beta2[i][j]);
            }
            fprintf(fout, "\n");
        }
    }

    fclose(fout);
    return 0;
}

int model::save_model_others(string filename, int language) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
        printf("Cannot open file %s to save!\n", filename.c_str());
        return 1;
    }
    if (language == 1) {
        fprintf(fout, "alpha=%f\n", alpha);
        fprintf(fout, "beta=%f\n", eta);
        fprintf(fout, "ntopics=%d\n", K);
        fprintf(fout, "ndocs=%d\n", M);
        fprintf(fout, "nwords=%d\n", V);
    } else if (language == 2) {
        fprintf(fout, "alpha=%f\n", alpha);
        fprintf(fout, "beta=%f\n", eta);
        fprintf(fout, "ntopics=%d\n", K);
        fprintf(fout, "ndocs=%d\n", M2);
        fprintf(fout, "nwords=%d\n", V2);
    }
    if (train_method == 1) {
        fprintf(fout, "niters=%d\n", niters);
    } else if (train_method == 2) {
        if (language == 1) {
            fprintf(fout, "niters=%d\n", niters);
        } else if (language == 2) {
            fprintf(fout, "niters=%d\n", niters2);
        }
    }

    fclose(fout);
    return 0;
}

int model::save_model_twords(string filename, int language) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
        printf("Cannot open file %s to save!\n", filename.c_str());
        return 1;
    }

    if (language == 1) {
        if (twords > V) {
            twords = V;
        }
        mapid2word::iterator it;
        for (int k = 0; k < K; k++) {
            vector<pair<int, double> > words_probs;
            pair<int, double> word_prob;
            for (int v = 0; v < V; v++) {
                word_prob.first = v;
                word_prob.second = beta[k][v];
                words_probs.push_back(word_prob);
            }
            // quick sort to sort topic_word probability
            utils::quicksort(words_probs, 0, words_probs.size() - 1);
            fprintf(fout, "Topic %dth:\n", k);
            for (int i = 0; i < twords; i++) {
                it = id2word.find(words_probs[i].first);
                if (it != id2word.end()) {
                    fprintf(fout, "\t%s   %f\n", (it->second).c_str(), words_probs[i].second);
                }
            }
        }
    } else if (language == 2) {
        if (twords > V2) {
            twords = V2;
        }
        mapid2word::iterator it;
        for (int k = 0; k < K; k++) {
            vector<pair<int, double> > words_probs;
            pair<int, double> word_prob;
            for (int v = 0; v < V2; v++) {
                word_prob.first = v;
                word_prob.second = beta2[k][v];
                words_probs.push_back(word_prob);
            }
            // quick sort to sort topic_word probability
            utils::quicksort(words_probs, 0, words_probs.size() - 1);
            fprintf(fout, "Topic %dth:\n", k);
            for (int i = 0; i < twords; i++) {
                it = id2word2.find(words_probs[i].first);
                if (it != id2word2.end()) {
                    fprintf(fout, "\t%s   %f\n", (it->second).c_str(), words_probs[i].second);
                }
            }
        }
    }

    fclose(fout);
    return 0;
}

model::~model() {
    if (phi) {
        for (int n = 0; n < ptrndata->N_max; n++) {
            delete[] phi[n];
        }
    }
    if (phi2) {
        for (int n2 = 0; n2 < ptrndata2->N_max; n2++) {
            delete[] phi2[n2];
        }
    }

    if (gamma) {
        for (int m = 0; m < M; m++) {
            delete[] gamma[m];
        }
    }

    if (dig_gamma) {
        for (int m = 0; m < M; m++) {
            delete[] dig_gamma[m];
        }
    }

    if (gamma_sum) {
        delete[] gamma_sum;
    }

    if (lambda) {
        for (int k = 0; k < K; k++) {
            delete[] lambda[k];
        }
    }
    if (lambda2) {
        for (int k2 = 0; k2 < K; k2++) {
            delete[] lambda2[k2];
        }
    }

    if (dig_lambda) {
        for (int k = 0; k < K; k++) {
            delete[] dig_lambda[k];
        }
    }
    if (dig_lambda2) {
        for (int k2 = 0; k2 < K; k2++) {
            delete[] dig_lambda2[k2];
        }
    }

    if (lambda_sum) {
        delete[] lambda_sum;
    }
    if (lambda_sum2) {
        delete[] lambda_sum2;
    }

    if (dig_lambda_sum) {
        delete[] dig_lambda_sum;
    }
    if (dig_lambda_sum2) {
        delete[] dig_lambda_sum2;
    }

    if (theta) {
        for (int m = 0; m < M; m++) {
            delete[] theta[m];
        }
    }

    if (beta) {
        for (int k = 0; k < K; k++) {
            delete[] beta[k];
        }
    }
    if (beta2) {
        for (int k2 = 0; k2 < K; k2++) {
            delete[] beta2[k2];
        }
    }
}
