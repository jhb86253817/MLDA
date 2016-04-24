#ifndef	_MODEL_H
#define	_MODEL_H

#include "dataset.h"

using namespace std;

// LDA model
class model {
public:
    // fixed options
    string wordmapfile; // file that contains wordmap [string -> integer id]
    string wordmapfile2; // file that contains wordmap [string -> integer id]
    string theta_suffix; // suffix for theta file
    string beta_suffix; // suffix for beta file
    string others_suffix; // suffix for file containing other parameters
    string twords_suffix; // suffix for file containing words-per-topics

    string dir; // model directory
    string dir2; // model directory
    string dfile; // data file
    string dfile2; // data file
    string model_name; // model name
    string model_name2; // model name
    dataset * ptrndata; // pointer to training dataset object
    dataset * ptrndata2; // pointer to training dataset object
    mapid2word id2word; //word map [int => string]
    mapid2word id2word2; //word map [int => string]

    int train_method; // 1: joint training 2: separate training
    int M; // number of documents
    int M2; // number of documents
    int V; // size of vocabulary
    int V2; // size of vocabulary
    int K; // number of topics
    double alpha, eta; // hyperparameters
    int niters; // number of iterations
    int niters2; // number of iterations of language 2, for training method 2
    int nested_iters; // number of iterations for nested computation in inference
    int twords; // print out top words for each topic

    bool greedy_init; // for training method 2, whether use greedy initialization for language 2

    double train_time_lang1; // training time of language 1
    double train_time_lang2; // training time of language 2

    // variational parameters
    double ** phi; // for latent topic z, size N_max x K
    double ** phi2; // for latent topic z, size N_max x K
    //double ** phi_sum; // sum of each latent topic's unnormalized topic probability, size M x N
    //double ** phi_sum2; // sum of each latent topic's unnormalized topic probability, size M x N
    double ** gamma; // for document parameter theta, size M x K
    double * gamma_sum; // sum of gamma on each topic, length M
    double ** lambda; // for topic parameter beta, size K x V
    double ** lambda2; // for topic parameter beta, size K x V
    double * lambda_sum; // sum of lambda on each topic, length K
    double * lambda_sum2; // sum of lambda on each topic, length K
    double ** dig_gamma; // digamma function of gamma
    double ** dig_lambda; // digamma function of lambda
    double ** dig_lambda2; // digamma function of lambda
    double * dig_lambda_sum; // digamma function of lambda_sum
    double * dig_lambda_sum2; // digamma function of lambda_sum

    // parameters
    double ** theta; // document-topic distribution, size M x K
    double ** beta; // topic-word distribution, size K x V
    double ** beta2; // topic-word distribution, size K x V

    ////////////////////////////////////////////////////////////

    model() {
        set_default_values();
    }
          
    ~model();

    void set_default_values();

    int parse_args(int argc, char ** argv);

    // initialize the model
    int init(int argc, char ** argv);

    int init_est();
    //int init_est_greedy();

    void estimate();

    void compute_parameter();

    //double compute_elbo();

    // save LDA model to files
    int save_model();
    int save_model_theta(string filename);
    int save_model_beta(string filename, int language);
    int save_model_others(string filename, int language);
    int save_model_twords(string filename, int language);
};

#endif
