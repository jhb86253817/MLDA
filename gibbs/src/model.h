#ifndef	_MODEL_H
#define	_MODEL_H

#include "constants.h"
#include "dataset.h"

using namespace std;

// LDA model
class model {
public:
    // fixed options
    string wordmapfile;		// file that contains word map [string -> integer id]
    string wordmapfile2;		// file that contains word map [string -> integer id]
    string trainlogfile;	// training log file
    string trainlogfile2;	// training log file
    string tassign_suffix;	// suffix for topic assignment file
    string theta_suffix;	// suffix for theta file
    string phi_suffix;		// suffix for phi file
    string others_suffix;	// suffix for file containing other parameters
    string twords_suffix;	// suffix for file containing words-per-topics

    string dir;			// model directory
    string dir2;			// model directory
    string dfile;		// data file    
    string dfile2;		// data file    
    string model_name;		// model name
    string model_name2;		// model name
    int model_status;		// model status:
				// MODEL_STATUS_UNKNOWN: unknown status
				// MODEL_STATUS_EST: estimating from scratch
				// MODEL_STATUS_INF: do inference
    int train_method; // 1: joint training 2: separate training 3: separate training and initialization
                      // 4: separate greedy training 5: separate training with loading previous trained model of language 1
    bool lang1_trained;

    dataset * ptrndata;	// pointer to training dataset object
    dataset * ptrndata2; // pointer to training dataset object of second language
    dataset * pnewdata; // pointer to new dataset object

    mapid2word id2word; // word map [int => string]
    mapid2word id2word2; // word map [int => string]
    
    double train_time_lang1; // training time of language 1
    double train_time_lang2; // training time of language 2

    // --- model parameters and variables ---    
    int M; // dataset size (i.e., number of docs)
    int M2; // dataset size (i.e., number of docs)
    int V; // vocabulary size
    int V2; // vocabulary size
    int K; // number of topics
    double alpha, beta; // LDA hyperparameters 
    int niters; // number of Gibbs sampling iterations
    int niters2; // number of Gibbs sampling iterations for train method 2 on language 2
    int liter; // last iteration where the model saved
    int twords; // print out top words per each topic
    int withrawstrs;


    double * p; // temp variable for sampling
    double * p2; // temp variable for sampling
    int ** z; // topic assignments for words, size M x doc.size()
    int ** z2; // topic assignments for words, size M x doc.size()
    int ** nw; // cwt[i][j]: number of instances of word/term i assigned to topic j, size V x K
    int ** nw2; // cwt[i][j]: number of instances of word/term i assigned to topic j, size V x K
    int ** nd; // na[i][j]: number of words in document i assigned to topic j, size M x K
    int ** nd2; // na[i][j]: number of words in document i assigned to topic j, size M x K
    int * nwsum; // nwsum[j]: total number of words assigned to topic j, size K
    int * nwsum2; // nwsum[j]: total number of words assigned to topic j, size K
    int * ndsum; // nasum[i]: total number of words in document i, size M
    int * ndsum2; // nasum[i]: total number of words in document i, size M
    double ** theta; // theta: document-topic distributions, size M x K
    double ** theta2; // theta: document-topic distributions, size M x K
    double ** phi; // phi: topic-word distributions, size K x V
    double ** phi2; // phi: topic-word distributions, size K x V
    
    // for inference only
    int inf_liter;
    int newM;
    int newV;
    int ** newz;
    int ** newnw;
    int ** newnd;
    int * newnwsum;
    int * newndsum;
    double ** newtheta;
    int language;

    int topics_num;
    int * top_topic_lang1; // for storing most popular topic of each document
    bool greedy_init; // for train method 2, whether use greedy initializaion for language 2 
    bool greedy_sampling; // for train method 2, whether use greedy sampling for language 2
    // --------------------------------------
    
    model() {
	    set_default_values();
    }
          
    ~model();
    
    // set default values for variables
    void set_default_values();   

    // parse command line to get options
    int parse_args(int argc, char ** argv);
    
    // initialize the model
    int init(int argc, char ** argv);
    
    // load LDA model to continue estimating or to do inference
    int load_model(string model_name);
    int load_phi(string model_name);
    int load_theta(string model_name);
    
    // save LDA model to files
    // model_name.theta: document-topic distributions
    // model_name.phi: topic-word distributions
    // model_name.others: containing other parameters of the model (alpha, beta, M, V, K)
    int save_model(string model_name);
    int save_model_theta(string filename, int language);
    int save_model_phi(string filename, int language);
    int save_model_others(string filename, int language);
    int save_model_twords(string filename, int language);
    
    // saving inference outputs
    int save_inf_model(string model_name);
    int save_inf_model_newtheta(string filename);
    int save_inf_model_others(string filename);
    
    // init for estimation
    int init_est();
    int init_est_lang1();
    int init_est_lang2();
    int init_est_greedy();
	
    // estimate LDA model using Gibbs sampling
    void estimate();
    int sampling(int m, int n);
    int sampling2(int m, int n);
    int sampling2_greedy(int m, int n);
    void compute_theta();
    void compute_phi();
    void summarize_topics_lang1();

    
    // init for inference
    int init_inf();
    // inference for new (unseen) data based on the estimated LDA model
    void inference();
    int inf_sampling(int m, int n);
    void compute_newtheta();
};

#endif

