#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "constants.h"
#include "strtokenizer.h"
#include "utils.h"
#include "dataset.h"
#include "model.h"
#include <set>

using namespace std;

void model::set_default_values() {
    wordmapfile = "wordmap.txt";
    trainlogfile = "trainlog.txt";
    theta_suffix = ".theta";
    phi_suffix = ".phi";
    others_suffix = ".others";
    twords_suffix = ".twords";
    
    dir = "./";
    dfile = "trndocs.dat";
    model_name = "model-final";    
    model_status = MODEL_STATUS_UNKNOWN;
    train_method = 1;
    lang1_trained = false;

    train_time_lang1 = 0.0;
    train_time_lang2 = 0.0;
    
    ptrndata = NULL;
    pnewdata = NULL;
    
    M = 0;
    V = 0;
    K = 100;
    alpha = 50.0 / K;
    beta = 0.1;
    niters = 500;
    niters2 = 500;
    liter = 0;
    twords = 0;
    withrawstrs = 0;
    
    p = NULL;
    z = NULL;
    nw = NULL;
    nd = NULL;
    nwsum = NULL;
    ndsum = NULL;
    theta = NULL;
    phi = NULL;
    
    newM = 0;
    newV = 0;
    newz = NULL;
    newnw = NULL;
    newnd = NULL;
    newnwsum = NULL;
    newndsum = NULL;
    newtheta = NULL;
    language = 0;

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    wordmapfile2 = "wordmap2.txt";
    trainlogfile2 = "trainlog2.txt";

    dir2 = "./";
    dfile2 = "trndocs.dat";
    model_name2 = "model-final2";    

    ptrndata2 = NULL;

    M2 = 0;
    V2 = 0;

    p2 = NULL;
    z2 = NULL;
    nw2 = NULL;
    nd2 = NULL;
    nwsum2 = NULL;
    ndsum2 = NULL;
    theta2 = NULL;
    phi2 = NULL;

    topics_num = 0;
    top_topic_lang1 = NULL;
    greedy_init = false;
    greedy_sampling = false;
}

int model::parse_args(int argc, char ** argv) {
    return utils::parse_args(argc, argv, this);
}

int model::init(int argc, char ** argv) {

    srandom(1); // initialize for random number generation

    // call parse_args
    if (parse_args(argc, argv)) {
	    return 1;
    }
    
    if (model_status == MODEL_STATUS_EST) {
	    // estimating the model from scratch
	    if (init_est()) {
	        return 1;
	    }
	
    } else if (model_status == MODEL_STATUS_INF) {
	    // do inference
	    if (init_inf()) {
	        return 1;
	    }
    }
    
    return 0;
}

int model::load_phi(string model_name) {
    int i, j, k, w;

    string filename = dir + model_name + phi_suffix;
    FILE * fin = fopen(filename.c_str(), "r");
    if (!fin) {
	    printf("Cannot open file %s to load model!\n", filename.c_str());
	    return 1;
    }

    char buff[BUFF_SIZE_LONG];
    string line;

    phi = new double*[K];
    for (k = 0; k < K; k++) {
        phi[k] = new double[V];
    }    

    for (k = 0; k < K; k++) {
	    char * pointer = fgets(buff, BUFF_SIZE_LONG, fin);
	    if (!pointer) {
	        printf("Invalid phi file, check the number of topics!\n");
	        return 1;
	    }
	    line = buff;
	    strtokenizer strtok(line, " \t\r\n");
	    int length = strtok.count_tokens();
        if (length != V) {
            printf("The vocabulary size is wrong!\n");
        }
        for (w = 0; w < V; w++) {
            phi[k][w] = atof(strtok.token(w).c_str());
        }
    }

    fclose(fin);
    
    return 0;
}

int model::load_theta(string model_name) {
    int i, j, k, m, n;

    string filename = dir + model_name + theta_suffix;
    FILE * fin = fopen(filename.c_str(), "r");
    if (!fin) {
	    return 1;
    }

    char buff[BUFF_SIZE_LONG];
    string line;

    theta = new double*[M];
    for (m = 0; m < M; m++) {
        theta[m] = new double[K];
    }

    for (m = 0; m < M; m++) {
	    char * pointer = fgets(buff, BUFF_SIZE_LONG, fin);
	    if (!pointer) {
	        printf("Invalid theta file, check the number of documents!\n");
	        return 1;
	    }
	    line = buff;
	    strtokenizer strtok(line, " \t\r\n");
	    int length = strtok.count_tokens();
        if (length != K) {
            printf("The topic number is wrong!\n");
        }
        for (k = 0; k < K; k++) {
            theta[m][k] = atof(strtok.token(k).c_str());
        }
    }

    fclose(fin);

    if (greedy_init) {
        // load top topic file
        filename = dir + "model-final-1.top_topic";
        fin = fopen(filename.c_str(), "r");
        if (!fin) {
            printf("Cannot load top_topic file!\n");
            return 1;
        }
        top_topic_lang1 = new int[M];
        for (m=0; m < M; m++) {
            fgets(buff, BUFF_SIZE_LONG - 1, fin);
            top_topic_lang1[m] = atoi(buff);
        }
        fclose(fin);
    }

    if (greedy_sampling) {
        // load topic set file
        ptrndata = new dataset(M);
        for (m=0; m < M; m++) {
            ptrndata->docs[m] = new document();
        }

        string filename = dir + "model-final-1.topic_set";
        fin = fopen(filename.c_str(), "r");
        if (!fin) {
            printf("Cannot load topic_set file!\n");
            return 1;
        }
        for (m=0; m < M; m++) {
            fgets(buff, BUFF_SIZE_LONG - 1, fin);
            line = buff;
	        strtokenizer strtok(line, " \t\r\n");
            int length = strtok.count_tokens();
            ptrndata->docs[m]->topic_num = length;
            set<int> topic_set;
            for (n=0; n < length; n++) {
                topic_set.insert(atoi(strtok.token(n).c_str()));
            }
            ptrndata->docs[m]->topic_set.assign(topic_set.begin(), topic_set.end());
            ptrndata->docs[m]->pp = new double[length];
        }
        fclose(fin);
    }
    
    return 0;
}

int model::init_est() {
    if (train_method == 1) {
        init_est_lang1();
        init_est_lang2();
    } else if (train_method == 2) {
        if (load_theta("model-final-1")) {
            init_est_lang1();
        } else {
            lang1_trained = true;
        }
        init_est_lang2();
    }

    return 0;
}

int model::init_est_lang1() {
    // for first language
    int m, n, w, k;

    p = new double[K];

    // + read training data
    ptrndata = new dataset;
    if (ptrndata->read_trndata(dir + dfile, dir + wordmapfile)) {
        printf("Fail to read training data of first language!\n");
        return 1;
    }
		
    // + allocate memory and assign values for variables
    M = ptrndata->M;
    V = ptrndata->V;
    // K: from command line or default value
    // alpha, beta: from command line or default values
    // niters, savestep: from command line or default values

    nw = new int*[V];
    for (w = 0; w < V; w++) {
        nw[w] = new int[K];
        for (k = 0; k < K; k++) {
    	    nw[w][k] = 0;
        }
    }

    nd = new int*[M];
    for (m = 0; m < M; m++) {
        nd[m] = new int[K];
        for (k = 0; k < K; k++) {
    	    nd[m][k] = 0;
        }
    }

    nwsum = new int[K];
    for (k = 0; k < K; k++) {
	    nwsum[k] = 0;
    }
	
    ndsum = new int[M];
    for (m = 0; m < M; m++) {
	    ndsum[m] = 0;
    }

    z = new int*[M];
    for (m = 0; m < ptrndata->M; m++) {
	    int N = ptrndata->docs[m]->length;
	    z[m] = new int[N];
	
        // initialize for z
        for (n = 0; n < N; n++) {
    	    int topic = (int)(((double)random() / RAND_MAX) * K);
    	    z[m][n] = topic;
    	    
    	    // number of instances of word i assigned to topic j
    	    nw[ptrndata->docs[m]->words[n]][topic] += 1;
    	    // number of words in document i assigned to topic j
    	    nd[m][topic] += 1;
    	    // total number of words assigned to topic j
    	    nwsum[topic] += 1;
        } 
        // total number of words in document i
        ndsum[m] = N;      
    }
    
    theta = new double*[M];
    for (m = 0; m < M; m++) {
        theta[m] = new double[K];
    }
	
    phi = new double*[K];
    for (k = 0; k < K; k++) {
        phi[k] = new double[V];
    } 
}

int model::init_est_lang2() {
    // for second language
    int m2, n2, w2, k2;

    p2 = new double[K];

    // + read training data of second language
    ptrndata2 = new dataset;
    if (ptrndata2->read_trndata(dir2 + dfile2, dir2 + wordmapfile2)) {
        printf("Fail to read training data of second language!\n");
        return 1;
    }
		
    // + allocate memory and assign values for variables
    M2 = ptrndata2->M;
    V2 = ptrndata2->V;
    // K: from command line or default value
    // alpha, beta: from command line or default values
    // niters, savestep: from command line or default values

    nw2 = new int*[V2];
    for (w2 = 0; w2 < V2; w2++) {
        nw2[w2] = new int[K];
        for (k2 = 0; k2 < K; k2++) {
    	    nw2[w2][k2] = 0;
        }
    }

    nd2 = new int*[M2];
    for (m2 = 0; m2 < M2; m2++) {
        nd2[m2] = new int[K];
        for (k2 = 0; k2 < K; k2++) {
    	    nd2[m2][k2] = 0;
        }
    }

    nwsum2 = new int[K];
    for (k2 = 0; k2 < K; k2++) {
	    nwsum2[k2] = 0;
    }
	
    ndsum2 = new int[M2];
    for (m2 = 0; m2 < M2; m2++) {
	    ndsum2[m2] = 0;
    }

    z2 = new int*[M2];
    for (m2 = 0; m2 < M2; m2++) {
	    int N2 = ptrndata2->docs[m2]->length;
	    z2[m2] = new int[N2];
	
        // initialize for z2
        for (n2 = 0; n2 < N2; n2++) {
    	    int topic2 = (int)(((double)random() / RAND_MAX) * K);
    	    z2[m2][n2] = topic2;
    	    
    	    // number of instances of word i assigned to topic j
    	    nw2[ptrndata2->docs[m2]->words[n2]][topic2] += 1;
    	    // number of words in document i assigned to topic j
    	    nd2[m2][topic2] += 1;
    	    // total number of words assigned to topic j
    	    nwsum2[topic2] += 1;
        } 
        // total number of words in document i
        ndsum2[m2] = N2;      
    }
    
    theta2 = new double*[M2];
    for (m2 = 0; m2 < M2; m2++) {
        theta2[m2] = new double[K];
    }
	
    phi2 = new double*[K];
    for (k2 = 0; k2 < K; k2++) {
        phi2[k2] = new double[V2];
    }
}

int model::init_est_greedy() {
    int m2, n2, w2, k2;

    M2 = ptrndata2->M;
    V2 = ptrndata2->V;

    for (w2 = 0; w2 < V2; w2++) {
        for (k2 = 0; k2 < K; k2++) {
    	    nw2[w2][k2] = 0;
        }
    }

    for (m2 = 0; m2 < M2; m2++) {
        for (k2 = 0; k2 < K; k2++) {
    	    nd2[m2][k2] = 0;
        }
    }

    for (k2 = 0; k2 < K; k2++) {
	    nwsum2[k2] = 0;
    }
	
    for (m2 = 0; m2 < M2; m2++) {
	    ndsum2[m2] = 0;
    }

    for (m2 = 0; m2 < M2; m2++) {
	    int N2 = ptrndata2->docs[m2]->length;
        // initialize z2 with each doc being the topic with largest probability
        int topic2 = top_topic_lang1[m2];
        for (n2 = 0; n2 < N2; n2++) {
    	    z2[m2][n2] = topic2;
    	    
    	    // number of instances of word i assigned to topic j
    	    nw2[ptrndata2->docs[m2]->words[n2]][topic2] += 1;
    	    // number of words in document i assigned to topic j
    	    nd2[m2][topic2] += 1;
    	    // total number of words assigned to topic j
    	    nwsum2[topic2] += 1;
        } 
        // total number of words in document i
        ndsum2[m2] = N2;      
    }
}

void model::estimate() {
    int m, n, w, k, m2, n2, w2, k2;
    if (twords > 0) {
	    // print out top words per topic
        if (!lang1_trained) {
	        dataset::read_wordmap(dir + wordmapfile, &id2word);
        }
	    dataset::read_wordmap(dir2 + wordmapfile2, &id2word2);
    }

    //////////////////////////////////////////////////////////////////////////////////
    if (train_method == 1) {
        printf("Using training method 1!\n");
        printf("Sampling %d iterations for 2 languages!\n", niters);

        // recording training time
        clock_t tStart = clock();

        for (liter = 1; liter <= niters; liter++) {
	        printf("Iteration %d ...\n", liter);
	        
	        // for all z_i of first language
	        for (m = 0; m < M; m++) {
	            for (n = 0; n < ptrndata->docs[m]->length; n++) {
	        	    // (z_i = z[m][n])
	        	    // sample from p(z_i|z_-i, w)
	        	    int topic = sampling(m, n);
	        	    z[m][n] = topic;
	            }
	        }

	        // for all z_i of second language
	        for (m2 = 0; m2 < M2; m2++) {
	            for (n2 = 0; n2 < ptrndata2->docs[m2]->length; n2++) {
	        	    // (z_i = z[m][n])
	        	    // sample from p(z_i|z_-i, w)
	        	    int topic2 = sampling2(m2, n2);
	        	    z2[m2][n2] = topic2;
	            }
	        }
        }
        train_time_lang1 = ((double)(clock()-tStart) / CLOCKS_PER_SEC) / 2.0;
        train_time_lang2 = train_time_lang1;

        compute_theta();

    } else if (train_method == 2) {
        printf("Using training method 2!\n");
        if (!lang1_trained) {
            printf("Sampling %d iterations for language 1!\n", niters);

            // recording training time
            clock_t tStart = clock();

            for (liter = 1; liter <= niters; liter++) {
	            printf("Iteration %d ...\n", liter);
	            
	            // for all z_i of first language
	            for (m = 0; m < M; m++) {
	                for (n = 0; n < ptrndata->docs[m]->length; n++) {
	            	    // (z_i = z[m][n])
	            	    // sample from p(z_i|z_-i, w)
	            	    int topic = sampling(m, n);
	            	    z[m][n] = topic;
	                }
	            }
            }
            train_time_lang1 = (double)(clock()-tStart) / CLOCKS_PER_SEC;

            compute_theta();
            summarize_topics_lang1();

        } else {
            printf("language 1 has been trained and loaded...\n");
            for (m = 0; m < M; m++) {
	            for (k = 0; k < K; k++) {
	                theta2[m][k] = theta[m][k];
	            }
            }
        }

        if (greedy_init) {
            printf("Using greedy initialization!\n");
            init_est_greedy();
        }

        printf("Sampling %d iterations for language 2!\n", niters2);

        // recording training time
        clock_t tStart = clock();

        if (greedy_sampling) {
            printf("Using greedy sampling!\n");
            for (liter = 1; liter <= niters2; liter++) {
	            printf("Iteration %d ...\n", liter);
	            // for all z_i of second language
	            for (m2 = 0; m2 < M2; m2++) {
	                for (n2 = 0; n2 < ptrndata2->docs[m2]->length; n2++) {
	            	    // (z_i = z[m][n])
	            	    // sample from p(z_i|z_-i, w)
	            	    int topic2 = sampling2_greedy(m2, n2);
	            	    z2[m2][n2] = topic2;
	                }
	            }
            }
        } else {
            for (liter = 1; liter <= niters2; liter++) {
	            printf("Iteration %d ...\n", liter);
	            // for all z_i of second language
	            for (m2 = 0; m2 < M2; m2++) {
	                for (n2 = 0; n2 < ptrndata2->docs[m2]->length; n2++) {
	            	    // (z_i = z[m][n])
	            	    // sample from p(z_i|z_-i, w)
	            	    int topic2 = sampling2(m2, n2);
	            	    z2[m2][n2] = topic2;
	                }
	            }
            }
        }

        train_time_lang2 = (double)(clock()-tStart) / CLOCKS_PER_SEC;

    } 

    printf("topics_num:%d\n", topics_num);
    
    printf("Gibbs sampling completed!\n");
    printf("Training time of language 1: %.2f\n", train_time_lang1);
    printf("Training time of language 2: %.2f\n", train_time_lang2);
    printf("Total training time: %.2f\n", train_time_lang1 + train_time_lang2);
    printf("Saving the final model!\n");
    compute_phi();
    liter--;
    save_model(utils::generate_model_name(-1));
}

int model::sampling(int m, int n) {
    // remove z_i from the count variables
    int topic = z[m][n];
    int w = ptrndata->docs[m]->words[n];
    nw[w][topic] -= 1;
    nd[m][topic] -= 1;
    nwsum[topic] -= 1;
    ndsum[m] -= 1;

    double Vbeta = V * beta;
    double Kalpha = K * alpha;    
    // do multinomial sampling via cumulative method
    if (train_method == 1) {
        for (int k = 0; k < K; k++) {
	        p[k] = (nw[w][k] + beta) / (nwsum[k] + Vbeta) *
	        	    (nd[m][k] + nd2[m][k] + alpha) / (ndsum[m] + ndsum2[m] + Kalpha);
        }
    } else {
        for (int k = 0; k < K; k++) {
	        p[k] = (nw[w][k] + beta) / (nwsum[k] + Vbeta) *
	        	    (2*nd[m][k] + alpha) / (2*ndsum[m] + Kalpha);
        }
    } 
    // cumulate multinomial parameters
    for (int k = 1; k < K; k++) {
	    p[k] += p[k - 1];
    }
    // scaled sample because of unnormalized p[]
    double u = ((double)random() / RAND_MAX) * p[K - 1];
    
    for (topic = 0; topic < K; topic++) {
	    if (p[topic] > u) {
	        break;
	    }
    }
    // avoid out of range
    if (topic == K) {
        topic = K - 1;
    }
    
    // add newly estimated z_i to count variables
    nw[w][topic] += 1;
    nd[m][topic] += 1;
    nwsum[topic] += 1;
    ndsum[m] += 1;    
    
    return topic;
}

int model::sampling2(int m2, int n2) {
    // remove z_i from the count variables
    int topic2 = z2[m2][n2];
    int w2 = ptrndata2->docs[m2]->words[n2];
    nw2[w2][topic2] -= 1;
    nd2[m2][topic2] -= 1;
    nwsum2[topic2] -= 1;
    ndsum2[m2] -= 1;

    double Vbeta = V2 * beta;
    double Kalpha = K * alpha;    
    // do multinomial sampling via cumulative method
    if (train_method == 1) {
        for (int k2 = 0; k2 < K; k2++) {
	        p2[k2] = (nw2[w2][k2] + beta) / (nwsum2[k2] + Vbeta) *
	        	    (nd2[m2][k2] + nd[m2][k2] + alpha) / (ndsum2[m2] + ndsum[m2] + Kalpha);
        }
    } else {
        for (int k2 = 0; k2 < K; k2++) {
	        p2[k2] = (nw2[w2][k2] + beta) / (nwsum2[k2] + Vbeta) *
	        	    theta2[m2][k2];
        }
    } 
    // cumulate multinomial parameters
    for (int k2 = 1; k2 < K; k2++) {
	    p2[k2] += p2[k2 - 1];
    }
    // scaled sample because of unnormalized p[]
    double u2 = ((double)random() / RAND_MAX) * p2[K - 1];
    
    for (topic2 = 0; topic2 < K; topic2++) {
	    if (p2[topic2] > u2) {
	        break;
	    }
    }
    // avoid out of range
    if (topic2 == K) {
        topic2 = K - 1;
    }
    
    // add newly estimated z_i to count variables
    nw2[w2][topic2] += 1;
    nd2[m2][topic2] += 1;
    nwsum2[topic2] += 1;
    ndsum2[m2] += 1;    
    
    return topic2;
}

int model::sampling2_greedy(int m2, int n2) {
    // remove z_i from the count variables
    int topic2 = z2[m2][n2];
    int w2 = ptrndata2->docs[m2]->words[n2];
    nw2[w2][topic2] -= 1;
    nd2[m2][topic2] -= 1;
    nwsum2[topic2] -= 1;
    ndsum2[m2] -= 1;

    double Vbeta = V2 * beta;
    double Kalpha = K * alpha;    

    ///////////////////////////////////////////////////////////////////////////////////////////////
    int i;
    vector<int> * topic_set = &ptrndata->docs[m2]->topic_set;
    int size = ptrndata->docs[m2]->topic_num;
    double * pp = ptrndata->docs[m2]->pp;
    for (i=0; i<size; ++i) {
	    pp[i] = (nw2[w2][(*topic_set)[i]] + beta) / (nwsum2[(*topic_set)[i]] + Vbeta) *
	    	    theta2[m2][(*topic_set)[i]];
    }

    // cumulate multinomial parameters
    for (i=1; i<size; ++i) {
	    pp[i] += pp[i-1];
    }
    // scaled sample because of unnormalized p[]
    double u2 = ((double)random() / RAND_MAX) * pp[size-1];

    for (i=0; i < size; ++i) {
	    if (pp[i] > u2) {
	        break;
	    }
    }
    // avoid out of range
    if (i == size) {
        i = size - 1;
    }
    topic2 = (*topic_set)[i];

    // add newly estimated z_i to count variables
    nw2[w2][topic2] += 1;
    nd2[m2][topic2] += 1;
    nwsum2[topic2] += 1;
    ndsum2[topic2] += 1;    
    
    return topic2;
}

void model::compute_theta() {
    if (train_method == 1) {
        for (int m = 0; m < M; m++) {
	        for (int k = 0; k < K; k++) {
	            theta[m][k] = (nd[m][k] + nd2[m][k] + alpha) / (ndsum[m] + ndsum2[m] + K * alpha);
	            theta2[m][k] = theta[m][k];
	        }
        }
    } else if (train_method == 2) {
        for (int m = 0; m < M; m++) {
	        for (int k = 0; k < K; k++) {
	            theta[m][k] = (2 * nd[m][k] + alpha) / (2 * ndsum[m] + K * alpha);
	            theta2[m][k] = theta[m][k];
	        }
        }
    }
}

void model::compute_phi() {
    if (!lang1_trained) {
        for (int k = 0; k < K; k++) {
	        for (int w = 0; w < V; w++) {
	            phi[k][w] = (nw[w][k] + beta) / (nwsum[k] + V * beta);
	        }
        }
    }
    for (int k2 = 0; k2 < K; k2++) {
	    for (int w2 = 0; w2 < V2; w2++) {
	        phi2[k2][w2] = (nw2[w2][k2] + beta) / (nwsum2[k2] + V2 * beta);
	    }
    }
}

void model::summarize_topics_lang1() {
    // find most popular topic for each document of language 1
    top_topic_lang1 = new int[M];
    for (int m=0; m < M; m++) {
        // sort topic probs to choose the largest one
        vector<pair<int, double> > topic_probs;
        pair<int, double> topic_prob;
        for (int k=0; k < K; k++) {
            topic_prob.first = k;
            topic_prob.second = theta[m][k];
            topic_probs.push_back(topic_prob);
        }
        utils::quicksort(topic_probs, 0, topic_probs.size() - 1);
        top_topic_lang1[m] = topic_probs[0].first;
    }

    // save the most popular topics
    string filename = dir + "model-final-1.top_topic";
    FILE * fout = fopen(filename.c_str(), "w");
    for (int m=0; m < M; m++) {
        fprintf(fout, "%d\n", top_topic_lang1[m]);
    }
    fclose(fout);

    // find topic set for each document of language 1
    for (int m=0; m < M; m++) {
        set<int> topic_set;
        for (int n=0; n < ptrndata->docs[m]->length; n++) {
            topic_set.insert(z[m][n]);
        }
        topics_num += topic_set.size();
        ptrndata->docs[m]->topic_num = topic_set.size();
        ptrndata->docs[m]->topic_set.assign(topic_set.begin(), topic_set.end());
        ptrndata->docs[m]->pp = new double[ptrndata->docs[m]->topic_num];
    }

    // save the topic sets
    filename = dir + "model-final-1.topic_set";
    fout = fopen(filename.c_str(), "w");

    for (int m=0; m < M; m++) {
        for (int n=0; n < ptrndata->docs[m]->topic_num; n++) {
            fprintf(fout, "%d ", ptrndata->docs[m]->topic_set[n]);
        }
        fprintf(fout, "\n");
    }
    fclose(fout);
}

int model::init_inf() {
    // estimating the model from a previously estimated one
    int m, n, w, k;

    p = new double[K];

    // load moel, i.e., read z and ptrndata
    if (load_phi(model_name)) {
	    printf("Fail to load phi file of the model!\n");
	    return 1;
    }
    
    // read new data for inference
    pnewdata = new dataset;
    if (withrawstrs) {
        if (language == 1) {
	        if (pnewdata->read_newdata_withrawstrs(dir + dfile, dir + wordmapfile)) {
                printf("Fail to read new data!\n");
            	return 1;
	        }    
        } else {
	        if (pnewdata->read_newdata_withrawstrs(dir + dfile, dir + wordmapfile2)) {
                printf("Fail to read new data!\n");
            	return 1;
	        }    
        }
    } else {
        if (language == 1) {
	        if (pnewdata->read_newdata(dir + dfile, dir + wordmapfile)) {
                printf("Fail to read new data!\n");
            	return 1;
	        }    
        } else {
	        if (pnewdata->read_newdata(dir + dfile, dir + wordmapfile2)) {
                printf("Fail to read new data!\n");
            	return 1;
	        }    
        }
    }
    
    newM = pnewdata->M;
    newV = pnewdata->V;
    
    newnw = new int*[newV];
    for (w = 0; w < newV; w++) {
        newnw[w] = new int[K];
        for (k = 0; k < K; k++) {
    	    newnw[w][k] = 0;
        }
    }
	
    newnd = new int*[newM];
    for (m = 0; m < newM; m++) {
        newnd[m] = new int[K];
        for (k = 0; k < K; k++) {
    	    newnd[m][k] = 0;
        }
    }
	
    newnwsum = new int[K];
    for (k = 0; k < K; k++) {
	    newnwsum[k] = 0;
    }
    
    newndsum = new int[newM];
    for (m = 0; m < newM; m++) {
	    newndsum[m] = 0;
    }

    newz = new int*[newM];
    for (m = 0; m < pnewdata->M; m++) {
	    int N = pnewdata->docs[m]->length;
	    newz[m] = new int[N];

	    // assign values for nw, nd, nwsum, and ndsum	
        for (n = 0; n < N; n++) {
    	    int w = pnewdata->docs[m]->words[n];
    	    int _w = pnewdata->_docs[m]->words[n];
    	    int topic = (int)(((double)random() / RAND_MAX) * K);
    	    newz[m][n] = topic;
    	    
    	    // number of instances of word i assigned to topic j
    	    newnw[_w][topic] += 1;
    	    // number of words in document i assigned to topic j
    	    newnd[m][topic] += 1;
    	    // total number of words assigned to topic j
    	    newnwsum[topic] += 1;
        } 
        // total number of words in document i
        newndsum[m] = N;      
    }    
    
    newtheta = new double*[newM];
    for (m = 0; m < newM; m++) {
        newtheta[m] = new double[K];
    }
	
    return 0;        
}

void model::inference() {
    int m, n, k;
    if (twords > 0) {
	    // print out top words per topic
        if (language == 1) {
	        dataset::read_wordmap(dir + wordmapfile, &id2word);
            printf("read wordmap for language 1!\n");
        }
        else {
	        dataset::read_wordmap(dir + wordmapfile2, &id2word);
            printf("read wordmap for language 2!\n");
        }
    }

    printf("Sampling %d iterations for inference!\n", niters);
    
    for (inf_liter = 1; inf_liter <= niters; inf_liter++) {
	    printf("Iteration %d ...\n", inf_liter);
	    
	    // for all newz_i
	    for (m = 0; m < newM; m++) {
	        for (n = 0; n < pnewdata->docs[m]->length; n++) {
	    	    // (newz_i = newz[m][n])
	    	    // sample from p(z_i|z_-i, w)
	    	    int topic = inf_sampling(m, n);
	    	    newz[m][n] = topic;
	        }
	    }
    }
    
    printf("Gibbs sampling for inference completed!\n");
    printf("Saving the inference outputs!\n");
    compute_newtheta();
    inf_liter--;
    save_inf_model(dfile);
}

int model::inf_sampling(int m, int n) {
    // remove z_i from the count variables
    int topic = newz[m][n];
    int w = pnewdata->docs[m]->words[n];
    int _w = pnewdata->_docs[m]->words[n];
    newnw[_w][topic] -= 1;
    newnd[m][topic] -= 1;
    newnwsum[topic] -= 1;
    newndsum[m] -= 1;
    
    double Vbeta = V * beta;
    double Kalpha = K * alpha;
    // do multinomial sampling via cumulative method
    for (int k = 0; k < K; k++) {
	    p[k] = phi[k][w] *
	    	    (newnd[m][k] + alpha) / (newndsum[m] + Kalpha);
    }
    // cumulate multinomial parameters
    for (int k = 1; k < K; k++) {
	    p[k] += p[k - 1];
    }
    // scaled sample because of unnormalized p[]
    double u = ((double)random() / RAND_MAX) * p[K - 1];
    
    for (topic = 0; topic < K; topic++) {
	    if (p[topic] > u) {
	        break;
	    }
    }
    
    // add newly estimated z_i to count variables
    newnw[_w][topic] += 1;
    newnd[m][topic] += 1;
    newnwsum[topic] += 1;
    newndsum[m] += 1;    
    
    return topic;
}

void model::compute_newtheta() {
    for (int m = 0; m < newM; m++) {
	    for (int k = 0; k < K; k++) {
	        newtheta[m][k] = (newnd[m][k] + alpha) / (newndsum[m] + K * alpha);
	    }
    }
}

int model::save_model(string model_name) {
    if (!lang1_trained) {
        if (save_model_others(dir + model_name + "-1" + others_suffix, 1)) {
	        return 1;
        }
        if (save_model_theta(dir + model_name + "-1" + theta_suffix, 1)) {
	        return 1;
        }
        if (save_model_phi(dir + model_name + "-1" + phi_suffix, 1)) {
	        return 1;
        }
    }

    if (save_model_others(dir2 + model_name + "-2" + others_suffix, 2)) {
	    return 1;
    }
    
    if (save_model_theta(dir2 + model_name + "-2" + theta_suffix, 2)) {
	    return 1;
    }
    
    if (save_model_phi(dir2 + model_name + "-2" + phi_suffix, 2)) {
	    return 1;
    }
    
    if (twords > 0) {
        if (!lang1_trained) {
	        if (save_model_twords(dir + model_name + "-1" + twords_suffix, 1)) {
	            return 1;
	        }
        }
	    if (save_model_twords(dir2 + model_name + "-2" + twords_suffix, 2)) {
	        return 1;
	    }
    }
    
    return 0;
}

int model::save_model_theta(string filename, int language) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	    printf("Cannot open file %s to save!\n", filename.c_str());
	    return 1;
    }
    
    if (language == 1) {
        for (int i = 0; i < M; i++) {
	        for (int j = 0; j < K; j++) {
	            fprintf(fout, "%f ", theta[i][j]);
	        }
	        fprintf(fout, "\n");
        }
    } else {
        for (int i = 0; i < M2; i++) {
	        for (int j = 0; j < K; j++) {
	            fprintf(fout, "%f ", theta2[i][j]);
	        }
	        fprintf(fout, "\n");
        }
    }
    
    fclose(fout);
    
    return 0;
}

int model::save_model_phi(string filename, int language) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	    printf("Cannot open file %s to save!\n", filename.c_str());
	    return 1;
    }
    
    if (language == 1) {
        for (int i = 0; i < K; i++) {
	        for (int j = 0; j < V; j++) {
	            fprintf(fout, "%f ", phi[i][j]);
	        }
	        fprintf(fout, "\n");
        }
    } else {
        for (int i = 0; i < K; i++) {
	        for (int j = 0; j < V2; j++) {
	            fprintf(fout, "%f ", phi2[i][j]);
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

    fprintf(fout, "trnmethod=%d\n", train_method);
    fprintf(fout, "alpha=%f\n", alpha);
    fprintf(fout, "beta=%f\n", beta);
    fprintf(fout, "ntopics=%d\n", K);
    if (language == 1) {
        fprintf(fout, "ndocs=%d\n", M);
        fprintf(fout, "nwords=%d\n", V);
        fprintf(fout, "trntime=%.2f\n", train_time_lang1);
    } else if (language == 2) {
        fprintf(fout, "ndocs=%d\n", M2);
        fprintf(fout, "nwords=%d\n", V2);
        fprintf(fout, "trntime=%.2f\n", train_time_lang2);
    }
    fprintf(fout, "niters=%d\n", niters);
    if (train_method == 2) {
        fprintf(fout, "niters2=%d\n", niters2);
        fprintf(fout, "greedyinit=%d\n", greedy_init);
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
	        for (int w = 0; w < V; w++) {
	            word_prob.first = w;
	            word_prob.second = phi[k][w];
	            words_probs.push_back(word_prob);
	        }
            
                // quick sort to sort word-topic probability
	        utils::quicksort(words_probs, 0, words_probs.size() - 1);
	        
	        fprintf(fout, "Topic %dth:\n", k);
	        for (int i = 0; i < twords; i++) {
	            it = id2word.find(words_probs[i].first);
	            if (it != id2word.end()) {
	        	    fprintf(fout, "\t%s   %f\n", (it->second).c_str(), words_probs[i].second);
	            }
	        }
        }
    } else {
        if (twords > V2) {
	        twords = V2;
        }
        mapid2word::iterator it;
        
        for (int k = 0; k < K; k++) {
	        vector<pair<int, double> > words_probs;
	        pair<int, double> word_prob;
	        for (int w = 0; w < V2; w++) {
	            word_prob.first = w;
	            word_prob.second = phi2[k][w];
	            words_probs.push_back(word_prob);
	        }
            
                // quick sort to sort word-topic probability
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

int model::save_inf_model(string model_name) {
    if (save_inf_model_others(dir + model_name + others_suffix)) {
	    return 1;
    }
    
    if (save_inf_model_newtheta(dir + model_name + theta_suffix)) {
	    return 1;
    }
    
    return 0;
}

int model::save_inf_model_newtheta(string filename) {
    int i, j;

    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	    printf("Cannot open file %s to save!\n", filename.c_str());
	    return 1;
    }
    
    for (i = 0; i < newM; i++) {
	    for (j = 0; j < K; j++) {
	        fprintf(fout, "%f ", newtheta[i][j]);
	    }
	    fprintf(fout, "\n");
    }
    
    fclose(fout);
    
    return 0;
}

int model::save_inf_model_others(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	    printf("Cannot open file %s to save!\n", filename.c_str());
	    return 1;
    }

    fprintf(fout, "alpha=%f\n", alpha);
    fprintf(fout, "beta=%f\n", beta);
    fprintf(fout, "ntopics=%d\n", K);
    fprintf(fout, "ndocs=%d\n", newM);
    fprintf(fout, "nwords=%d\n", newV);
    fprintf(fout, "niters=%d\n", niters);
    
    fclose(fout);    
    
    return 0;
}

model::~model() {
    if (p) {
	    delete p;
    }
    if (p2) {
	    delete p2;
    }

    if (ptrndata) {
	    delete ptrndata;
    }
    if (ptrndata2) {
	    delete ptrndata2;
    }
    
    if (pnewdata) {
	    delete pnewdata;
    }

    if (z) {
	    for (int m = 0; m < M; m++) {
	        if (z[m]) {
	    	    delete z[m];
	        }
	    }
    }
    if (z2) {
	    for (int m = 0; m < M2; m++) {
	        if (z2[m]) {
	    	    delete z2[m];
	        }
	    }
    }
    
    if (nw) {
	    for (int w = 0; w < V; w++) {
	        if (nw[w]) {
	    	    delete nw[w];
	        }
	    }
    }
    if (nw2) {
	    for (int w = 0; w < V2; w++) {
	        if (nw2[w]) {
	    	    delete nw2[w];
	        }
	    }
    }

    if (nd) {
	    for (int m = 0; m < M; m++) {
	        if (nd[m]) {
	    	    delete nd[m];
	        }
	    }
    } 
    if (nd2) {
	    for (int m = 0; m < M2; m++) {
	        if (nd2[m]) {
	    	    delete nd2[m];
	        }
	    }
    } 
    
    if (nwsum) {
	    delete nwsum;
    }   
    if (nwsum2) {
	    delete nwsum2;
    }   
    
    if (ndsum) {
	    delete ndsum;
    }
    if (ndsum2) {
	    delete ndsum2;
    }
    
    if (theta) {
	    for (int m = 0; m < M; m++) {
	        if (theta[m]) {
	    	    delete theta[m];
	        }
	    }
    }
    if (theta2) {
	    for (int m = 0; m < M2; m++) {
	        if (theta2[m]) {
	    	    delete theta2[m];
	        }
	    }
    }
    
    if (phi) {
	    for (int k = 0; k < K; k++) {
	        if (phi[k]) {
	    	    delete phi[k];
	        }
	    }
    }
    if (phi2) {
	    for (int k = 0; k < K; k++) {
	        if (phi2[k]) {
	    	    delete phi2[k];
	        }
	    }
    }

    // only for inference
    if (newz) {
	    for (int m = 0; m < newM; m++) {
	        if (newz[m]) {
	    	    delete newz[m];
	        }
	    }
    }
    
    if (newnw) {
	    for (int w = 0; w < newV; w++) {
	        if (newnw[w]) {
	    	    delete newnw[w];
	        }
	    }
    }

    if (newnd) {
	    for (int m = 0; m < newM; m++) {
	        if (newnd[m]) {
	    	    delete newnd[m];
	        }
	    }
    } 
    
    if (newnwsum) {
	    delete newnwsum;
    }   
    
    if (newtheta) {
	    for (int m = 0; m < newM; m++) {
	        if (newtheta[m]) {
	    	    delete newtheta[m];
	        }
	    }
    }

    if (top_topic_lang1) {
        delete top_topic_lang1;
    }
}

