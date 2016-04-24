#include <stdio.h>
#include <string>
#include <map>
#include "strtokenizer.h"
#include "utils.h"
#include "model.h"
#include <cstdlib>

using namespace std;

int utils::parse_args(int argc, char ** argv, model * pmodel) {
    string dir = "";
    string dir2 = "";
    string model_name = "";
    string dfile = "";
    string dfile2 = "";
    double alpha = -1.0;
    double eta = -1.0;
    int K = 0;
    int niters = 0;
    int niters2 = 0;
    int twords = 0;
    int train_method = 1;
    int greedy_init = 0;
    int nested_iters = 0;

    int i = 0; 
    while (i < argc) {
	    string arg = argv[i];
	    
	    if (arg == "-dir") {
	        dir = argv[++i];	    
	    } else if (arg == "-dir2") {
	        dir2 = argv[++i];	    
	    } else if (arg == "-dfile") {
	        dfile = argv[++i];	    
	    } else if (arg == "-dfile2") {
	        dfile2 = argv[++i];	    
	    } else if (arg == "-model") {
	        model_name = argv[++i];	    	    
	    } else if (arg == "-alpha") {
	        alpha = atof(argv[++i]);	    
	    } else if (arg == "-eta") {
	        eta = atof(argv[++i]);	    
	    } else if (arg == "-ntopics") {
	        K = atoi(argv[++i]);	    
	    } else if (arg == "-niters") {
	        niters = atoi(argv[++i]);	    
	    } else if (arg == "-niters2") {
	        niters2 = atoi(argv[++i]);	    
	    } else if (arg == "-twords") {
	        twords = atoi(argv[++i]);
	    } else if (arg == "-trnmethod") {
	        train_method = atoi(argv[++i]);
	    } else if (arg == "-greedyinit") {
	        greedy_init = atoi(argv[++i]);
	    } else if (arg == "-nestediters") {
	        nested_iters = atoi(argv[++i]);
	    } else {
	        // any more?
	    }	
	    i++;
    }
    
	if (dfile == "") {
	    printf("Please specify the input data file for model estimation!\n");
	    return 1;
	}
	
	if (K > 0) {
	    pmodel->K = K;
	}
	if (alpha >= 0.0) {
	    pmodel->alpha = alpha;
	} 
	if (eta >= 0.0) {
	    pmodel->eta = eta;
	}
	if (niters > 0) {
	    pmodel->niters = niters;
	} 
	if (niters2 > 0) {
	    pmodel->niters2 = niters2;
	} 
	if (twords > 0) {
	    pmodel->twords = twords;
	}
    if (train_method > 0 && train_method < 3) {
        pmodel->train_method = train_method;
    }
    if (greedy_init >= 0 && greedy_init <= 1) {
        pmodel->greedy_init = (greedy_init != 0);
    }
    if (nested_iters > 0) {
        pmodel->nested_iters = nested_iters;
    }
	
	pmodel->dfile = dfile;
	pmodel->dfile2 = dfile2;
	
    // for language 1
	string::size_type idx = dfile.find_last_of("/");			
	if (idx == string::npos) {
	    pmodel->dir = "./";
	} else {
	    pmodel->dir = dfile.substr(0, idx + 1);
	    pmodel->dfile = dfile.substr(idx + 1, dfile.size() - pmodel->dir.size());
	    printf("dir = %s\n", pmodel->dir.c_str());
	    printf("dfile = %s\n", pmodel->dfile.c_str());
	}

    // for language 2
	idx = dfile2.find_last_of("/");			
	if (idx == string::npos) {
	    pmodel->dir2 = "./";
	} else {
	    pmodel->dir2 = dfile2.substr(0, idx + 1);
	    pmodel->dfile2 = dfile2.substr(idx + 1, dfile2.size() - pmodel->dir2.size());
	    printf("dir2 = %s\n", pmodel->dir2.c_str());
	    printf("dfile2 = %s\n", pmodel->dfile2.c_str());
	}
    
    return 0;
}

//int utils::read_and_parse(string filename, model * pmodel) {
//    // open file <model>.others to read:
//    // alpha=?
//    // beta=?
//    // ntopics=?
//    // ndocs=?
//    // nwords=?
//    // niter=? 
//    // trntime=?
//    
//    FILE * fin = fopen(filename.c_str(), "r");
//    if (!fin) {
//	    printf("Cannot open file: %s\n", filename.c_str());
//	    return 1;
//    }
//    
//    char buff[BUFF_SIZE_SHORT];
//    string line;
//    
//    while (fgets(buff, BUFF_SIZE_SHORT - 1, fin)) {
//	    line = buff;
//	    strtokenizer strtok(line, "= \t\r\n");
//	    int count = strtok.count_tokens();
//	    
//	    if (count != 2) {
//	        // invalid, ignore this line
//	        continue;
//	    }
//
//	    string optstr = strtok.token(0);
//	    string optval = strtok.token(1);
//	    
//	    if (optstr == "alpha") {
//	        pmodel->alpha = atof(optval.c_str());
//	    } else if (optstr == "beta") {	    
//	        pmodel->beta = atof(optval.c_str());
//	    } else if (optstr == "ntopics") {
//	        pmodel->K = atoi(optval.c_str());
//	    } else if (optstr == "ndocs") {	   
//	        pmodel->M = atoi(optval.c_str());
//	    } else if (optstr == "nwords") {
//	        pmodel->V = atoi(optval.c_str());
//	    } else if (optstr == "trntime") {
//            pmodel->train_time_lang1 = atof(optval.c_str());
//	    } else {
//	        // any more?
//	    }
//    }
//    
//    fclose(fin);
//    
//    return 0;
//}

//string utils::generate_model_name(int iter) {
//    string model_name = "model-";
//
//    char buff[BUFF_SIZE_SHORT];
//    
//    if (0 <= iter && iter < 10) {
//	    sprintf(buff, "0000%d", iter);
//    } else if (10 <= iter && iter < 100) {
//	    sprintf(buff, "000%d", iter);
//    } else if (100 <= iter && iter < 1000) {
//	    sprintf(buff, "00%d", iter);
//    } else if (1000 <= iter && iter < 10000) {
//	    sprintf(buff, "0%d", iter);
//    } else {
//	    sprintf(buff, "%d", iter);
//    }
//    
//    if (iter >= 0) {
//	    model_name += buff;
//    } else {
//	    model_name += "final";
//    }
//    
//    return model_name;
//}

void utils::sort(vector<double> & probs, vector<int> & words) {
    for (int i = 0; i < probs.size() - 1; i++) {
	    for (int j = i + 1; j < probs.size(); j++) {
	        if (probs[i] < probs[j]) {
	    	    double tempprob = probs[i];
	    	    int tempword = words[i];
	    	    probs[i] = probs[j];
	    	    words[i] = words[j];
	    	    probs[j] = tempprob;
	    	    words[j] = tempword;
	        }
	    }
    }
}

void utils::quicksort(vector<pair<int, double> > & vect, int left, int right) {
    int l_hold, r_hold;
    pair<int, double> pivot;
    
    l_hold = left;
    r_hold = right;    
    int pivotidx = left;
    pivot = vect[pivotidx];

    while (left < right) {
	    while (vect[right].second <= pivot.second && left < right) {
	        right--;
	    }
	    if (left != right) {
	        vect[left] = vect[right];
	        left++;
	    }
	    while (vect[left].second >= pivot.second && left < right) {
	        left++;
	    }
	    if (left != right) {
	        vect[right] = vect[left];
	        right--;
	    }
    }

    vect[left] = pivot;
    pivotidx = left;
    left = l_hold;
    right = r_hold;
    
    if (left < pivotidx) {
	    quicksort(vect, left, pivotidx - 1);
    }
    if (right > pivotidx) {
	    quicksort(vect, pivotidx + 1, right);
    }    
}

