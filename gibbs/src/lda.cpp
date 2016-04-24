#include "model.h"
#include <cstdio>

using namespace std;

void show_help();

int main(int argc, char ** argv) {
    model lda;

    if (lda.init(argc, argv)) {
	    show_help();
	    return 1;
    }
    
    if (lda.model_status == MODEL_STATUS_EST) {
	    // parameter estimation
	    lda.estimate();
    }
    
    if (lda.model_status == MODEL_STATUS_INF) {
	    // do inference
	    lda.inference();
    }

    return 0;
}

void show_help() {
    printf("Command line usage:\n");
    printf("\tlda -est -ntopics <int> -niters <int> -niters2 <int> -twords <int> -trnmethod <int> -greedyinit <int> -greedysampling <int> -dfile <string> -dfile2 <string>\n");
    printf("\tlda -inf -dir <string> -model <string> -niters <int> -twords <int> -lang <int> -dfile <string>\n");
}

