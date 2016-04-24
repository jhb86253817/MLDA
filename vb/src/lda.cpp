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
    
	// parameter estimation
	lda.estimate();
    
    return 0;
}

void show_help() {
    printf("Command line usage:\n");
    printf("\tlda -ntopics <int> -niters <int> -niters2 <int> -twords <int> -trnmethod <int> -greedyinit <int> -nestediters <int> -dfile <string> -dfile2 <string>\n");
}

