#include <stdio.h>
#include <stdlib.h>
#include "constants.h"
#include "strtokenizer.h"
#include "dataset.h"

using namespace std;

int dataset::write_wordmap(string wordmapfile, mapword2id * pword2id) {
    FILE * fout = fopen(wordmapfile.c_str(), "w");
    if (!fout) {
	    printf("Cannot open file %s to write!\n", wordmapfile.c_str());
	    return 1;
    }    
    
    mapword2id::iterator it;
    fprintf(fout, "%d\n", pword2id->size());
    for (it = pword2id->begin(); it != pword2id->end(); it++) {
	    fprintf(fout, "%s %d\n", (it->first).c_str(), it->second);
    }
    
    fclose(fout);
    
    return 0;
}

int dataset::read_wordmap(string wordmapfile, mapword2id * pword2id) {
    pword2id->clear();
    
    FILE * fin = fopen(wordmapfile.c_str(), "r");
    if (!fin) {
	    printf("Cannot open file %s to read!\n", wordmapfile.c_str());
	    return 1;
    }    
    
    char buff[BUFF_SIZE_SHORT];
    string line;
    
    fgets(buff, BUFF_SIZE_SHORT - 1, fin);
    int nwords = atoi(buff);
    
    for (int i = 0; i < nwords; i++) {
	    fgets(buff, BUFF_SIZE_SHORT - 1, fin);
	    line = buff;
	    
	    strtokenizer strtok(line, " \t\r\n");
	    if (strtok.count_tokens() != 2) {
	        continue;
	    }
	    
	    pword2id->insert(pair<string, int>(strtok.token(0), atoi(strtok.token(1).c_str())));
    }
    
    fclose(fin);
    
    return 0;
}

int dataset::read_wordmap(string wordmapfile, mapid2word * pid2word) {
    pid2word->clear();
    
    FILE * fin = fopen(wordmapfile.c_str(), "r");
    if (!fin) {
	    printf("Cannot open file %s to read!\n", wordmapfile.c_str());
	    return 1;
    }    
    
    char buff[BUFF_SIZE_SHORT];
    string line;
    
    fgets(buff, BUFF_SIZE_SHORT - 1, fin);
    int nwords = atoi(buff);
    
    for (int i = 0; i < nwords; i++) {
	    fgets(buff, BUFF_SIZE_SHORT - 1, fin);
	    line = buff;
	    
	    strtokenizer strtok(line, " \t\r\n");
	    if (strtok.count_tokens() != 2) {
	        continue;
	    }
	    
	    pid2word->insert(pair<int, string>(atoi(strtok.token(1).c_str()), strtok.token(0)));
    }
    
    fclose(fin);
    
    return 0;
}

int dataset::read_trndata(string dfile, string wordmapfile) {
    mapword2id word2id;
    
    FILE * fin = fopen(dfile.c_str(), "r");
    if (!fin) {
	    printf("Cannot open file %s to read!\n", dfile.c_str());
	    return 1;
    }   
    
    mapword2id::iterator it;    
    char buff[BUFF_SIZE_LONG];
    string line;
    
    // get the number of documents
    fgets(buff, BUFF_SIZE_LONG - 1, fin);
    M = atoi(buff);
    if (M <= 0) {
	    printf("No document available!\n");
	    return 1;
    }
    
    // allocate memory for corpus
    if (docs) {
	    deallocate();
    } else {
	    docs = new document*[M];
    }
    
    // set size of vocabulary to zero
    V = 0;
    
    for (int i = 0; i < M; i++) {
	    fgets(buff, BUFF_SIZE_LONG - 1, fin);
	    line = buff;
	    strtokenizer strtok(line, " \t\r\n");
	    int length = strtok.count_tokens();

	    if (length <= 0) {
            printf("i: %d\n", i);
            printf("length: %d\n", length);
	        printf("Invalid (empty) document!\n");
	        deallocate();
	        M = V = 0;
	        return 1;
	    }
	    
	    // allocate new document
	    document * pdoc = new document(length);
	    
	    for (int j = 0; j < length; j++) {
	        it = word2id.find(strtok.token(j));
	        if (it == word2id.end()) {
	    	    // word not found, i.e., new word
	    	    pdoc->words[j] = word2id.size();
	    	    word2id.insert(pair<string, int>(strtok.token(j), word2id.size()));
	        } else {
	    	    pdoc->words[j] = it->second;
	        }
	    }
	    
	    // add new doc to the corpus
	    add_doc(pdoc, i);

        // update N_max
        if (length > this->N_max) {
            this->N_max = length;
        }
    }
    
    fclose(fin);
    
    // write word map to file
    if (write_wordmap(wordmapfile, &word2id)) {
	    return 1;
    }
    
    // update size of vocabulary
    V = word2id.size();
    
    return 0;
}

