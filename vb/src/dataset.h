#ifndef	_DATASET_H
#define	_DATASET_H

#include <string>
#include <vector>
#include <map>

using namespace std;

// map of words/terms [string => int]
typedef map<string, int> mapword2id;
// map of words/terms [int => string]
typedef map<int, string> mapid2word;

class document {
public:
    int * words;
    int length;
    
    document() {
	words = NULL;
	length = 0;	
    }
    
    document(int length) {
	this->length = length;
	this->words = new int[length];	
    }
    
    document(int length, int * words) {
	this->length = length;
	this->words = new int[length];
	for (int i = 0; i < length; i++) {
	    this->words[i] = words[i];
	}
    }

    ~document() {
	if (words) {
	    delete [] words;
	}
    }
};

class dataset {
public:
    document ** docs;
    int M; // number of documents
    int V; // size of vocabulary
    int N_max; // max number of words of a doc
    
    dataset() {
	docs = NULL;
	M = 0;
	V = 0;
    N_max = 0;
    }
    
    dataset(int M) {
	this->M = M;
	this->V = 0;
    this->N_max = 0;
	docs = new document*[M];	
    }   
    
    ~dataset() {
	if (docs) {
	    for (int i = 0; i < M; i++) {
		delete docs[i];
	    }
	}
	delete docs;
    }
    
    void deallocate() {
	if (docs) {
	    for (int i = 0; i < M; i++) {
		delete docs[i];
	    }
	}
	delete docs;
	docs = NULL;
    }
    
    void add_doc(document * doc, int idx) {
	if (0 <= idx && idx < M) {
	    docs[idx] = doc;
	}
    }   

    static int write_wordmap(string wordmapfile, mapword2id * pword2id);
    static int read_wordmap(string wordmapfile, mapword2id * pword2id);
    static int read_wordmap(string wordmapfile, mapid2word * pid2word);
    
    int read_trndata(string dfile, string wordmapfile);
};

#endif

