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
    string rawstr;
    int length;
    vector<int> topic_set;
    int topic_num;
    double * pp;
    
    document() {
	words = NULL;
	rawstr = "";
	length = 0;	
    topic_num = 0;
    pp = NULL;
    }
    
    document(int length) {
	this->length = length;
	rawstr = "";
	this->words = new int[length];	
    topic_num = 0;
    pp = NULL;
    }
    
    document(int length, int * words) {
	this->length = length;
	rawstr = "";
	this->words = new int[length];
	for (int i = 0; i < length; i++) {
	    this->words[i] = words[i];
	}
    topic_num = 0;
    pp = NULL;
    }

    document(int length, int * words, string rawstr) {
	this->length = length;
	this->rawstr = rawstr;
	this->words = new int[length];
	for (int i = 0; i < length; i++) {
	    this->words[i] = words[i];
	}
    topic_num = 0;
    pp = NULL;
    }
    
    document(vector<int> & doc) {
	this->length = doc.size();
	rawstr = "";
	this->words = new int[length];
	for (int i = 0; i < length; i++) {
	    this->words[i] = doc[i];
	}
    topic_num = 0;
    pp = NULL;
    }

    document(vector<int> & doc, string rawstr) {
	this->length = doc.size();
	this->rawstr = rawstr;
	this->words = new int[length];
	for (int i = 0; i < length; i++) {
	    this->words[i] = doc[i];
	}
    topic_num = 0;
    pp = NULL;
    }
    
    ~document() {
	if (words) {
	    delete [] words;
	}
     if (pp) {
         delete pp;
     }
    }
};

class dataset {
public:
    document ** docs;
    document ** _docs; // used only for inference
    map<int, int> _id2id; // also used only for inference
    int M; // number of documents
    int V; // number of words
    
    dataset() {
	docs = NULL;
	_docs = NULL;
	M = 0;
	V = 0;
    }
    
    dataset(int M) {
	this->M = M;
	this->V = 0;
	docs = new document*[M];	
	_docs = NULL;
    }   
    
    ~dataset() {
	if (docs) {
	    for (int i = 0; i < M; i++) {
		delete docs[i];
	    }
	}
	delete docs;
	
	if (_docs) {
	    for (int i = 0; i < M; i++) {
		delete _docs[i];		
	    }
	}
	delete _docs;	
    }
    
    void deallocate() {
	if (docs) {
	    for (int i = 0; i < M; i++) {
		delete docs[i];
	    }
	}
	delete docs;
	docs = NULL;

	if (_docs) {
	    for (int i = 0; i < M; i++) {
		delete _docs[i];
	    }
	}
	delete _docs;
	_docs = NULL;
    }
    
    void add_doc(document * doc, int idx) {
	if (0 <= idx && idx < M) {
	    docs[idx] = doc;
	}
    }   
    
    void _add_doc(document * doc, int idx) {
	if (0 <= idx && idx < M) {
	    _docs[idx] = doc;
	}
    }       

    static int write_wordmap(string wordmapfile, mapword2id * pword2id);
    static int read_wordmap(string wordmapfile, mapword2id * pword2id);
    static int read_wordmap(string wordmapfile, mapid2word * pid2word);
    
    int read_trndata(string dfile, string wordmapfile);
    int read_newdata(string dfile, string wordmapfile);
    int read_newdata_withrawstrs(string dfile, string wordmapfile);
};

#endif

