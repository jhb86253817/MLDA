lda: lda.cpp dataset.o strtokenizer.o utils.o model.o
	g++ -o lda lda.cpp dataset.o strtokenizer.o utils.o model.o
dataset.o: dataset.h dataset.cpp
	g++ -c -o dataset.o dataset.cpp
strtokenizer.o: strtokenizer.h strtokenizer.cpp
	g++ -c -o strtokenizer.o strtokenizer.cpp
utils.o: utils.h utils.cpp
	g++ -c -o utils.o utils.cpp
model.o: model.h model.cpp
	g++ -c -o model.o model.cpp
