build: main.cpp
	g++ -o align main.cpp

sample:
	./align sample_seq_x.txt sample_seq_y.txt

run:
	./align seq_x.txt seq_y.txt
