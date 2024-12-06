all:
	g++ ./main.cpp ./minisat/core/Solver.cc ./minisat/utils/*.cc -I./minisat/ -fpermissive