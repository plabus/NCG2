CXX=clang++
CXXFLAGS= -std=c++14 -O2
WARNINGS= -Wall -Wpedantic
LIBRARY=

CXXFLAGS+=${WARNINGS}

.PHONY: all
all: ncg++


ncg++: markov_chain.o
	${CXX} -o $@ $^ ${CXXFLAGS} ${LIBRARY}



%.o: %.c
	${CXX} -c ${CXXFLAGS} $< -o $@


clean:
	rm -f *.o
