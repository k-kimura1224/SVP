SRCS=$(wildcard *.cpp)
DEST=$(patsubst %.cpp,%,${SRCS})

CXXFLAGS	=	-I/Users/kimurakeiji/ntl-11.3.0/include/NTL
LIB=-L./ntl-11.3.0/src

all: ${DEST}

clean:
	rm ${DEST} *~

%: %.cpp *.h
#	g++ -g3 -O0 -Wall -static -o $@ $< -lntl -lgmp -lm
	g++ -g3 -O0 -Wall -static -o $(CXXFLAGS) $@ $< ${LIB} -lntl -lgmp -lm

