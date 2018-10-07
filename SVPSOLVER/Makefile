# Makefile

OPT=opt
PG=false
SRCDIR	=	./src
OBJDIR	=	./obj/$(OPT)
BINDIR	=	./bin

MAINOBJ		=	main.o read_data.o qpsolver.o vector.o \
					probdata.o svpsolver.o solution.o svps_solve.o \
					heur_init.o solution_pool.o node.o stopwatch.o \
					relax.o branch.o svps_parallel.o \
					heur.o heur_unitsphere.o Schmidt_manager.o \
					heur_quadratic.o cut.o cut_pool.o cut_oa.o qpdata.o \
					testwatch.o nodelist.o
MAINSRC		=	$(addprefix $(SRCDIR)/,$(MAINOBJ:.o=.cpp))
MAINOBJFILES	=	$(addprefix $(OBJDIR)/,$(MAINOBJ))

DEPFILE		=	$(SRCDIR)/svp.depend

NAME		=	SVPSOLVER
ifeq ($(OPT),dbg)
	MAIN		=	$(NAME).dbg
else
	MAIN		=	$(NAME).opt
endif

MAINFILE	=	$(BINDIR)/$(MAIN)

CXX	=	/usr/local/bin/g++-7

CXXFLAGS	=	-Wall -llapack -lpthread -fopenmp -O2 -std=c++14

#Mac
#CXXFLAGS	+=	-lblas

#Linux
CXXFLAGS 	+=	-lopenblas

ifeq ($(OPT),opt)
	CXXFLAGS	+=	-DNDEBUG
endif

ifeq ($(PG),true)
	CXXFLAGS	+=	-pg
endif

#Mac
#CXXFLAGS	+=	-I/usr/local/Cellar/openblas/0.2.18_2/include

#Linux
CXXLDFLAGS	=	-L/usr/lib64/libopenblas.so

.PHONY: all
all: $(MAIN)

$(MAIN):	$(MAINOBJFILES)
		@$(CXX) -o $(MAINFILE) $(CXXLDFLAGS) $(MAINOBJFILES) $(CXXFLAGS)


$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp
		@echo "-> compiling $@"
		@$(CXX)  -c -g $(CXXFLAGS) $<
		@mv *.o $(OBJDIR)

-include		$(DEPFILE)

.PHONY: clean
clean:
	@rm -rf $(MAINFILE) $(MAINOBJFILES)
	@echo "-> cleaning bin/"
	@echo "-> cleaning obj/"

#
