# Makefile

SRCDIR	=	./src
OBJDIR	=	./obj
BINDIR	=	./bin

MAINOBJ		=	main.o read_data.o qpsolver.o vector.o \
					probdata.o svpsolver.o solution.o svps_solve.o \
					heur_init.o solution_pool.o node.o stopwatch.o \
					select_node.o relax.o branch.o svps_parallel.o \
					heur.o heur_unitsphere.o Schmidt_manager.o \
					heur_quadratic.o cut.o cut_pool.o cut_oa.o
MAINSRC		=	$(addprefix $(SRCDIR)/,$(MAINOBJ:.o=.cpp))
MAINOBJFILES	=	$(addprefix $(OBJDIR)/,$(MAINOBJ))

MAIN		=	SVPSOLVER
MAINFILE	=	$(BINDIR)/$(MAIN)

#Mac
CXX	=	/usr/local/bin/g++-4.9

#Mac
CXXFLAGS	=	-Wall -lblas -llapack -lpthread -fopenmp

#Linux
#CXXFLAGS	=	-Wall -lopenblas -llapack -lpthread -fopenmp

CXXFLAGS	+=	-DNDEBUG

#Mac
CXXFLAGS	+=	-I/usr/local/Cellar/openblas/0.2.18_2/include

#Linux
#CXXLDFLAGS	=	-L/usr/lib64/libopenblas.so

.PHONY: all
all: $(MAIN)

$(MAIN):	$(MAINOBJFILES)
		@$(CXX) -O2 -o $(MAINFILE) $(CXXLDFLAGS) $(MAINOBJFILES) $(CXXFLAGS)

$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp
		@echo "-> compiling $@"
		@$(CXX) $(CXXFLAGS) -c -g $< 
		@mv *.o $(OBJDIR)

.PHONY: clean
clean:
	@rm -rf $(MAINFILE) $(MAINOBJFILES)
	@echo "-> cleaning bin/"
	@echo "-> cleaning obj/"
#
