
#compiler
cxx = g++

# executable name
program = run

# source files
srcs = main.cpp invariants.cpp

# directory to sources and compiled obj files
srcdir = src
objdir = obj
#pchdir = pch

objs = $(srcs:.cpp=.o)
bins = $(src:.cpp=)

# global cpp flags
cppflags = -std=c++17 -O3 -fopenmp

# specialized cpp flags - per compile basis
cppflags += $(cppflags-$@)
# flags needed only to compile executable at final step
#cppflags-run += -lGLU -lGL -lglut -lboost_program_options -lsfml-graphics -lsfml-system -lsfml-window
cppflags-run += -lmpfr -lgmp -lGLU -lGL -lglut -fopenmp
# linker flags:
#linkflags = -lGLU -lGL -lglut -lboost_program_options


# giving folder prefix to the executable files.
exeobjs = $(addprefix $(objdir)/, $(objs))


# recipe to make
all: $(program)

$(program): $(exeobjs)
	$(cxx) -o $(program) $(exeobjs) $(cppflags) $(linkflags)

$(objdir)/%.o: $(srcdir)/%.cpp | objdirmk
	$(cxx) $(cppflags) -c $< -o $@


# check to see if obj/ is already a directory. if not, make it.
objdirmk:
	@mkdir -p $(objdir)

.PHONY: clean
# removes the executable and the compiled obj .o files
clean:
	@echo Removing executable and cleaning object directories.
	$(RM) $(program) $(EXEOBJS) 
	$(RM) -r $(objdir)


