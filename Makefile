#
# Makefile for vdb-gpu
#

# ----- Make Macros -----

# OPTFLAGS  =   -O2
DEFINES   =
INCLUDES  =
STD = -std=c++11 -stdlib=libc++
CXXFLAGS  =	$(STD) -g -Wall -Wextra -pedantic $(DEFINES) $(OPTFLAGS) $(INCLUDES)
CXX	  =	clang++

TARGETS = testing
TESTING_TARGETS = testing.o

# ----- Make Rules -----

all:	$(TARGETS)

testing:    $(TESTING_TARGETS)
	$(CXX) $(CXXFLAGS) -o testing $(TESTING_TARGETS)

clean:
	rm -f $(TARGETS) *.o

# ------ Dependences (.cpp -> .o using default Makefile rule) -----
testing.o: testing.cpp coord.hpp utils.hpp vdb.hpp vdb-private.hpp
