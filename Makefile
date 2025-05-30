#***************************************************************************
#       user configuration
# set mpi = yes for the parallel version
# set debug = yes to debug
DEBUG ?= no
PROFILE ?= no
# end user configuration;
#***************************************************************************

#SHELL = /bin/sh

CC = g++

ifeq ($(strip $(DEBUG)), yes)
    CFLAGS += -ggdb -DDEBUG -D THPOOL_DEBUG 
else
    CFLAGS += -O3 -march=native #-fopenmp
endif

ifeq ($(strip $(PROFILE)), yes)
    PFLAGS += -g -pg -O0 
endif
    
CFLAGS += -Wall -m64 -fno-math-errno  -ffast-math #-fsignaling-nans


HTSINC?=-I/usr/local/include/htslib #-I/usr/local/include/eigen3/ -I/usr/include/boost/math/distributions
HTSLIB?=/usr/local/lib  #-L/usr/local/boost/lib 
CFLAGS += $(HTSINC) -g
LDFLAGS = -lm -lhts -lz -llzma -lbz2 #-llapack -lblas -lboost_iostreams
#LDFLAGS += -lpthread 


OBJS = poe.o 
#OBJS = poe.o thpool.o

all:	poe 

.cpp.o : ;  $(CC) $(CFLAGS) $(PFLAGS) -c $<

$(OBJS): %.o: %.cpp

poe: poe.o; $(CC) $(OBJS) $(LDFLAGS) $(PFLAGS) -o poe
#poe: poe.o thpool.o; $(CC) $(OBJS) $(LDFLAGS) $(PFLAGS) -o idul

#idul.o: idul.cpp 
#thpool.o: thpool.cpp thpool.h


clean:
	rm -f poe poe.o 

