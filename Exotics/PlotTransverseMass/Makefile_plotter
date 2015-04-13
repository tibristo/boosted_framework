SHELL = /bin/sh
UNAME = $(shell uname)

# Root variables
ROOTCFLAGS   = -L$(ROOTSYS)/lib $(shell root-config --cflags)
ROOTLIBS     = -L$(ROOTSYS)/lib $(shell root-config --libs)
ROOTGLIBS    = -L$(ROOTSYS)/lib $(shell root-config --glibs)

################
##linux
CXX          = g++
RM           = rm -f
AR           = ar cru

##Flags
CXXFLAGS     = -O0 -g -fPIC -funroll-loops -std=c++11 #-Wall 

EXENAME		= plotter
SRCEXT   	= cpp
SRCDIR  	= src
INCDIR   	= include
OBJDIR   	= build
EXEDIR  	= bin
SRCS    	:= src/plotter.cpp
OBJS    	:= build/plotter.o

GARBAGE  = $(OBJDIR)/plotter.o $(EXEDIR)/$(EXENAME)
INCDIR+=-I/home/tim/Programs/fastjet-install/include
#################
##Dependencies
# Linux
ifeq "$(UNAME)" "Linux"
RANLIB       = ranlib
CXXFLAGS    += -I$(INCDIR) $(ROOTCFLAGS) #-I$(GSLINC)
LINKFLAGS    = -g $(shell root-config --nonew) $(shell root-config --ldflags) -Wl,--no-as-needed
endif

# OS X
ifeq "$(UNAME)" "Darwin"
RANLIB       = ranlib
CXXFLAGS    += -I$(INCDIR) $(ROOTCFLAGS) `/home/tim/Programs/fastjet-install/bin/fastjet-config --cxxflags`
#
#-I$(GSLINC)
LINKFLAGS    =
endif

##Libraries
LIBS       += $(ROOTLIBS) -lThread `/home/tim/Programs/fastjet-install/bin/fastjet-config --libs`

##Targets
all : $(EXEDIR)/$(EXENAME)

$(EXEDIR)/$(EXENAME) : $(OBJS)
	$(CXX) -o $@ $(OBJS) $(LINKFLAGS) $(LIBS)

$(OBJDIR)/plotter.o : $(SRCDIR)/plotter.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean   :
	$(RM) $(GARBAGE)

cleanall:
	$(RM) $(GARBAGE)


##fastjet
#   g++ short-example.cc -o short-example \
      `fastjet-install/bin/fastjet-config --cxxflags --libs --plugins`
