.PHONY: all clean

.DEFAULT_GOAL = all

PROJECT = barnes

C_SOURCES = $(shell find . -name "*.c")
CC_SOURCES = $(shell find . -name "*.cc")

VERBOSE = @ 
OBJDIR = ./build
DEPDIR = ./dep
#CXX = clang++
#CC = clang
CXX = g++
CC = gcc

WARNFLAGS = -Wall -Wextra
COMMONFLAGS = -fopenmp `sdl-config --cflags`
OPTFLAGS = -O3
#OPTFLAGS = -O0 -ggdb
LDFLAGS = `sdl-config --cflags --libs` -lGL -lglut -lGLU
CXXFLAGS = $(OPTFLAGS) $(WARNFLAGS) $(COMMONFLAGS)
CCFLAGS = $(OPTFLAGS) $(WARNFLAGS) $(COMMONFLAGS)

VPATH = $(sort $(dir $(C_SOURCES) $(CC_SOURCES)))

OBJECTS = $(notdir $(CC_SOURCES:.cc=.o)) $(notdir $(C_SOURCES:.c=.o))
DEP_FILES = $(patsubst %.o,$(DEPDIR)/%.d,$(OBJECTS))
OBJPRE = $(addprefix $(OBJDIR)/,$(OBJECTS))
MAKEFILE_LIST = Makefile

all: $(OBJDIR)/$(PROJECT)

$(OBJDIR)/$(PROJECT): $(OBJPRE) $(MAKEFILE_LIST)
	@echo "LD		$@"
	@if test \( ! \( -d $(@D) \) \) ;then mkdir -p $(@D);fi
	$(VERBOSE) $(CXX) -o $(OBJDIR)/$(PROJECT) $(COMMONFLAGS) $(LDFLAGS) $(OBJPRE)

$(DEPDIR)/%.d : %.c $(MAKEFILE_LIST)
	@echo "dep		$@"
	@if test \( ! \( -d $(@D) \) \) ;then mkdir -p $(@D);fi
	$(VERBOSE) $(CC) $(CCFLAGS) -MM -MT $(OBJDIR)/$*.o -MF $@ $<

$(DEPDIR)/%.d : %.cc $(MAKEFILE_LIST)
	@echo "dep		$@"
	@if test \( ! \( -d $(@D) \) \) ;then mkdir -p $(@D);fi
	$(VERBOSE) $(CXX) $(CXXFLAGS) -MM -MT $(OBJDIR)/$*.o -MF $@ $<

$(OBJDIR)/%.o : %.c $(MAKEFILE_LIST)
	@echo "cc		$@"
	@if test \( ! \( -d $(@D) \) \) ;then mkdir -p $(@D);fi
	$(VERBOSE) $(CC) -c $(CCFLAGS) -o $@ $<

$(OBJDIR)/%.o : %.cc $(MAKEFILE_LIST)
	@echo "cxx		$@"
	@if test \( ! \( -d $(@D) \) \) ;then mkdir -p $(@D);fi
	$(VERBOSE) $(CXX) -c $(CXXFLAGS) -o $@ $<

clean:	
	@echo "rm		$(OBJDIR)"
	$(VERBOSE) rm -rf $(OBJDIR)
	@echo "rm		$(DEPDIR)"
	$(VERBOSE) rm -rf $(DEPDIR)

ifneq ($(MAKECMDGOALS),clean)
-include $(DEP_FILES)
endif

