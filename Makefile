# Makefile assumes PREFIX (root path for installing files) to be set.
# If it is not set, it will point to /usr/.
# The path to TMV is taken from TMV_PREFIX, which defaults to PREFIX.
# Additional compilation flags can be specified with SPECIAL_FLAGS, 
# linker flags with SPECIAL_LIBS

LIBNAME = shapelens
INCLPATH = include
SRCPATH = src
LIBPATH = lib
DOCPATH = doc
PROGSRCPATH = progs
PROGPATH = bin
TESTSRCPATH = test

########
# setting defaults and other environment variables
########
ifndef PREFIX
	PREFIX = /usr
endif

ifndef TMV_PREFIX
	TMV_PREFIX = $(PREFIX)
endif

# which OS
UNAME := $(shell uname)

# compilation flags
CC = g++
CFLAGS = -ansi $(SPECIAL_FLAGS)
ifneq ($(UNAME),Linux)
	CFLAGS += -bind_at_load
endif

# linker flags
TMVLINK := $(shell cat ${TMV_PREFIX}/share/tmv/tmv-link)
LIBS = $(TMVLINK) -lcfitsio

ifneq (,$(findstring HAS_WCSLIB,$(SPECIAL_FLAGS)))
	LIBS += -lwcs
endif

# archiver flags
AR = ar
ARFLAGS = -sr
ifeq ($(UNAME),Linux)
	SHAREDFLAGS = -shared -fPIC 
	LIBEXT = so
else
	SHAREDFLAGS = -dynamiclib -fPIC
	LIBEXT = dylib
endif

########
# target section
########
SRC = $(wildcard $(SRCPATH)/*.cc)
OBJECTS = $(SRC:$(SRCPATH)/%.cc=$(LIBPATH)/%.o)
PROGS = $(wildcard $(PROGSRCPATH)/*.cc)
PROGSOBJECTS = $(PROGS:$(PROGSRCPATH)/%.cc=$(PROGPATH)/%)
TESTS = $(wildcard $(TESTSRCPATH)/*.cc)
TESTOBJECTS = $(TESTS:$(TESTSRCPATH)/%.cc=$(TESTSRCPATH)/%)
HEADERS = $(wildcard $(INCLPATH)/*.h)

all: $(LIBPATH) $(DOCPATH) library shared

$(LIBPATH):
	mkdir -p $(LIBPATH)

$(DOCPATH):
	mkdir -p $(DOCPATH)

.PHONY: clean

clean: 
	rm -f $(OBJECTS)
	rm -f $(LIBPATH)/lib$(LIBNAME).a
	rm -f $(LIBPATH)/lib$(LIBNAME).$(LIBEXT)	

cleandocs:
	rm -rf $(DOCPATH)/*

cleanprogs:
	rm -f $(PROGSOBJECTS)

library: $(LIBPATH)/lib$(LIBNAME).a

shared: $(LIBPATH)/lib$(LIBNAME).$(LIBEXT)

install: library shared
	mkdir -p $(PREFIX)/lib
	cp $(LIBPATH)/lib$(LIBNAME).a $(PREFIX)/lib/
	cp $(LIBPATH)/lib$(LIBNAME).$(LIBEXT) $(PREFIX)/lib/
	mkdir  -p $(PREFIX)/include/$(LIBNAME)
	cd $(INCLPATH) && find . -type f -name '*.h' -exec  cp --parents {} $(PREFIX)/include/$(LIBNAME)/ \; && cd ../

progs: $(PROGSOBJECTS)

installprogs: progs
	mkdir -p $(PREFIX)/bin
	cp $(PROGSOBJECTS) $(PREFIX)/bin/

test: $(TESTOBJECTS)

docs: $(HEADERS)
	doxygen doc/Doxyfile


$(LIBPATH)/lib$(LIBNAME).a: $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $?

$(LIBPATH)/lib$(LIBNAME).$(LIBEXT): $(OBJECTS)
	rm -f $(LIBPATH)/lib$(LIBNAME).$(LIBEXT)	
ifeq ($(UNAME),Linux)
	$(CC) $(SHAREDFLAGS) -o $@ $^
else
	$(CC) $(SHAREDFLAGS) $(SPECIAL_LIBS) -o $@ $^ $(LIBS)
endif

$(LIBPATH)/%.o: $(SRCPATH)/%.cc
	$(CC) $(CFLAGS) -c $< -o $@

$(PROGPATH)/%: $(PROGSRCPATH)/%.cc
	$(CC) $(CFLAGS) $(SPECIAL_LIBS) -L$(LIBPATH) $< -o $@ -l$(LIBNAME) $(LIBS)

$(TESTSRCPATH)/%: $(TESTSRCPATH)/%.cc
	$(CC) $(CFLAGS) $(SPECIAL_LIBS) -L$(LIBPATH) $< -o $@ -l$(LIBNAME) $(LIBS)
