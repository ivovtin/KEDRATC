ifeq (,$(CC_LOCAL))
CC           := gcc
else
CC           := $(CC_LOCAL)
endif
ifeq (,$(CXX_LOCAL))
##CXX           := g++-4.4
CXX           := g++
else
CXX           := $(CXX_LOCAL)
endif

CPPFLAGS+=$(shell root-config --cflags) -I.
#LDFLAGS+=$(shell root-config --libs) -Wl,-rpath,$(shell root-config --libdir)
LDFLAGS+=$(shell root-config --libs --glibs) -lGui -lMinuit -lNew -Wl,-rpath -Wl,$(shell root-config --libdir)

.PHONY: all clean

all: npe_run_single npe_run thick_cnt_npe_run

%.o: %.cc
	$(CXX) -c $(CPPFLAGS) $^       
npe_run_single: npe_run_single.o
	$(CXX) -o $@ $(LDFLAGS) $^
npe_run: npe_run.o
	$(CXX) -o $@ $(LDFLAGS) $^
thick_cnt_npe_run: thick_cnt_npe_run.o
	$(CXX) -o $@ $(LDFLAGS) $^
clean:
	rm -f npe_run_single.o npe_run_single npe_run.o npe_run thick_cnt_npe_run.o thick_cnt_npe_run
     
