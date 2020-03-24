ifeq (,$(CC_LOCAL))
CC           := gcc
else
CC           := $(CC_LOCAL)
endif
ifeq (,$(CXX_LOCAL))
CXX           := g++
else
CXX           := $(CXX_LOCAL)
endif

CPPFLAGS+=$(shell root-config --cflags) -I.
#LDFLAGS+=$(shell root-config --libs) -Wl,-rpath,$(shell root-config --libdir)
##LDFLAGS+=$(shell root-config --libs --glibs) -lGui -lMinuit -lNew -Wl,-rpath -Wl,$(shell root-config --libdir)
LDFLAGS+=$(shell root-config --libs --glibs) -lGui -lMinuit -Wl,-rpath -Wl,$(shell root-config --libdir)

.PHONY: all clean

all: eff_single_cnt eff_cnt_all eff_single_cnt_endcap.C 

%.o: %.cc
	$(CXX) -c $(CPPFLAGS) $^
        
eff_cnt_single: eff_cnt_single.o
	$(CXX) -o $@ $(LDFLAGS) $^
eff_cnt_all: eff_cnt_all.o
	$(CXX) -o $@ $(LDFLAGS) $^
eff_single_cnt_endcap: eff_single_cnt_endcap.o
	$(CXX) -o $@ $(LDFLAGS) $^
clean:
	rm -f eff_single_cnt.o eff_single_cnt eff_cnt_all.o eff_cnt_all eff_single_cnt_endcap.o eff_single_cnt_endcap
       
