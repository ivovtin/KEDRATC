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
##LDFLAGS+=$(shell root-config --libs --glibs) -lGui -lMinuit -lNew -Wl,-rpath -Wl,$(shell root-config --libdir)
LDFLAGS+=$(shell root-config --libs --glibs) -lGui -lMinuit -Wl,-rpath -Wl,$(shell root-config --libdir)

.PHONY: all clean

all: light_collection light_col_coef_cnt  

%.o: %.cc
	$(CXX) -c $(CPPFLAGS) $^
        
light_collection: light_collection.o
	$(CXX) -o $@ $(LDFLAGS) $^
light_col_coef_cnt: light_col_coef_cnt.o
	$(CXX) -o $@ $(LDFLAGS) $^       
clean:
	rm -f light_collection.o light_collection light_col_coef_cnt.o light_col_coef_cnt
       
