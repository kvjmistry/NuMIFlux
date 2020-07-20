CXX=g++

CXXFLAGS=-std=c++17 -Wall -g

CPPFLAGS=-I $(BOOST_INC) \
	 -I $(CANVAS_ROOT_IO_INC) \
	 -I $(CANVAS_INC) \
	 -I $(CETLIB_INC) \
	 -I $(CETLIB_EXCEPT_INC) \
	 -I $(FHICLCPP_INC) \
	 -I $(HEP_CONCURRENCY_INC) \
	 -I $(GALLERY_INC) \
	 -I $(LARCOREOBJ_INC) \
	 -I $(LARDATAOBJ_INC) \
	 -I $(NUSIMDATA_INC) \
	 -I $(LARSIM_INC) \
	 -I $(DK2NUDATA_INC) \
	 -I $(ROOT_INC)

LDFLAGS=$(shell root-config --libs) \
	-L $(CANVAS_ROOT_IO_LIB) -l canvas_root_io\
	-L $(CANVAS_LIB) -l canvas\
	-L $(CETLIB_LIB) -l cetlib \
	-L $(CETLIB_EXCEPT_LIB) -l cetlib_except \
	-L $(GALLERY_LIB) -l gallery \
	-L $(NUSIMDATA_LIB) -l nusimdata_SimulationBase \
	-L $(LARCOREOBJ_LIB) -l larcoreobj_SummaryData \
	-L $(LARDATAOBJ_LIB) -l lardataobj_RecoBase \
	-L $(LARSIM_LIB) -l larsim_EventWeight_Base \
	-L $(DK2NUDATA_LIB) -l  dk2nuTree


#UNAME := $(shell uname -s)

EXEC1=bin/makehist
EXEC2=bin/checkfiles
EXEC3=bin/makehist_v2

all: $(EXEC1) $(EXEC2) $(EXEC3)
	
$(EXEC1): src/geo/GeoVector.cxx src/geo/GeoAABox.cxx src/geo/GeoHalfLine.cxx src/geo/GeoLine.cxx src/geo/GeoLineSegment.cxx src/geo/GeoCone.cxx src/geo/GeoSphere.cxx src/geo/GeoTrajectory.cxx src/geo/GeoAlgo.cxx src/makehist.cpp src/functions_makehist.h
	@echo Building $(EXEC1)
	@mkdir -p bin
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^


$(EXEC2):  src/geo/GeoVector.cxx src/geo/GeoAABox.cxx src/geo/GeoHalfLine.cxx src/geo/GeoLine.cxx src/geo/GeoLineSegment.cxx src/geo/GeoCone.cxx src/geo/GeoSphere.cxx src/geo/GeoTrajectory.cxx src/geo/GeoAlgo.cxx src/checkfiles.cpp
	@echo Building $(EXEC2)
	@mkdir -p bin
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

$(EXEC3): src/makehist_v2.cpp src/functions_makehist_v2.h
	@echo Building $(EXEC3)
	@mkdir -p bin
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

clean:
	rm bin/makehist bin/checkfiles bin/makehist_v2
