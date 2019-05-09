CXX=g++

CXXFLAGS=-std=c++17 -Wall -g

CPPFLAGS=-I $(BOOST_INC) \
	 -I $(CANVAS_ROOT_IO_INC) \
	 -I $(CANVAS_INC) \
	 -I $(CETLIB_INC) \
	 -I $(CETLIB_EXCEPT_INC) \
	 -I $(FHICLCPP_INC) \
	 -I $(GALLERY_INC) \
	 -I $(LARCOREOBJ_INC) \
	 -I $(LARDATAOBJ_INC) \
	 -I $(NUSIMDATA_INC) \
	 -I $(LARSIM_INC) \
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
	-L $(LARSIM_LIB) -l larsim_EventWeight_Base


#UNAME := $(shell uname -s)

EXEC1=bin/makehist
EXEC2=bin/makehist_uboone
EXEC3=bin/makehist_parent
EXEC4=bin/makehist_uboone_2D
EXEC5=bin/checkfiles
EXEC6=bin/makehist_uboone_detweights
EXEC7=bin/makehist_uboone_detweights_2D
EXEC8=bin/makehist_nova

all: $(EXEC1) $(EXEC2) $(EXEC3) $(EXEC4) $(EXEC5) $(EXEC6) $(EXEC7) $(EXEC8)
	
$(EXEC1): src/geo/GeoVector.cxx src/geo/GeoAABox.cxx src/geo/GeoHalfLine.cxx src/geo/GeoLine.cxx src/geo/GeoLineSegment.cxx src/geo/GeoCone.cxx src/geo/GeoSphere.cxx src/geo/GeoTrajectory.cxx src/geo/GeoAlgo.cxx src/makehist.cpp src/functions_makehist.h
	@echo Building $(EXEC1)
	@mkdir -p bin
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

$(EXEC2): src/geo/GeoVector.cxx src/geo/GeoAABox.cxx src/geo/GeoHalfLine.cxx src/geo/GeoLine.cxx src/geo/GeoLineSegment.cxx src/geo/GeoCone.cxx src/geo/GeoSphere.cxx src/geo/GeoTrajectory.cxx src/geo/GeoAlgo.cxx src/makehist_uboone.cpp
	@echo Building $(EXEC2)
	@mkdir -p bin
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

$(EXEC3): src/geo/GeoVector.cxx src/geo/GeoAABox.cxx src/geo/GeoHalfLine.cxx src/geo/GeoLine.cxx src/geo/GeoLineSegment.cxx src/geo/GeoCone.cxx src/geo/GeoSphere.cxx src/geo/GeoTrajectory.cxx src/geo/GeoAlgo.cxx src/makehist_parent.cpp src/functions_makehist.h
	@echo Building $(EXEC3)
	@mkdir -p bin
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

$(EXEC4):  src/geo/GeoVector.cxx src/geo/GeoAABox.cxx src/geo/GeoHalfLine.cxx src/geo/GeoLine.cxx src/geo/GeoLineSegment.cxx src/geo/GeoCone.cxx src/geo/GeoSphere.cxx src/geo/GeoTrajectory.cxx src/geo/GeoAlgo.cxx src/makehist_uboone_2D.cpp
	@echo Building $(EXEC4)
	@mkdir -p bin
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

$(EXEC5):  src/geo/GeoVector.cxx src/geo/GeoAABox.cxx src/geo/GeoHalfLine.cxx src/geo/GeoLine.cxx src/geo/GeoLineSegment.cxx src/geo/GeoCone.cxx src/geo/GeoSphere.cxx src/geo/GeoTrajectory.cxx src/geo/GeoAlgo.cxx src/checkfiles.cpp
	@echo Building $(EXEC5)
	@mkdir -p bin
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

$(EXEC6): src/makehist_uboone_detweights.cpp
	@echo Building $(EXEC6)
	@mkdir -p bin
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

$(EXEC7): src/makehist_uboone_detweights_2D.cpp
	@echo Building $(EXEC7)
	@mkdir -p bin
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

$(EXEC8): src/geo/GeoVector.cxx src/geo/GeoAABox.cxx src/geo/GeoHalfLine.cxx src/geo/GeoLine.cxx src/geo/GeoLineSegment.cxx src/geo/GeoCone.cxx src/geo/GeoSphere.cxx src/geo/GeoTrajectory.cxx src/geo/GeoAlgo.cxx src/makehist_nova.cpp src/functions_makehist.h
	@echo Building $(EXEC8)
	@mkdir -p bin
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

clean:
	rm bin/makehist bin/makehist_nova bin/makehist_uboone bin/makehist_uboone_parent bin/makehist_uboone_2D bin/checkfiles bin/makehist_uboone_detweights bin/makehist_uboone_detweights_2D
