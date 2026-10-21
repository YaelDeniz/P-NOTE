#MACROS:  Expression which is likely to be used repeatedl
CC= g++

# Source and build output
SRC_DIR = ./src
MAIN_DIR= ./main
OBJ_DIR = ./obj
BIN_DIR = .

# ROOT: resolved at build time via root-config (root-config must be on PATH)
ROOT_LIB    := $(shell root-config --libdir)
ROOT_INCDIR := $(shell root-config --incdir)

# OscProb: override from command line or env:  make OSCPROB_PATH=/your/path
OSCPROB_PATH          ?= $(HOME)/OscProb
OSCPROB_LIB            = $(OSCPROB_PATH)/lib
OSCPROB_INCLUDE        = $(OSCPROB_PATH)/inc
OSCPROB_EIGEN_INCLUDE  = $(OSCPROB_PATH)/eigen

export OSCPROB_PATH

# This project's own headers (relative to repo root, works on any machine)
PNOTE_INCLUDE = $(CURDIR)/inc


#Flags
ROOT_FLAGS = $(shell root-config --cflags --glibs --libs)
CPPFLAGS = -L$(ROOT_LIB) -I$(ROOT_INCDIR) -I$(OSCPROB_INCLUDE) -I$(OSCPROB_EIGEN_INCLUDE) -I$(PNOTE_INCLUDE)
CFLAGS   := -Wall
LDFLAGS  := -L$(OSCPROB_LIB) -Wl,-rpath,$(OSCPROB_LIB) -L$(ROOT_LIB)
# Full libs for targets that use MathMore (requires ROOT built with mathmore=ON + libgsl-dev)
LDLIBS        := -lOscProb -lMathMore -lGeom
# Libs for Earth3D — does not use MathMore, avoids the missing-library link error
LDLIBS_NOMATH := -lOscProb -lGeom


#OBJECTS = $(SOURCES:.c=.o)

OBJECTS_AsivTrue = MathTools.o PhyTools.o Earth3DModel.o AsimovDataTrue.o 
OBJECTS_SensitivityTrue = MathTools.o PhyTools.o Earth3DModel.o AsimovDataTrue.o ExportData.o
OBJECTS_AsivObs  = MathTools.o PhyTools.o Earth3DModel.o AsimovDataTrue.o AsimovDataObs.o 
OBJECTS_SensitivityObs  = MathTools.o PhyTools.o Earth3DModel.o AsimovDataTrue.o AsimovDataObs.o ExportData.o

OBJECTS_Res  = MathTools.o PhyTools.o  DetectorResolution.o 
OBJECTS_OscProb  = MathTools.o PhyTools.o  OscProbEarth.o 
OBJECTS_Flux = PhyTools.o
OBJECTS_EARTHMODEL = Earth3DModel.o #OscProbEarth.o 


#OBJECTSProb = GetProbPREM.o


.PHONY: all clean

#TARGETS

#all: ObsN ObsSenv IntN IntSenv EarthSenv OscProbEarth HondaFlux Earth3D MatterN   #DetectorResolution   # All targets #Add StatsAnaTrue laters

all: Earth3D

MatterN: $(MAIN_DIR)/MatterEffect.cpp $(OBJECTS_AsivTrue) 
	@echo " "
	@echo "Generating Application for True Events"
	$(CC) -g $^ $(ROOT_FLAGS) $(LDFLAGS) $(LDLIBS) $(CPPFLAGS) -o $@
	@echo " "
	@echo "Matter effect example generator Done"

IntN: $(MAIN_DIR)/IntOscGram.cpp $(OBJECTS_AsivTrue) 
	@echo " "
	@echo "Generating Application for True Events"
	$(CC) -g $^ $(ROOT_FLAGS) $(LDFLAGS) $(LDLIBS) $(CPPFLAGS) -o $@
	@echo " "
	@echo "True event generator Done"

IntSenv: $(MAIN_DIR)/TrueSensitivity.cpp $(OBJECTS_SensitivityTrue) 
	@echo " "
	@echo "Generating application for true sensitivity simulations"
	$(CC) -g $^ $(ROOT_FLAGS) $(LDFLAGS) $(LDLIBS) $(CPPFLAGS) -o $@
	@echo " "
	@echo "Application for true sensitivity simulation ready"

ObsN: $(MAIN_DIR)/ObsOscGram.cpp $(OBJECTS_AsivObs) 
	@echo " "
	@echo "Generating Application for True Events"
	$(CC) -g $^ $(ROOT_FLAGS) $(LDFLAGS) $(LDLIBS) $(CPPFLAGS) -o $@
	@echo " "
	@echo "Observed event generator Done"

ObsSenv: $(MAIN_DIR)/ObsSensitivity.cpp $(OBJECTS_SensitivityObs) 
	@echo " "
	@echo "Generating application for observed sensitivity simulations"
	$(CC) -g $^ $(ROOT_FLAGS) $(LDFLAGS) $(LDLIBS) $(CPPFLAGS) -o $@
	@echo " "
	@echo "Application for observed sensitivity simulation ready"

EarthSenv: $(MAIN_DIR)/EarthTrueSensitivity.cpp $(OBJECTS_AsivTrue) 
	@echo " "
	@echo "Generating Application for True Events"
	$(CC) -g $^ $(ROOT_FLAGS) $(LDFLAGS) $(LDLIBS) $(CPPFLAGS) -o $@
	@echo " "
	@echo "True event generator Done"
#StatsAnaTrue: $(MAIN_DIR)/StatsAnaTrue.cpp $(OBJECTS_AsivTrue)
#	@echo " "
#	@echo "Generating Application for Statistical Analysis of True Events"
#	$(CC) -g $^ $(ROOT_FLAGS) $(LDFLAGS) $(LDLIBS) $(CPPFLAGS) -o $@
#	@echo " "
#	@echo "StatsTrue Done"

#StatsAnaObs: $(MAIN_DIR)/StatsAnaObs.cpp $(OBJECTS_AsivObs)
#	@echo " "
#	@echo "Generating Application for Statistical Analysis of Observed Events"
#	$(CC) -g $^ $(ROOT_FLAGS) $(LDFLAGS) $(LDLIBS) $(CPPFLAGS) -o $@
#	@echo " "
#	@echo "StatsObs Done"

OscProbEarth: $(MAIN_DIR)/OscProbEarth.cpp $(OBJECTS_OscProb)
	@echo " "
	@echo "Generating Application for Oscillation Probabilities"
	$(CC) -g $^ $(ROOT_FLAGS) $(LDFLAGS) $(LDLIBS) $(CPPFLAGS) -o $@
	@echo " "
	@echo "Probabilities Done"

HondaFlux: $(MAIN_DIR)/HondaFlux.cpp $(OBJECTS_Flux) 
	@echo " "
	@echo "Generating Application for Neutrino Flux visualization"
	$(CC) -g $^ $(ROOT_FLAGS) $(LDFLAGS) $(LDLIBS) $(CPPFLAGS) -o $@
	@echo " "
	@echo "Neutrino Flux vis generator Done"

Earth3D: $(MAIN_DIR)/Earth3DVis.cpp $(OBJECTS_EARTHMODEL)
	@echo " "
	@echo "Generating Application for Model visualization"
	$(CC) -g $^ $(ROOT_FLAGS) $(LDFLAGS) $(LDLIBS_NOMATH) $(CPPFLAGS) -o $@
	@echo " "
	@echo "3D Earth Modelvis generator Done"

#Build Objects

AsimovDataTrue.o: $(SRC_DIR)/AsimovDataTrue.C MathTools.o PhyTools.o 
	@echo "  "
	@echo "Creating Asimov True Events Objects"
	$(CC) -g -c $^ $(ROOT_FLAGS) $(LDFLAGS) $(LDLIBS) $(CPPFLAGS) -o $@


AsimovDataObs.o: $(SRC_DIR)/AsimovDataObs.C MathTools.o PhyTools.o 
	@echo "  "
	@echo "Creating Asimov Observed Events Objects"
	$(CC) -g -c $^ $(ROOT_FLAGS) $(LDFLAGS) $(LDLIBS) $(CPPFLAGS) -o $@

OscProbEarth.o: $(SRC_DIR)/OscProbEarth.C MathTools.o PhyTools.o
	@echo "  "
	@echo "Creating Oscillation Probabilities  Objects"
	$(CC) -g -c $^ $(ROOT_FLAGS) $(LDFLAGS) $(LDLIBS) $(CPPFLAGS) -o $@

MathTools.o: $(SRC_DIR)/MathTools.C
	@echo "  "
	@echo "Creating MathTools Objects"
	$(CC) -c -g $^ $(ROOT_FLAGS) $(CPPFLAGS) -o $@


PhyTools.o: $(SRC_DIR)/PhyTools.C
	@echo "  "
	@echo "Creating PhyTools Objects"
	$(CC) -c -g $^ $(ROOT_FLAGS) $(CPPFLAGS) -o $@
	

Earth3DModel.o: $(SRC_DIR)/Earth3DModel.C
	@echo "  "
	@echo "Creating a 3D Earth from a Geological models Objects"
	$(CC) -c -g $^ $(ROOT_FLAGS) $(CPPFLAGS) -o $@
	
ExportData.o: $(SRC_DIR)/ExportData.C
	@echo "  "
	@echo "Creating Data exportation tools"
	$(CC) -c -g $^ $(ROOT_FLAGS) $(CPPFLAGS) -o $@




clean:
	@echo "Removing Objects"
	rm -f *.o ObsN ObsSenv IntN IntSenv EarthSenv OscProbEarth HondaFlux Earth3D MatterN
