.PHONY: all clean

all: RivetFlav.so RivetFlavLHCb.so

RivetFlav.so:
	rivet-build RivetFlav.so FlavAlgAnalysis.cc -Dhepmc3 -I.. -lIFNPlugin -lCMPPlugin -lGHSAlgo -lSDFlavPlugin -L../IFNPlugin/ -L../CMPPlugin/ -L../GHSAlgo/ -L../SDFlavPlugin/

RivetFlavLHCb.so:
	rivet-build RivetFlavLHCb.so FlavAlgAnalysisLHCb.cc -Dhepmc3 -I.. -lIFNPlugin -lCMPPlugin -lGHSAlgo -lSDFlavPlugin -L../IFNPlugin/ -L../CMPPlugin/ -L../GHSAlgo/ -L../SDFlavPlugin/

clean:
	rm -rf *.so
