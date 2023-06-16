.PHONY: RivetFlav.so
RivetFlav.so:
	rivet-build RivetFlav.so FlavAlgAnalysis.cc -I.. -lIFNPlugin -lCMPPlugin -lGHSAlgo -lSDFlavPlugin -L../IFNPlugin/ -L../CMPPlugin/ -L../GHSAlgo/ -L../SDFlavPlugin/
