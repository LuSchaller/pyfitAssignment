SHELL := /usr/bin/bash
pyfitAssignment:
	c++ -ggdb3 -O3 -Wall -shared -std=c++17 -fPIC $$(python -m pybind11 --includes) pyfitAssignment.cc -o pyfitAssignment$$(python-config --extension-suffix)  -L/usr/lib64 -L/usr/lib64/root -L/cvmfs/cms.cern.ch/slc7_amd64_gcc530/cms/cmssw/CMSSW_8_0_27/lib/slc7_amd64_gcc530 -L/cvmfs/cms.cern.ch/slc7_amd64_gcc530/cms/cmssw/CMSSW_8_0_27/external/slc7_amd64_gcc530/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame -pthread -lm -ldl -rdynamic -lPhysicsToolsKinFitter -ltbb



