CODA decoder for Hall A Moller polarimeter data based on R. Michael's
simpleAna standalone code. -Don Jones, Dec. 2019

This folder contains all the files required to read the CODA binary and
interpret it. Requires that two environment variables be set
   MOLLER_DATA_DIR directory where CODA binary .dat files live
   MOLLER_ROOTFILE_DIR directory where the created ROOT files are to be stored

MollerDecodeCODA.C contains main and compiles to the executable molana.
  Arguments: int run_number (number of the run used to construct the filename)
             int n_events (optional number of events to analyze)
	     int start_event (optional event to start analysis)

newrun is a bash script that calls molana and requires at least a run number
  ./newrun 204 10 6 will analyze 10 events from run 204 starting at event 6

Makefile
    to compile, type "make"

Classes required to build the analyzer.
   et.h
   evio.C
   evio.h
   swap_util.C
   THaCodaData.C
   THaCodaData.h
   THaCodaFile.C
   THaCodaFile.h
   THaEtClient.C
   THaEtClient.h

