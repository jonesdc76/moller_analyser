CODA decoder for Hall A Moller polarimeter data based on R. Michael's
simpleAna standalone code. -Don Jones, Dec. 2019

This folder contains all the files required to read the CODA binary and
interpret it. User manuals for the TDC and scalers are included. For the TDC in
particular this was critical in understanding the structure of the TDC data.

Documentation for decoding the CODA binary is included in the slide show DecodeMollerCODA.pdf(pptx).

Requires that two environment variables be set

    MOLLER_DATA_DIR directory where CODA binary .dat files live
   
    MOLLER_ROOTFILE_DIR directory where the created ROOT files are to be stored

MollerDecodeCODA.C contains main and compiles to the executable molana.

  Arguments: 
            
	    1. Either int run_number (number of the run used to construct the filename)
	       or char* full name of file excluding directory path. Program will first
	       see if you input an integer before trying to interpret input as a filename.
	    
	    2. int n_events (optional number of events to analyze)
	     
	    3. int start_event (optional event to start analysis)

newrun is a bash script that calls molana and requires at least a run number

   ./newrun 204 10 6 will analyze 10 events from run 204 starting at event 6
   ./newrun moller_data_17000.dat will analyze the full run 17000.

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

