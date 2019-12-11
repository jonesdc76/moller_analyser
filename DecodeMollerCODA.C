// Analysis of data in crate "mycrate" and bank "mybank".
// See the parameters "mycrate" and "mybank" below.
// R. Michaels, Oct 2016
// Don Jones-- Dec 2019 Added argument to analyze starting at event number
#define MAXLEN  200000
#define MAXROC    50
#define MAXBANK   20
#define MAXRAW  6000
#define MAXDAT   500
#define NADC_CH 24
#define NTDC_CH 32
#define NSCA_CH 32
#define NTRIG 8
#define ERRORFLAG -9999

#include <iostream>
#include <string>
#include <vector>
#include "THaCodaFile.h"
#include "THaEtClient.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TNtuple.h"
#include "TRandom.h"
#include <vector>
#include <TGClient.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TRandom.h>
#include <TGButton.h>
#include <TGFrame.h>
#include <TRootEmbeddedCanvas.h>
#include <RQ_OBJECT.h>

using namespace std;
Int_t parseADCevent(int *data, int start_pos, int mol_type);
Int_t parseScalerEvent(int *data, int start_pos);
void usage();
int decode(int* data);
void clear();

// Global data 

const int  MAXEVT = 100000000;
Int_t *irn, *rocpos, *roclen;
Int_t *bankpos, *banklen;
const Int_t myroc1=7, myroc2=16;
Int_t myroc=-1;
Int_t evlen, evnum;
Int_t rdata[MAXDAT];
Int_t Scal[NSCA_CH];
Short_t ADC[NADC_CH], TDC_ch[NTDC_CH], TDC_ed[NTDC_CH], TDC[NTDC_CH];
Short_t coda_data_type, coda_event_type, run = 0;
Int_t ret, nadc, ntdc, nsca, tick;
static Int_t ntrig = NTRIG;
Short_t trig[NTRIG];


int main(int argc, char* argv[])
{
  //  Int_t load_status = gSystem->Load("./UTEvent_minimal_h.so");
  //  cout << "Status: " << load_status << endl;
  if(!gSystem->Getenv("MOLLER_DATA_DIR")){
    cout<<"Set environment variable MOLLER_DATA_DIR\n"
	<<"to point at the directory where the binary\n"
	<<".dat files are then rerun."<<endl;
    exit(0);
  }
  if(!gSystem->Getenv("MOLLER_ROOTFILE_DIR")){
    cout<<"Set environment variable MOLLER_ROOTFILE_DIR\n"
	<<"to point at the directory where the .root files\n"
	<<"should be stored then rerun."<<endl;
    exit(0);
  }
  
  THaCodaData *coda;      

  Int_t maxevent = MAXEVT, start_event = 0;
  Int_t evlo = start_event, evhi = evlo + maxevent;

  irn = new Int_t[MAXROC];
  rocpos = new Int_t[MAXROC];
  roclen = new Int_t[MAXROC];
  TString filename, rootfname;
  TString default_name = "run.dat";//may be a link to CODA datafile on disk
  TFile *file;
  TTree *trSca, *trADC;


  
  //Process command line arguments
  ////////////////////////////////////////
  if (argc > 1){
    run = atoi(argv[1]);
    cout << "Processing run number "<<run<<endl;
    filename = Form("%s/moller_data_%i.dat",
		    gSystem->Getenv("MOLLER_DATA_DIR"), run);
    rootfname = Form("%s/moller_data_%i.root",
		     gSystem->Getenv("MOLLER_ROOTFILE_DIR"), run);
    cout << "Processing file "<<filename.Data()<<endl;
  }else{
    usage();
    exit(0);
  }

  if (argc > 2){
    maxevent = atoi(argv[2]);
    evhi = maxevent;
    cout << "Processing a maximum of "<<maxevent<<" events";
  }

  if (argc > 3){
    start_event = atoi(argv[3]);
    evlo = start_event;
    evhi += evlo;
    cout << " starting at event "<<start_event<<endl;
  }else{
    cout<<endl;
  }
    


  

  // Initialize root and output
  //////////////////////////////////////////
  file = new TFile(rootfname.Data(),"RECREATE","Hall A Moller data");
  coda = new THaCodaFile();
  if (coda->codaOpen(filename) != 0) {
    cout << "ERROR:  Cannot open CODA data" << endl;
    goto end1;
  }
    
  //Set up scaler tree
  //////////////////////////////////////////
  trSca = new TTree ("trSca","Moller scaler tree");
  trSca->Branch("irun", &run);
  trSca->Branch("idtype", &coda_data_type);
  trSca->Branch("ievtype", &coda_event_type);
  trSca->Branch("iret", &ret);
  trSca->Branch("ntrig", &ntrig);
  trSca->Branch("itrig", trig, "itrig[ntrig]/S");
  trSca->Branch("itick", &tick);
  trSca->Branch("nsca", &nsca);
  trSca->Branch("isca", Scal, "isca[nsca]/I");

  //Set up ADC/TDC tree
  //////////////////////////////////////////
  trADC = new TTree ("trADC","Moller ADC and TDC tree");
  trADC->Branch("irun", &run);
  trADC->Branch("idtype", &coda_data_type);
  trADC->Branch("ievtype", &coda_event_type);
  trADC->Branch("iret", &ret);
  trADC->Branch("ntrig", &ntrig);
  trADC->Branch("itrig", trig);
  trADC->Branch("itick", &tick);
  trADC->Branch("nadc", &nadc);
  trADC->Branch("iadc", ADC, "iadc[nadc]/S");
  trADC->Branch("ntdc", &ntdc);
  trADC->Branch("itim", TDC, "itim[ntdc]/S");
  trADC->Branch("itch", TDC_ch, "itch[ntdc]/S");
  trADC->Branch("ited", TDC_ed, "ited[ntdc]/S");

  
  // Loop over events
  //////////////////////////////////////////
  int status;

  for (int iev = 0; iev < evhi; iev++) {      
    if (iev > 0 && ((iev%1000)==0) ) printf("%d events\n",iev);

    clear();

    status = coda->codaRead();  

    if (status != 0) {  // EOF or no data
	  
      if ( status == -1)  {
	cout << "End of CODA file. Bye bye." << endl;
	evhi=iev;
      }else {
	cout << "ERROR: codaRread status = " << hex << status << endl;
      }
      goto end1;
  

    }else {   // have data ...
      if(iev >= evlo){
	int type = decode( coda->getEvBuffer() );//real work happens here
	if(type >= 35 && type <= 37){
	  trADC->Fill();
	}
	if(type == 32) {
	  trSca->Fill();
	}
      }
    }
  }
  
 end1:
  
  coda->codaClose();

  file->Write();
  file->Close();
  return 0;

}; //end of main function




void usage() {  
  cout << "Usage:  molana(int run, int n_events = "<<MAXEVT<<", "
       <<"int start_event = 0) " << endl;
}; 

void clear() {
  memset(rdata,0,MAXDAT*sizeof(Int_t));
}


int decode (int* data) {
//Decodes CODA binary files
//Can analyze multiple ROCs
///////////////////////////////
  int mol_type = 0;
  evlen = data[0] + 1;
  coda_event_type = Short_t(data[1]>>16);
  coda_data_type = Short_t((data[1]>>8) & 0xff);
  evnum = data[4];
  static int printout = 0;  // print the raw data
  static int debug = 0;     // debug the decoding

  if (printout) { // Event dump
    cout << "\n\n Event number " << dec << evnum;
    cout << " length " << evlen << " CODA data type " << coda_data_type
	 << " CODA event type " << coda_event_type << endl;
    
    int ipt = 0;
    for (int j = 0; j < (evlen/5); j++) {
      cout << dec << "\n evbuffer[" << ipt << "] = ";
      for (int k=j; k<j+5; k++) {
	cout << hex << data[ipt++] << " ";
      }
      cout << endl;
    }
    if (ipt < evlen) {
      cout << dec << "\n evbuffer[" << ipt << "] = ";
      for (int k=ipt; k<evlen; k++) {
	cout << hex << data[ipt++] << " ";
      }
      cout << endl;
    }
    
  }

  // Decoding for "physics" events

  if (coda_data_type == 16 && coda_event_type <15) {

    // First find pointers to ROCs.
    // Useful for multi-crate analysis.
    // For Moller only ROC=7 used.
    ///////////////////////////////////
    int pos = data[2]+3;
    int nroc = 0;
    while( pos+1 < evlen ) {
      int len  = data[pos];
      int iroc = (data[pos+1]&0xff0000)>>16;
      if(iroc>=MAXROC || nroc >= MAXROC){
	cout << "Decoder error:  ";
	cout << "  ROC num " <<dec<<iroc;
	cout << "  nroc "<<nroc<<endl;
	return -1;
      }
      // Save position and length of each found ROC data block
      rocpos[iroc]  = pos;
      roclen[iroc]  = len;
      irn[nroc++] = iroc;
      if(debug)
	std::cout<<"ROC No. "<<iroc<<" with data length "<<len<<std::endl;
      pos += len+1;
    }
    if(debug)
      std::cout<<nroc<<" ROCs found"<<std::endl;
    Int_t found = 0;
    for (int j=0; j<nroc; j++){
      if(myroc1 == irn[j]){
	found = 1; myroc = myroc1;
	break;
      }else if (myroc2 == irn[j]){
	found = 1; myroc = myroc1;
	break;
      }
    }
    if (!found) {
      cout << "ERROR:  ROCs "<<dec<<myroc1<<" and "<<myroc2
	   <<" not in datastream !!"<<endl;
      return -1;
    }

    if (debug) {
      cout << "Roc info "<<nroc<<endl;
      for (int j=0; j<nroc; j++) 
	{
	  Int_t iroc = irn[j];
	  cout << "irn "<<dec<<iroc<<"   pos "<<rocpos[iroc];
	  cout <<"   len "<<roclen[iroc]<<endl;
	}
    }

    // Go over the data in myroc

    pos = rocpos[myroc]+2;
    if(pos+2 >= evlen){
      cout<<"ERROR: event length incorrect!"<<endl;
      return -1;
    }
    mol_type = data[pos++];// Moller data type scaler=32, adc=35, 36 or 37

    if(mol_type >= 35 && mol_type <= 37){
      
      parseADCevent(data, pos, mol_type);

    }else if(mol_type == 32){
      
      parseScalerEvent(data, pos);
      
    }else{
      
      cout<<"Error. Neither scalar nor ADC: mol_type="<<mol_type<<endl;
      
    }
  }
  return mol_type;
}


Int_t parseADCevent(int *data, int start_pos, int mol_type){
  ntdc = nadc = 0;
  //Parse Moller ADC readout
  ////////////////////////////
  Bool_t debug = 0;
  Int_t pos = start_pos;
  Int_t nchan = data[pos++];   // number of channels to read
  if(nchan > NADC_CH){
    cout<<"ERROR: too many ADC hits"<<endl;
    exit(0);
  }
  for(int i=0;i<NADC_CH;++i)ADC[i] = ERRORFLAG;
  for(int i=0; i<nchan; ++i){
    if(pos >= evlen){
      cout<<"Error. End of event encountered while reading ADC."<<endl;
      exit(0);
    }
    ADC[i] = (Short_t)data[pos++];
    ++nadc;
  }
      
  //Parse Moller TDC readout
  ////////////////////////////
  for(int i=0;i<NADC_CH;++i){
    TDC[i] = ERRORFLAG; TDC_ed[i] = ERRORFLAG; TDC_ch[i] = ERRORFLAG;
  }

  nchan = data[pos++];
  if(nchan > NTDC_CH){
    cout<<"ERROR: too many TDC hits"<<endl;
    exit(0);
  }
  for(int i=0; i<nchan; ++i){
    if(pos >= evlen){
      cout<<"Error. End of event encountered while reading TDC."<<endl;
      exit(0);
    }
    Int_t tmp = data[pos++];
    TDC[i] = Short_t(tmp & 0xffff);
    tmp >>= 16;
    TDC_ed[i] = Short_t(tmp & 1);//TDC trigger phase.
    //Leading edge=1, Trailing edge=0
    tmp >>= 1;
    TDC_ch[i] = Short_t(tmp & 0x1f);
    ++ntdc;
  }
      
  //Parse 2nd Moller TDC readout if it exists
  ///////////////////////////////////////////
  if(mol_type == 37){
    nchan = data[pos++];
    if(nchan > NTDC_CH){
      cout<<"Error. Too many TDC hits."<<endl;
      exit(0);
    }
    for(int i=0; i<nchan; ++i){
      if(pos >= evlen){ 
	cout<<"Error. End of event encountered while reading TDC 2."<<endl;
	exit(0);
      }
      Int_t tmp = data[pos++];
      TDC[i] = Short_t(tmp & 0xffff);
      tmp >>= 16;
      TDC_ed[i] = Short_t(tmp & 1);//TDC trigger phase.
      //Leading edge=1, Trailing edge=0
      tmp >>= 1;
      TDC_ch[i] = Short_t((tmp & 0x1f) + 32);
    }
    ++ntdc;
  }


  //Parse Moller status record
  ////////////////////////////
  tick = ERRORFLAG, ret = ERRORFLAG;
  for(int i=0;i<ntrig;++i)trig[i] = ERRORFLAG;
  nchan = data[pos++];
  for(int i=0; i<nchan; ++i){
    if(pos >= evlen){
      cout<<"Error. End of event encountered while reading status record."
	  <<endl;
      if(debug){
	cout<<"nchan = "<<nchan<<" pos= "<<pos<<" data["<<pos<<"]= "<<data[pos]
	    <<endl;
      }
      cout<<endl;
      exit(0);
    }
    int tmp = data[pos++];
    if(i==0){
      for(int j=0;j<ntrig;++j)
	trig[j] = Bool_t((1<<j) & tmp);
    }
    else if(i==2){
      tick = tmp;
    }
    else if(i==4)
      ret = tmp;
  }
  return pos;
}

Int_t parseScalerEvent(int *data, int start_pos){
  nsca = 0;
  //Parse Moller scaler readout
  /////////////////////////////
  Int_t pos = start_pos;
  Int_t  nchan = data[pos++];
  for(int i=0;i<NSCA_CH;++i)Scal[i] = ERRORFLAG;
  if(nchan > NSCA_CH){
    cout<<"ERROR: too many scaler hits."<<endl;
    exit(0);
  }
  for(int i=0; i<nchan; ++i){
    if(pos+1 >= evlen){
      cout<<"Error. End of event encountered while reading scalers.";
      cout<<endl;
      exit(0);
    }
    Scal[i] = (Int_t)data[pos++];
    ++nsca;
  }


  
  //It appears that there is a 0 length ADC + 0 length TDC event + status record
  //appended to the end of each scaler event. Perhaps this is just the way CODA
  //forces the writeout of a status record. So let's parse an ADC event...
  int mol_type = data[pos++];
  return parseADCevent(data, pos, mol_type);
}
