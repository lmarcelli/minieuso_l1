#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <string.h>
#include <stdlib.h>
#include <Riostream.h>
#include <stddef.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TH1.h>
#include <TSystem.h>


UInt_t calc_timeoffset(UInt_t time_x){  
  UInt_t timeoff;
  UInt_t time_0 = 1565611200;     //unix_time corresponding to August 12, 2019, time 12:00 pm 
  Float_t dailyoff = 2.6;
  timeoff=((time_x-time_0)/86400)*dailyoff;   //timeoff: "dailyoff" seconds per day
  //printf("timeoff da funzione: %f   %i", dailyoff, timeoff);
  return timeoff;
}


Int_t  main(int argc, char *argv[]){ 
  
  gSystem->Load("libTree");
  
  //Get raw data file from command line
  if ( argc == 1 ) {
    printf( "No arguments passed. You have to pass 2 arguments: Level0_filename and tle_filename. A third argument is optional: cpu_time_offset to be added to timestamp_unix\n" ); 
  }
  else {
    TFile *oldfile = new TFile(argv[1]);
    TFile *tlefile = new TFile(argv[2]);
    UInt_t timeoff_arg = 0;
    if (argc==4) timeoff_arg = atoi(argv[3]);
    //printf("timeoff_arg:  %i", timeoff_arg);

    //Create a new file: l1 file
    int len = strlen(argv[1]);
    int corr_len= len-5;
    char outputfile[corr_len];
    strncpy(outputfile, argv[1], corr_len);
    //strncat(outputfile,"_l1_compiled.root", 17);
    strncat(outputfile,"_l1.root", 8);
    TFile *newfile = new TFile(outputfile,"recreate");
    
    typedef struct {
      Double_t abstime;
      Double_t timeunix;
      Double_t lat;
      Double_t lon;
      Double_t alt;
      Double_t x;
      Double_t y;
      Double_t z;
      Double_t vx;
      Double_t vy;
      Double_t vz;
      Double_t CdLon;
      Double_t CdLat;
      Double_t X_geo;				
      Double_t Y_geo;				
      Double_t Z_geo;							
    } positionDiss;
    positionDiss PosizioneISS;
    
    typedef struct {			
      Double_t Sun_from_ISS_isDay;		
      Double_t Sun_from_ISS_alt;		
      Double_t Sun_from_ISS_az;		
      Double_t Sun_from_ISS_dist;		
      Double_t Sun_from_Earth_isDay;		
      Double_t Sun_from_Earth_alt;		
      Double_t Sun_from_Earth_az;		
      Double_t Sun_from_Earth_dist;
    } positionDSun;
    positionDSun PosizioneSun;
    
    typedef struct {		
      Double_t Moon_from_ISS_isDay;		
      Double_t Moon_from_ISS_alt;		
      Double_t Moon_from_ISS_az;		
      Double_t Moon_from_ISS_dist;       	
      Double_t Moon_from_Earth_isDay;		
      Double_t Moon_from_Earth_alt;		
      Double_t Moon_from_Earth_az;		
      Double_t Moon_from_Earth_dist;
      Double_t FractionMoon;
    } positionDMoon;
    positionDMoon PosizioneMoon;
    
    
    //Clone TTrees that do not need to be corrected in the new file
    if(oldfile->Get("televent")){
      TTree *televent_old = (TTree*)oldfile->Get("televent");
      TTree *televent_new = televent_old->CloneTree();
    }
    if(oldfile->Get("thv")){
      TTree *thv_old = (TTree*)oldfile->Get("thv");
      TTree *thv_new = thv_old->CloneTree();
    }
    if(oldfile->Get("thk")){
      TTree *thk_old = (TTree*)oldfile->Get("thk");
      TTree *thk_new = thk_old->CloneTree();
    }
    if(oldfile->Get("texp")){
      TTree *texp_old = (TTree*)oldfile->Get("texp");
      TTree *texp_new = texp_old->CloneTree();
    }
    if(oldfile->Get("ttherm")){
      TTree *ttherm_old = (TTree*)oldfile->Get("ttherm");
      TTree *ttherm_new = ttherm_old->CloneTree();
    }
    
    
    ///////////////////////////////////////////////////////////
    //tevent
    ///////////////////////////////////////////////////////////
    TTree *tevent_old = (TTree*)oldfile->Get("tevent");
    
    tevent_old->SetBranchStatus("*",1); 
    tevent_old->SetBranchStatus("photon_count_data",0); 
    TTree *tevent_new = tevent_old->CloneTree(0); 
    tevent_new->CopyEntries(tevent_old);  // Here we copy the branches */
    
    tevent_old->SetBranchStatus("*",0); 
    tevent_old->SetBranchStatus("photon_count_data",1); 
    Long64_t nentries = tevent_old->GetEntries();
    UChar_t ph[1][1][48][48];
    tevent_old->SetBranchAddress("photon_count_data",&ph);
    TBranch *corr_branch = tevent_new->Branch("photon_count_data", &ph, "photon_count_data[1][1][48][48]/b");
    
    for (Long64_t i=0;i<nentries; i++) {
      tevent_old->GetEntry(i);
      for (Long64_t j=40; j<48; j++) {
	for (Long64_t k=32; k<40; k++) {
	  ph[0][0][j][k]=ph[0][0][j][k]>>1;
	}
      }
      for (Long64_t j=8; j<16; j++) {
	for (Long64_t k=8; k<16; k++) {
	  ph[0][0][j][k]=ph[0][0][j][k]>>1;
	}
      }
      corr_branch->Fill();
    }
    
    ///////////////////////////////////////////////////////////
    //tevent_1st_integral
    ///////////////////////////////////////////////////////////
    TTree *tevent_1st_integral_old = (TTree*)oldfile->Get("tevent_1st_integral");
    
    tevent_1st_integral_old->SetBranchStatus("*",1); 
    tevent_1st_integral_old->SetBranchStatus("photon_count_data",0); 
    TTree *tevent_1st_integral_new = tevent_1st_integral_old->CloneTree(0); 
    tevent_1st_integral_new->CopyEntries(tevent_1st_integral_old);  // Here we copy the branches */
    
    tevent_1st_integral_old->SetBranchStatus("*",0); 
    tevent_1st_integral_old->SetBranchStatus("photon_count_data",1); 
    Long64_t nentries_1st = tevent_1st_integral_old->GetEntries();
    Float_t ph_1st[1][1][48][48];
    tevent_1st_integral_old->SetBranchAddress("photon_count_data",&ph_1st);
    TBranch *corr_branch_1st = tevent_1st_integral_new->Branch("photon_count_data", &ph_1st, "photon_count_data[1][1][48][48]/F");
    
    for (Long64_t i=0;i<nentries_1st; i++) {
      tevent_1st_integral_old->GetEntry(i);
      for (Long64_t j=40; j<48; j++) {
	for (Long64_t k=32; k<40; k++) {
	  ph_1st[0][0][j][k]=ph_1st[0][0][j][k]/2;
	}
      }
      for (Long64_t j=8; j<16; j++) {
	for (Long64_t k=8; k<16; k++) {
	  ph_1st[0][0][j][k]=ph_1st[0][0][j][k]/2;
	}
      }
      corr_branch_1st->Fill();
    }
    
    
    ///////////////////////////////////////////////////////////
    //tevent_2nd_integral
    ///////////////////////////////////////////////////////////
    TTree *tle_tree = (TTree*)tlefile->Get("Position");
    Int_t tletotev = tle_tree->GetEntries();
    TBranch *b_posiss = tle_tree->GetBranch("positionISS");
    b_posiss->SetAddress(&PosizioneISS.abstime);
    TBranch *b_possun = tle_tree->GetBranch("positionSun");
    b_possun->SetAddress(&PosizioneSun.Sun_from_ISS_isDay);
    TBranch *b_posmoon = tle_tree->GetBranch("positionMoon");
    b_posmoon->SetAddress(&PosizioneMoon.Moon_from_ISS_isDay);
    Int_t tick;
    tle_tree->SetBranchAddress("positionInt",&tick);
    
    TTree *tevent_2nd_integral_old = (TTree*)oldfile->Get("tevent_2nd_integral");   
    tevent_2nd_integral_old->SetBranchStatus("*",1); 
    tevent_2nd_integral_old->SetBranchStatus("photon_count_data",0); 
    TTree *tevent_2nd_integral_new = tevent_2nd_integral_old->CloneTree(0); 
    tevent_2nd_integral_new->CopyEntries(tevent_2nd_integral_old);  // Here we copy the branches */
    
    tevent_2nd_integral_old->SetBranchStatus("*",0); 
    tevent_2nd_integral_old->SetBranchStatus("timestamp_unix",1);
    tevent_2nd_integral_old->SetBranchStatus("photon_count_data",1); 
    Long64_t nentries_2nd = tevent_2nd_integral_old->GetEntries();
    Float_t ph_2nd[1][1][48][48];
    tevent_2nd_integral_old->SetBranchAddress("photon_count_data",&ph_2nd);
    TBranch *corr_branch_2nd = tevent_2nd_integral_new->Branch("photon_count_data", &ph_2nd, "photon_count_data[1][1][48][48]/F");
    
    UInt_t timestamp_unix;
    tevent_2nd_integral_old->SetBranchAddress("timestamp_unix",&timestamp_unix);
    UInt_t time_off;  //time offset to add to timestamp_unix (in seconds)
    TBranch *time_o = tevent_2nd_integral_new->Branch("time_offset", &time_off, "time_offset/i");
    UInt_t timestamp_unix_off;  //time offset to add to timestamp_unix (in seconds)
    TBranch *timestamp_unix_o = tevent_2nd_integral_new->Branch("timestamp_unix_offset", &timestamp_unix_off, "timestamp_unix_offset/i");
    TBranch *b_tle_iss_2nd = tevent_2nd_integral_new->Branch("positionISS", &PosizioneISS.abstime, "abstime/D:timeunix/D:lat/D:lon/D:alt/D:x/D:y/D:z/D:vx/D:vy/D:vz/D:CdLon/D:CdLat/D:X_geo/D:Y_geo/D:Z_geo/D");
    TBranch *b_tle_sun_2nd = tevent_2nd_integral_new->Branch("positionSun", &PosizioneSun.Sun_from_ISS_isDay, "Sun_from_ISS_isDay/D:Sun_from_ISS_alt/D:Sun_from_ISS_az/D:Sun_from_ISS_dist/D:Sun_from_Earth_isDay/D:Sun_from_Earth_alt/D:Sun_from_Earth_az/D:Sun_from_Earth_dist/D");
    TBranch *b_tle_moon_2nd = tevent_2nd_integral_new->Branch("positionMoon", &PosizioneMoon.Moon_from_ISS_isDay, "Moon_from_ISS_isDay/D:Moon_from_ISS_alt/D:Moon_from_ISS_az/D:Moon_from_ISS_dist/D:Moon_from_Earth_isDay/D:Moon_from_Earth_alt/D:Moon_from_Earth_az/D:Moon_from_Earth_dist/D:FractionMoon/D");
    TBranch *tick_2nd = tevent_2nd_integral_new->Branch("positionInt", &tick, "tick/I");
    
    
    if(timeoff_arg!=0) time_off=timeoff_arg;
    else {
      tevent_2nd_integral_old->GetEntry(0);
      time_off=calc_timeoffset(timestamp_unix); 
      //printf("\n first_timestamp_unix - timeoff: %i  %i ", timestamp_unix,  time_off);
    }
    // printf("\n time_off:  %i \n", time_off);
    Int_t cont = 0;
    Float_t diff1=0;
    Float_t diff2=0;
    
    for (Long64_t i=0;i<nentries_2nd; i++) {
      tevent_2nd_integral_old->GetEntry(i);
      timestamp_unix_off=timestamp_unix-time_off;
      for (Long64_t j=cont;j<tletotev; j++) {
	b_posiss->GetEntry(j);
	b_possun->GetEntry(j);
	b_posmoon->GetEntry(j);
	tle_tree->GetEntry(j);
	if(timestamp_unix_off-PosizioneISS.timeunix<=1){
	  diff1=timestamp_unix_off-PosizioneISS.timeunix;
	  b_posiss->GetEntry(j+1);  
	  tle_tree->GetEntry(j+1); 
	  diff2=timestamp_unix_off-PosizioneISS.timeunix;
	  if(fabs(diff1)<fabs(diff2)){
	    b_posiss->GetEntry(j);  
	    tle_tree->GetEntry(j); 
	    cont=j;
	  }else{
	    b_possun->GetEntry(j+1);  
	    b_posmoon->GetEntry(j+1);	    
	    cont=j+1;
	  }
	  j=tletotev;
	} 
      } 
      diff1=0;
      diff2=0;
      time_o->Fill();
      timestamp_unix_o->Fill();
      b_tle_iss_2nd->Fill();
      b_tle_sun_2nd->Fill();
      b_tle_moon_2nd->Fill();
      tick=tick;
      tick_2nd->Fill();
      //printf("\n fillato: %i  %i  %i  %i ", i, timestamp_unix,  time_off, timestamp_unix_off);
    }    
    
    
    for (Long64_t i=0;i<nentries_2nd; i++) {
      tevent_2nd_integral_old->GetEntry(i);
      for (Long64_t j=40; j<48; j++) {
	for (Long64_t k=32; k<40; k++) {
	  ph_2nd[0][0][j][k]=ph_2nd[0][0][j][k]/2;
	}
      }
      for (Long64_t j=8; j<16; j++) {
	for (Long64_t k=8; k<16; k++) {
	  ph_2nd[0][0][j][k]=ph_2nd[0][0][j][k]/2;
	}
      }
      corr_branch_2nd->Fill();
    }
    
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    
    newfile->Write(); 
  }
  
  return 0;
  
}
    
