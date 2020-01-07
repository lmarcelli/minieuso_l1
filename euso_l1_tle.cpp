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


Float_t calc_timeoffset(UInt_t time_x){  
  Double_t timeoff;
  UInt_t time_0 = 1565611200;     //unix_time corresponding to August 12, 2019, time 12:00 pm 
  Float_t dailyoff = 2.6;
  timeoff=((time_x-time_0)/86400.)*dailyoff;   //timeoff: "dailyoff" seconds per day
  printf("timeoff da funzione: % i  %f   %f", time_x,  dailyoff, timeoff);
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
    
    
    //Clone TTrees that do not need to be corrected in the new file
    if(oldfile->Get("thk")){
      TTree *thk_old = (TTree*)oldfile->Get("thk");
      TTree *thk_new = thk_old->CloneTree();
    }
    if(oldfile->Get("texp")){
      TTree *texp_old = (TTree*)oldfile->Get("texp");
      TTree *texp_new = texp_old->CloneTree();
      Char_t l1_sw_vers[5]="v6.0"; 
      TBranch *l1_soft_vers_br = texp_new->Branch("l1_soft_vers", &l1_sw_vers, "l1_soft_vers/C");
      l1_soft_vers_br->Fill();
    }
    if(oldfile->Get("ttherm")){
      TTree *ttherm_old = (TTree*)oldfile->Get("ttherm");
      TTree *ttherm_new = ttherm_old->CloneTree();
    }
    

    Int_t cont = 0;
    Float_t diff1=0;
    Float_t diff2=0;
    Double_t time_off;  //time offset to add to timestamp_unix (in seconds)
    Double_t timestamp_unix_off;  //time offset to add to timestamp_unix (in seconds)
    Double_t gtu_time_off;  //time offset to add to timestamp_unix (in seconds)

    ///////////////////////////////////////////////////////////
    //tevent_2nd_integral
    ///////////////////////////////////////////////////////////
    TTree *tevent_2nd_integral_old = (TTree*)oldfile->Get("tevent_2nd_integral");   
    tevent_2nd_integral_old->SetBranchStatus("*",1); 
    tevent_2nd_integral_old->SetBranchStatus("photon_count_data",0); 
    TTree *tevent_2nd_integral_new = tevent_2nd_integral_old->CloneTree(0); 
    tevent_2nd_integral_new->CopyEntries(tevent_2nd_integral_old);  // Here we copy the branches */
    
    tevent_2nd_integral_old->SetBranchStatus("*",0); 
    tevent_2nd_integral_old->SetBranchStatus("timestamp_unix",1);
    tevent_2nd_integral_old->SetBranchStatus("gtu_time",1);
    tevent_2nd_integral_old->SetBranchStatus("photon_count_data",1); 
    Long64_t nentries_2nd = tevent_2nd_integral_old->GetEntries();
    printf("\n D3: nentries  %i  ", nentries_2nd);
    Float_t ph_2nd[1][1][48][48];
    tevent_2nd_integral_old->SetBranchAddress("photon_count_data",&ph_2nd);
    TBranch *corr_branch_2nd = tevent_2nd_integral_new->Branch("photon_count_data", &ph_2nd, "photon_count_data[1][1][48][48]/F");
    
    UInt_t timestamp_unix;
    tevent_2nd_integral_old->SetBranchAddress("timestamp_unix",&timestamp_unix);
    Double_t gtu_time;
    tevent_2nd_integral_old->SetBranchAddress("gtu_time",&gtu_time);
    TBranch *time_o_2nd = tevent_2nd_integral_new->Branch("time_offset", &time_off, "time_offset/D");
    TBranch *timestamp_unix_o_2nd = tevent_2nd_integral_new->Branch("timestamp_unix_offset", &timestamp_unix_off, "timestamp_unix_offset/D");
    TBranch *gtu_time_o_2nd = tevent_2nd_integral_new->Branch("gtu_time_offset", &gtu_time_off, "gtu_time_offset/D");
    TBranch *b_tle_iss_2nd = tevent_2nd_integral_new->Branch("positionISS", &PosizioneISS.abstime, "abstime/D:timeunix/D:lat/D:lon/D:alt/D:x/D:y/D:z/D:vx/D:vy/D:vz/D:CdLon/D:CdLat/D:X_geo/D:Y_geo/D:Z_geo/D");
    TBranch *b_tle_sun_2nd = tevent_2nd_integral_new->Branch("positionSun", &PosizioneSun.Sun_from_ISS_isDay, "Sun_from_ISS_isDay/D:Sun_from_ISS_alt/D:Sun_from_ISS_az/D:Sun_from_ISS_dist/D:Sun_from_Earth_isDay/D:Sun_from_Earth_alt/D:Sun_from_Earth_az/D:Sun_from_Earth_dist/D");
    TBranch *b_tle_moon_2nd = tevent_2nd_integral_new->Branch("positionMoon", &PosizioneMoon.Moon_from_ISS_isDay, "Moon_from_ISS_isDay/D:Moon_from_ISS_alt/D:Moon_from_ISS_az/D:Moon_from_ISS_dist/D:Moon_from_Earth_isDay/D:Moon_from_Earth_alt/D:Moon_from_Earth_az/D:Moon_from_Earth_dist/D:FractionMoon/D");
    TBranch *tick_2nd = tevent_2nd_integral_new->Branch("positionInt", &tick, "tick/I");
    
    time_off=0;
    if(nentries_2nd>0){
      if(timeoff_arg!=0) time_off=timeoff_arg;
      else {
    	tevent_2nd_integral_old->GetEntry(0);
    	time_off=calc_timeoffset(timestamp_unix); 
    	//printf("\n time_offset from D3: first_timestamp_unix - timeoff: %i  %f ", timestamp_unix,  time_off);
      }
      printf("\n time_off:  %f \n", time_off);
    }

    cont = 0;
    diff1=0;
    diff2=0;
    
    for (Long64_t i=0;i<nentries_2nd; i++) {
      tevent_2nd_integral_old->GetEntry(i);
      // if(timeoff_arg!=0) time_off=timeoff_arg;
      // else {
      // 	time_off=calc_timeoffset(timestamp_unix); 
      // }
      timestamp_unix_off=timestamp_unix-time_off;
      gtu_time_off=gtu_time-time_off;
      //printf("\n tevent_2nd: timestamp_unix -timeoff - timestamp_unix_off: %i    %f   %f", timestamp_unix,  time_off, timestamp_unix_off);
      //printf("\n tevent_2nd: gtu_time - timeoff - gtu_time_off: %f    %f   %f", gtu_time,  time_off, gtu_time_off);
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
      time_o_2nd->Fill();
      timestamp_unix_o_2nd->Fill();
      gtu_time_o_2nd->Fill();
      b_tle_iss_2nd->Fill();
      b_tle_sun_2nd->Fill();
      b_tle_moon_2nd->Fill();
      tick=tick;
      tick_2nd->Fill();
    }    
    
    
    for (Long64_t i=0;i<nentries_2nd; i++) {
      tevent_2nd_integral_old->GetEntry(i);
      for (Long64_t j=0; j<8; j++) {
	for (Long64_t k=32; k<40; k++) {
	  ph_2nd[0][0][j][k]=ph_2nd[0][0][j][k]/2;
	}
      }
      for (Long64_t j=32; j<40; j++) {
	for (Long64_t k=8; k<16; k++) {
	  ph_2nd[0][0][j][k]=ph_2nd[0][0][j][k]/2;
	}
      }
      corr_branch_2nd->Fill();
    }
    
    ///////////////////////////////////////////////////////////
    //thv
    ///////////////////////////////////////////////////////////
    if(oldfile->Get("thv")){
      TTree *thv_old = (TTree*)oldfile->Get("thv");
      TTree *thv_new = thv_old->CloneTree();
      
      thv_old->SetBranchStatus("timestamp_unix",1);
      Long64_t nentries_hv = thv_old->GetEntries();
      printf("\n hv: nentries  %i  ", nentries_hv);
      thv_old->SetBranchAddress("timestamp_unix",&timestamp_unix);
      TBranch *time_o_hv = thv_new->Branch("time_offset", &time_off, "time_offset/D");
      TBranch *timestamp_unix_o_hv = thv_new->Branch("timestamp_unix_offset", &timestamp_unix_off, "timestamp_unix_offset/D");
      TBranch *b_tle_iss_hv = thv_new->Branch("positionISS", &PosizioneISS.abstime, "abstime/D:timeunix/D:lat/D:lon/D:alt/D:x/D:y/D:z/D:vx/D:vy/D:vz/D:CdLon/D:CdLat/D:X_geo/D:Y_geo/D:Z_geo/D");
      TBranch *b_tle_sun_hv = thv_new->Branch("positionSun", &PosizioneSun.Sun_from_ISS_isDay, "Sun_from_ISS_isDay/D:Sun_from_ISS_alt/D:Sun_from_ISS_az/D:Sun_from_ISS_dist/D:Sun_from_Earth_isDay/D:Sun_from_Earth_alt/D:Sun_from_Earth_az/D:Sun_from_Earth_dist/D");
      TBranch *b_tle_moon_hv = thv_new->Branch("positionMoon", &PosizioneMoon.Moon_from_ISS_isDay, "Moon_from_ISS_isDay/D:Moon_from_ISS_alt/D:Moon_from_ISS_az/D:Moon_from_ISS_dist/D:Moon_from_Earth_isDay/D:Moon_from_Earth_alt/D:Moon_from_Earth_az/D:Moon_from_Earth_dist/D:FractionMoon/D");
      TBranch *tick_hv = thv_new->Branch("positionInt", &tick, "tick/I");
      
      cont = 0;
      diff1=0;
      diff2=0;

      for (Long64_t i=0;i<nentries_hv; i++) {
	thv_old->GetEntry(i);
	if(timestamp_unix>10000) {
	  timestamp_unix_off=timestamp_unix-time_off;
	  // if(timeoff_arg!=0) time_off=timeoff_arg;
	  // else {
	  //   time_off=calc_timeoffset(timestamp_unix); 
	  // }
	  //printf("\n hv: timestamp_unix - timeoff - timestamp_unix_off: %i    %f   %f", timestamp_unix,  time_off, timestamp_unix_off);
	  //printf("\n hv: timestamp_gtu -  timeoff - timestamp_gtu_off: %i    %f   %f", timestamp_gtu,  time_off, timestamp_gtu_off);
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
	  time_o_hv->Fill();
	  timestamp_unix_o_hv->Fill();
	  b_tle_iss_hv->Fill();
	  b_tle_sun_hv->Fill();
	  b_tle_moon_hv->Fill();
	  tick=tick;
	  tick_hv->Fill();
	} else {
	  PosizioneISS.abstime=0;
	  PosizioneISS.timeunix=0;
	  PosizioneISS.lat=0;
	  PosizioneISS.lon=0;
	  PosizioneISS.alt=0;
	  PosizioneISS.x=0;
	  PosizioneISS.y=0;
	  PosizioneISS.z=0;
	  PosizioneISS.vx=0;
	  PosizioneISS.vy=0;
	  PosizioneISS.vz=0;
	  PosizioneISS.CdLon=0;
	  PosizioneISS.CdLat=0;
	  PosizioneISS.X_geo=0;
	  PosizioneISS.Y_geo=0;
	  PosizioneISS.Z_geo=0;
	  b_tle_iss_hv->Fill();
	  PosizioneSun.Sun_from_ISS_isDay=0;
	  PosizioneSun.Sun_from_ISS_alt=0;
	  PosizioneSun.Sun_from_ISS_az=0;
	  PosizioneSun.Sun_from_ISS_dist=0;
	  PosizioneSun.Sun_from_Earth_isDay=0;
	  PosizioneSun.Sun_from_Earth_alt=0;
	  PosizioneSun.Sun_from_Earth_az=0;
	  PosizioneSun.Sun_from_Earth_dist=0;
	  b_tle_sun_hv->Fill();
	  PosizioneMoon.Moon_from_ISS_isDay=0;
	  PosizioneMoon.Moon_from_ISS_alt=0;
	  PosizioneMoon.Moon_from_ISS_az=0;
	  PosizioneMoon.Moon_from_ISS_dist=0;
	  PosizioneMoon.Moon_from_Earth_isDay=0;
	  PosizioneMoon.Moon_from_Earth_alt=0;
	  PosizioneMoon.Moon_from_Earth_az=0;
	  PosizioneMoon.Moon_from_Earth_dist=0;
	  b_tle_moon_hv->Fill();
	  tick=0;
	  tick_hv->Fill();
	  time_o_hv->Fill();
	  timestamp_unix_off=0;
	  timestamp_unix_o_hv->Fill();
	}
      }
    }

    ///////////////////////////////////////////////////////////
    //televent
    ///////////////////////////////////////////////////////////
     if(oldfile->Get("televent")){
      TTree *televent_old = (TTree*)oldfile->Get("televent");
      TTree *televent_new = televent_old->CloneTree();
 
      televent_old->SetBranchStatus("timestamp_unix",1);
      Long64_t nentries_tel = televent_old->GetEntries();
      printf("\n televent: nentries  %i  ", nentries_tel);
      televent_old->SetBranchAddress("timestamp_unix",&timestamp_unix);
      TBranch *time_o_tel = televent_new->Branch("time_offset", &time_off, "time_offset/D");
      TBranch *timestamp_unix_o_tel = televent_new->Branch("timestamp_unix_offset", &timestamp_unix_off, "timestamp_unix_offset/D");
      TBranch *b_tle_iss_tel = televent_new->Branch("positionISS", &PosizioneISS.abstime, "abstime/D:timeunix/D:lat/D:lon/D:alt/D:x/D:y/D:z/D:vx/D:vy/D:vz/D:CdLon/D:CdLat/D:X_geo/D:Y_geo/D:Z_geo/D");
      TBranch *b_tle_sun_tel = televent_new->Branch("positionSun", &PosizioneSun.Sun_from_ISS_isDay, "Sun_from_ISS_isDay/D:Sun_from_ISS_alt/D:Sun_from_ISS_az/D:Sun_from_ISS_dist/D:Sun_from_Earth_isDay/D:Sun_from_Earth_alt/D:Sun_from_Earth_az/D:Sun_from_Earth_dist/D");
      TBranch *b_tle_moon_tel = televent_new->Branch("positionMoon", &PosizioneMoon.Moon_from_ISS_isDay, "Moon_from_ISS_isDay/D:Moon_from_ISS_alt/D:Moon_from_ISS_az/D:Moon_from_ISS_dist/D:Moon_from_Earth_isDay/D:Moon_from_Earth_alt/D:Moon_from_Earth_az/D:Moon_from_Earth_dist/D:FractionMoon/D");
      TBranch *tick_tel = televent_new->Branch("positionInt", &tick, "tick/I");
      
      cont = 0;
      diff1=0;
      diff2=0;

      for (Long64_t i=0;i<nentries_tel; i++) {
	televent_old->GetEntry(i);
	// if(timeoff_arg!=0) time_off=timeoff_arg;
	// else {
	//   time_off=calc_timeoffset(timestamp_unix); 
	// }
	timestamp_unix_off=timestamp_unix-time_off;
	if(timestamp_unix!=0) {
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
	  time_o_tel->Fill();
	  timestamp_unix_o_tel->Fill();
	  b_tle_iss_tel->Fill();
	  b_tle_sun_tel->Fill();
	  b_tle_moon_tel->Fill();
	  tick=tick;
	  tick_tel->Fill();
	} else {
	  PosizioneISS.abstime=0;
	  PosizioneISS.timeunix=0;
	  PosizioneISS.lat=0;
	  PosizioneISS.lon=0;
	  PosizioneISS.alt=0;
	  PosizioneISS.x=0;
	  PosizioneISS.y=0;
	  PosizioneISS.z=0;
	  PosizioneISS.vx=0;
	  PosizioneISS.vy=0;
	  PosizioneISS.vz=0;
	  PosizioneISS.CdLon=0;
	  PosizioneISS.CdLat=0;
	  PosizioneISS.X_geo=0;
	  PosizioneISS.Y_geo=0;
	  PosizioneISS.Z_geo=0;
	  b_tle_iss_tel->Fill();
	  PosizioneSun.Sun_from_ISS_isDay=0;
	  PosizioneSun.Sun_from_ISS_alt=0;
	  PosizioneSun.Sun_from_ISS_az=0;
	  PosizioneSun.Sun_from_ISS_dist=0;
	  PosizioneSun.Sun_from_Earth_isDay=0;
	  PosizioneSun.Sun_from_Earth_alt=0;
	  PosizioneSun.Sun_from_Earth_az=0;
	  PosizioneSun.Sun_from_Earth_dist=0;
	  b_tle_sun_tel->Fill();
	  PosizioneMoon.Moon_from_ISS_isDay=0;
	  PosizioneMoon.Moon_from_ISS_alt=0;
	  PosizioneMoon.Moon_from_ISS_az=0;
	  PosizioneMoon.Moon_from_ISS_dist=0;
	  PosizioneMoon.Moon_from_Earth_isDay=0;
	  PosizioneMoon.Moon_from_Earth_alt=0;
	  PosizioneMoon.Moon_from_Earth_az=0;
	  PosizioneMoon.Moon_from_Earth_dist=0;
	  b_tle_moon_tel->Fill();
	  tick=0;
	  tick_tel->Fill();
	  time_o_tel->Fill();
	  timestamp_unix_off=0;
	  timestamp_unix_o_tel->Fill();
	}
      }
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
    tevent_old->SetBranchStatus("timestamp_unix",1);
    tevent_old->SetBranchStatus("gtu_time",1);
    tevent_old->SetBranchStatus("photon_count_data",1); 
    Long64_t nentries = tevent_old->GetEntries();
    printf("\n D1: nentries %i  ", nentries);
    UChar_t ph[1][1][48][48];
    tevent_old->SetBranchAddress("photon_count_data",&ph);
    TBranch *corr_branch = tevent_new->Branch("photon_count_data", &ph, "photon_count_data[1][1][48][48]/b");

    tevent_old->SetBranchAddress("timestamp_unix",&timestamp_unix);
    tevent_old->SetBranchAddress("gtu_time",&gtu_time);
    TBranch *time_o = tevent_new->Branch("time_offset", &time_off, "time_offset/D");
    TBranch *timestamp_unix_o = tevent_new->Branch("timestamp_unix_offset", &timestamp_unix_off, "timestamp_unix_offset/D");
    TBranch *gtu_time_o = tevent_new->Branch("gtu_time_offset", &gtu_time_off, "gtu_time_offset/D");
    TBranch *b_tle_iss = tevent_new->Branch("positionISS", &PosizioneISS.abstime, "abstime/D:timeunix/D:lat/D:lon/D:alt/D:x/D:y/D:z/D:vx/D:vy/D:vz/D:CdLon/D:CdLat/D:X_geo/D:Y_geo/D:Z_geo/D");
    TBranch *b_tle_sun = tevent_new->Branch("positionSun", &PosizioneSun.Sun_from_ISS_isDay, "Sun_from_ISS_isDay/D:Sun_from_ISS_alt/D:Sun_from_ISS_az/D:Sun_from_ISS_dist/D:Sun_from_Earth_isDay/D:Sun_from_Earth_alt/D:Sun_from_Earth_az/D:Sun_from_Earth_dist/D");
    TBranch *b_tle_moon = tevent_new->Branch("positionMoon", &PosizioneMoon.Moon_from_ISS_isDay, "Moon_from_ISS_isDay/D:Moon_from_ISS_alt/D:Moon_from_ISS_az/D:Moon_from_ISS_dist/D:Moon_from_Earth_isDay/D:Moon_from_Earth_alt/D:Moon_from_Earth_az/D:Moon_from_Earth_dist/D:FractionMoon/D");
    TBranch *tick_0 = tevent_new->Branch("positionInt", &tick, "tick/I");
    
    cont = 0;
    diff1=0;
    diff2=0;
    
   for (Long64_t i=0;i<nentries; i++) {
     tevent_old->GetEntry(i);
     // if(timeoff_arg!=0) time_off=timeoff_arg;
     // else {
     //   time_off=calc_timeoffset(timestamp_unix); 
     // }
     timestamp_unix_off=timestamp_unix-time_off;
     gtu_time_off=gtu_time-time_off;
     //printf("\n tevent: timestamp_unix -  timeoff -  timestamp_unix_off: %i    %f   %f", timestamp_unix,  time_off, timestamp_unix_off);
     //printf("\n tevent: gtu_time -  timeoff -  gtu_time_off: %f    %f   %f", gtu_time,  time_off, gtu_time_off);
     if(timestamp_unix!=0) {
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
       gtu_time_o->Fill();
       b_tle_iss->Fill();
       b_tle_sun->Fill();
       b_tle_moon->Fill();
       tick=tick;
       tick_0->Fill();
     }
     else {
       PosizioneISS.abstime=0;
       PosizioneISS.timeunix=0;
       PosizioneISS.lat=0;
       PosizioneISS.lon=0;
       PosizioneISS.alt=0;
       PosizioneISS.x=0;
       PosizioneISS.y=0;
       PosizioneISS.z=0;
       PosizioneISS.vx=0;
       PosizioneISS.vy=0;
       PosizioneISS.vz=0;
       PosizioneISS.CdLon=0;
       PosizioneISS.CdLat=0;
       PosizioneISS.X_geo=0;
       PosizioneISS.Y_geo=0;
       PosizioneISS.Z_geo=0;
       b_tle_iss->Fill();
       PosizioneSun.Sun_from_ISS_isDay=0;
       PosizioneSun.Sun_from_ISS_alt=0;
       PosizioneSun.Sun_from_ISS_az=0;
       PosizioneSun.Sun_from_ISS_dist=0;
       PosizioneSun.Sun_from_Earth_isDay=0;
       PosizioneSun.Sun_from_Earth_alt=0;
       PosizioneSun.Sun_from_Earth_az=0;
       PosizioneSun.Sun_from_Earth_dist=0;
       b_tle_sun->Fill();
       PosizioneMoon.Moon_from_ISS_isDay=0;
       PosizioneMoon.Moon_from_ISS_alt=0;
       PosizioneMoon.Moon_from_ISS_az=0;
       PosizioneMoon.Moon_from_ISS_dist=0;
       PosizioneMoon.Moon_from_Earth_isDay=0;
       PosizioneMoon.Moon_from_Earth_alt=0;
       PosizioneMoon.Moon_from_Earth_az=0;
       PosizioneMoon.Moon_from_Earth_dist=0;
       b_tle_moon->Fill();
       tick=0;
       tick_0->Fill();
       time_o->Fill();
       timestamp_unix_off=0;
       timestamp_unix_o->Fill();
       gtu_time_off=0;
       gtu_time_o->Fill();
     }
   }
    
    for (Long64_t i=0;i<nentries; i++) {
      tevent_old->GetEntry(i);
      for (Long64_t j=0; j<8; j++) {
	for (Long64_t k=32; k<40; k++) {
	  ph[0][0][j][k]=ph[0][0][j][k]>>1;
	}
      }
      for (Long64_t j=32; j<40; j++) {
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
    tevent_1st_integral_old->SetBranchStatus("timestamp_unix",1);
    tevent_1st_integral_old->SetBranchStatus("gtu_time",1);
    tevent_1st_integral_old->SetBranchStatus("photon_count_data",1); 
    Long64_t nentries_1st = tevent_1st_integral_old->GetEntries();
    printf("\n D2: nentries  %i  ", nentries_1st);
    Float_t ph_1st[1][1][48][48];
    tevent_1st_integral_old->SetBranchAddress("photon_count_data",&ph_1st);
    TBranch *corr_branch_1st = tevent_1st_integral_new->Branch("photon_count_data", &ph_1st, "photon_count_data[1][1][48][48]/F");


    tevent_1st_integral_old->SetBranchAddress("timestamp_unix",&timestamp_unix);
    tevent_1st_integral_old->SetBranchAddress("gtu_time",&gtu_time);
    TBranch *time_o_1st = tevent_1st_integral_new->Branch("time_offset", &time_off, "time_offset/D");
    TBranch *timestamp_unix_o_1st = tevent_1st_integral_new->Branch("timestamp_unix_offset", &timestamp_unix_off, "timestamp_unix_offset/D");
    TBranch *gtu_time_o_1st = tevent_1st_integral_new->Branch("gtu_time_offset", &gtu_time_off, "gtu_time_offset/D");
    TBranch *b_tle_iss_1st = tevent_1st_integral_new->Branch("positionISS", &PosizioneISS.abstime, "abstime/D:timeunix/D:lat/D:lon/D:alt/D:x/D:y/D:z/D:vx/D:vy/D:vz/D:CdLon/D:CdLat/D:X_geo/D:Y_geo/D:Z_geo/D");
    TBranch *b_tle_sun_1st = tevent_1st_integral_new->Branch("positionSun", &PosizioneSun.Sun_from_ISS_isDay, "Sun_from_ISS_isDay/D:Sun_from_ISS_alt/D:Sun_from_ISS_az/D:Sun_from_ISS_dist/D:Sun_from_Earth_isDay/D:Sun_from_Earth_alt/D:Sun_from_Earth_az/D:Sun_from_Earth_dist/D");
    TBranch *b_tle_moon_1st = tevent_1st_integral_new->Branch("positionMoon", &PosizioneMoon.Moon_from_ISS_isDay, "Moon_from_ISS_isDay/D:Moon_from_ISS_alt/D:Moon_from_ISS_az/D:Moon_from_ISS_dist/D:Moon_from_Earth_isDay/D:Moon_from_Earth_alt/D:Moon_from_Earth_az/D:Moon_from_Earth_dist/D:FractionMoon/D");
    TBranch *tick_1st = tevent_1st_integral_new->Branch("positionInt", &tick, "tick/I");

   cont = 0;
   diff1=0;
   diff2=0;

   for (Long64_t i=0;i<nentries_1st; i++) {
     tevent_1st_integral_old->GetEntry(i);
     // if(timeoff_arg!=0) time_off=timeoff_arg;
     // else {
     //   time_off=calc_timeoffset(timestamp_unix); 
     // }
     timestamp_unix_off=timestamp_unix-time_off;
     gtu_time_off=gtu_time-time_off;
     //printf("\n tevent_1st: timestamp_unix - timeoff - timestamp_unix_off: %i    %f   %f", timestamp_unix,  time_off, timestamp_unix_off);
     //printf("\n tevent_1st: timestamp_gtu - timeoff - timestamp_gtu_off: %i    %f   %f", timestamp_gtu,  time_off, timestamp_gtu_off);
     //printf("\n tevent_1st: gtu_time - timeoff - gtu_time_off: %f    %f   %f", gtu_time,  time_off, gtu_time_off);
     if(timestamp_unix!=0) {
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
       time_o_1st->Fill();
       timestamp_unix_o_1st->Fill();
       gtu_time_o_1st->Fill();
       b_tle_iss_1st->Fill();
       b_tle_sun_1st->Fill();
       b_tle_moon_1st->Fill();
       tick=tick;
       tick_1st->Fill();
     } else {
       PosizioneISS.abstime=0;
       PosizioneISS.timeunix=0;
       PosizioneISS.lat=0;
       PosizioneISS.lon=0;
       PosizioneISS.alt=0;
       PosizioneISS.x=0;
       PosizioneISS.y=0;
       PosizioneISS.z=0;
       PosizioneISS.vx=0;
       PosizioneISS.vy=0;
       PosizioneISS.vz=0;
       PosizioneISS.CdLon=0;
       PosizioneISS.CdLat=0;
       PosizioneISS.X_geo=0;
       PosizioneISS.Y_geo=0;
       PosizioneISS.Z_geo=0;
       b_tle_iss_1st->Fill();
       PosizioneSun.Sun_from_ISS_isDay=0;
       PosizioneSun.Sun_from_ISS_alt=0;
       PosizioneSun.Sun_from_ISS_az=0;
       PosizioneSun.Sun_from_ISS_dist=0;
       PosizioneSun.Sun_from_Earth_isDay=0;
       PosizioneSun.Sun_from_Earth_alt=0;
       PosizioneSun.Sun_from_Earth_az=0;
       PosizioneSun.Sun_from_Earth_dist=0;
       b_tle_sun_1st->Fill();
       PosizioneMoon.Moon_from_ISS_isDay=0;
       PosizioneMoon.Moon_from_ISS_alt=0;
       PosizioneMoon.Moon_from_ISS_az=0;
       PosizioneMoon.Moon_from_ISS_dist=0;
       PosizioneMoon.Moon_from_Earth_isDay=0;
       PosizioneMoon.Moon_from_Earth_alt=0;
       PosizioneMoon.Moon_from_Earth_az=0;
       PosizioneMoon.Moon_from_Earth_dist=0;
       b_tle_moon_1st->Fill();
       tick=0;
       tick_1st->Fill();
       time_o_1st->Fill();
       timestamp_unix_off=0;
       timestamp_unix_o_1st->Fill();
       gtu_time_off=0;
       gtu_time_o_1st->Fill();
     }
   }
   
    
    for (Long64_t i=0;i<nentries_1st; i++) {
      tevent_1st_integral_old->GetEntry(i);
      for (Long64_t j=0; j<8; j++) {
	for (Long64_t k=32; k<40; k++) {
	  ph_1st[0][0][j][k]=ph_1st[0][0][j][k]/2;
	}
      }
      for (Long64_t j=32; j<40; j++) {
	for (Long64_t k=8; k<16; k++) {
	  ph_1st[0][0][j][k]=ph_1st[0][0][j][k]/2;
	}
      }
      corr_branch_1st->Fill();
    }
    
    
   
   
    /////////////////////////
    /////////////////////////
    /////////////////////////



    newfile->Write(); 
  }
  
  return 0;
  
}
    
