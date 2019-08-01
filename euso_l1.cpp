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


Int_t  main(int argc, char *argv[]){ 
  
  gSystem->Load("libTree");

  //Get raw data file from command line
  if ( argc == 1 ) {
    printf( "No arguments were passed. You have to pass filename: file.root\n" ); 
  }
  else {
    TFile *oldfile = new TFile(argv[1]);
    
    //Create a new file
    int len = strlen(argv[1]);
    int corr_len= len-5;
    char outputfile[corr_len];
    strncpy(outputfile, argv[1], corr_len);
    strncat(outputfile,"_l1.root", 8);

    TFile *newfile = new TFile(outputfile,"recreate");
    
    //Clone TTrees that do not need to be corrected in the new file
    TTree *televent_old = (TTree*)oldfile->Get("televent");
    TTree *thk_old = (TTree*)oldfile->Get("thk");
    TTree *ttherm_old = (TTree*)oldfile->Get("ttherm");
    TTree *texp_old = (TTree*)oldfile->Get("texp");
    
    TTree *televent_new = televent_old->CloneTree();
    TTree *thk_new = thk_old->CloneTree();
    TTree *ttherm_new = ttherm_old->CloneTree();
    TTree *texp_new = texp_old->CloneTree();
    
    
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
	for (Long64_t k=32; k<39; k++) {
	  ph[0][0][j][k]=ph[0][0][j][k]>>1;
	}
      }
      for (Long64_t j=8; j<15; j++) {
	for (Long64_t k=8; k<15; k++) {
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
	for (Long64_t k=32; k<39; k++) {
	  ph_1st[0][0][j][k]=ph_1st[0][0][j][k]/2;
	}
      }
      for (Long64_t j=8; j<15; j++) {
	for (Long64_t k=8; k<15; k++) {
	  ph_1st[0][0][j][k]=ph_1st[0][0][j][k]/2;
	}
      }
      corr_branch_1st->Fill();
    }
    
    
    ///////////////////////////////////////////////////////////
    //tevent_2nd_integral
    ///////////////////////////////////////////////////////////
    TTree *tevent_2nd_integral_old = (TTree*)oldfile->Get("tevent_2nd_integral");
    
    tevent_2nd_integral_old->SetBranchStatus("*",1); 
    tevent_2nd_integral_old->SetBranchStatus("photon_count_data",0); 
    TTree *tevent_2nd_integral_new = tevent_2nd_integral_old->CloneTree(0); 
    tevent_2nd_integral_new->CopyEntries(tevent_2nd_integral_old);  // Here we copy the branches */
    
    tevent_2nd_integral_old->SetBranchStatus("*",0); 
    tevent_2nd_integral_old->SetBranchStatus("photon_count_data",1); 
    Long64_t nentries_2nd = tevent_2nd_integral_old->GetEntries();
    Float_t ph_2nd[1][1][48][48];
    tevent_2nd_integral_old->SetBranchAddress("photon_count_data",&ph_2nd);
    TBranch *corr_branch_2nd = tevent_2nd_integral_new->Branch("photon_count_data", &ph_2nd, "photon_count_data[1][1][48][48]/F");
    
    for (Long64_t i=0;i<nentries_2nd; i++) {
      tevent_2nd_integral_old->GetEntry(i);
      for (Long64_t j=40; j<48; j++) {
	for (Long64_t k=32; k<39; k++) {
	  ph_2nd[0][0][j][k]=ph_2nd[0][0][j][k]/2;
	}
      }
      for (Long64_t j=8; j<15; j++) {
	for (Long64_t k=8; k<15; k++) {
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
    
