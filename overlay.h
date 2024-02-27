//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jul 13 17:55:04 2023 by ROOT version 6.28/02
// from TTree ECalBarrelCollectionTuple/columnwise ntuple with LCIO data
// found on file: /home/sestini/photons_Crilin_v2/simulation/ddsim/crilin_v2/175GeV_1.57rad/photonGun_175GeV_1.57rad_gen_crilin_v2.root
//////////////////////////////////////////////////////////

#ifndef overlay_h
#define overlay_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <string>

// Header file for the classes stored in the TTree if any.

class overlay {
public :

   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         primary_energy;
   Int_t           evevt;
   Int_t           evrun;
   Float_t         evwgt;
   Long64_t        evtim;
   Float_t         evsig;
   Float_t         evene;
   Float_t         evpoe;
   Float_t         evpop;
   Int_t           evnch;
   Char_t          evpro[1];   //[evnch]
   Int_t           nsch;
   Int_t           scori[25613];   //[nsch]
   Int_t           scci0[25613];   //[nsch]
   Int_t           scci1[25613];   //[nsch]
   Float_t         scpox[25613];   //[nsch]
   Float_t         scpoy[25613];   //[nsch]
   Float_t         scpoz[25613];   //[nsch]
   Float_t         scene[25613];   //[nsch]
   Int_t           scmcc[25613];   //[nsch]
   Float_t         sctim[25613][50];   //[nsch]

   // List of branches
   TBranch        *b_primary_energy;   //!
   TBranch        *b_evevt;   //!
   TBranch        *b_evrun;   //!
   TBranch        *b_evwgt;   //!
   TBranch        *b_evtim;   //!
   TBranch        *b_evsig;   //!
   TBranch        *b_evene;   //!
   TBranch        *b_evpoe;   //!
   TBranch        *b_evpop;   //!
   TBranch        *b_evnch;   //!
   TBranch        *b_evpro;   //!
   TBranch        *b_nsch;   //!
   TBranch        *b_scori;   //!
   TBranch        *b_scci0;   //!
   TBranch        *b_scci1;   //!
   TBranch        *b_scpox;   //!
   TBranch        *b_scpoy;   //!
   TBranch        *b_scpoz;   //!
   TBranch        *b_scene;   //!
   TBranch        *b_scmcc;   //!
   TBranch        *b_sctim;   //!

   overlay(TTree *tree=0);
   virtual ~overlay();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef overlay_cxx
overlay::overlay(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      std::string fname = "/lustre/cmsdata/fnardi/crilin_photon_data/uniform/photonGun_10_175GeV_1.57rad_gen_crilin_v2.root";
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(fname.c_str());
      if (!f || !f->IsOpen()) {
         f = new TFile(fname.c_str());
      }
      f->GetObject("ECalBarrelCollectionTuple",tree);

   }
   Init(tree);
}

overlay::~overlay()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t overlay::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t overlay::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void overlay::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);
   fChain->SetBranchAddress("primary_energy", &primary_energy, &b_primary_energy);
   fChain->SetBranchAddress("evevt", &evevt, &b_evevt);
   fChain->SetBranchAddress("evrun", &evrun, &b_evrun);
   fChain->SetBranchAddress("evwgt", &evwgt, &b_evwgt);
   fChain->SetBranchAddress("evtim", &evtim, &b_evtim);
   fChain->SetBranchAddress("evsig", &evsig, &b_evsig);
   fChain->SetBranchAddress("evene", &evene, &b_evene);
   fChain->SetBranchAddress("evpoe", &evpoe, &b_evpoe);
   fChain->SetBranchAddress("evpop", &evpop, &b_evpop);
   fChain->SetBranchAddress("evnch", &evnch, &b_evnch);
   fChain->SetBranchAddress("evpro", &evpro, &b_evpro);
   fChain->SetBranchAddress("nsch", &nsch, &b_nsch);
   fChain->SetBranchAddress("scori", scori, &b_scori);
   fChain->SetBranchAddress("scci0", scci0, &b_scci0);
   fChain->SetBranchAddress("scci1", scci1, &b_scci1);
   fChain->SetBranchAddress("scpox", scpox, &b_scpox);
   fChain->SetBranchAddress("scpoy", scpoy, &b_scpoy);
   fChain->SetBranchAddress("scpoz", scpoz, &b_scpoz);
   fChain->SetBranchAddress("scene", scene, &b_scene);
   fChain->SetBranchAddress("scmcc", scmcc, &b_scmcc);
   fChain->SetBranchAddress("sctim", sctim, &b_sctim);
   Notify();
}

Bool_t overlay::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void overlay::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t overlay::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef overlay_cxx
