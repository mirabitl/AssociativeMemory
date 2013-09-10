//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Mar 20 08:33:02 2013 by ROOT version 5.34/05
// from TTree L1TrackTrigger/L1TrackTrigger Analysis info
// found on file: C:/Users/laurent/Documents/SLHC_extr_POS_9.root
//////////////////////////////////////////////////////////

#ifndef L1TrackTrigger_h
#define L1TrackTrigger_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TApplication.h>
// Header file for the classes stored in the TTree if any.
#include <vector>
#include "DCHistogramHandler.h"
#undef USE_HV
#ifndef USE_HV
#define HOUGHLOCAL HoughLocal
#include "HoughLocal.h"
#else
#define HOUGHLOCAL HoughLocal1
#include "HoughLocal1.h"
#endif
using namespace std;
// Fixed size dimensions of array or collections stored in the TTree if any.

class L1TrackTrigger {
 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

  // Declaration of leaf types
  Int_t           evt;
  Int_t           CLUS_n;
  vector<float>   *CLUS_x;
  vector<float>   *CLUS_y;
  vector<float>   *CLUS_z;
  vector<float>   *CLUS_xmc;
  vector<float>   *CLUS_ymc;
  vector<float>   *CLUS_charge;
  vector<int>     *CLUS_layer;
  vector<int>     *CLUS_module;
  vector<int>     *CLUS_ladder;
  vector<int>     *CLUS_seg;
  vector<float>   *CLUS_strip;
  vector<int>     *CLUS_nstrip;
  vector<int>     *CLUS_nsat;
  vector<int>     *CLUS_match;
  vector<vector<int> > *CLUS_tp;
  vector<vector<int> > *CLUS_hits;
  vector<float>   *CLUS_zmc;
  Int_t           STUB_n;
  vector<float>   *STUB_pt;
  vector<float>   *STUB_ptMC;
  vector<float>   *STUB_pxGEN;
  vector<float>   *STUB_pyGEN;
  vector<float>   *STUB_etaGEN;
  vector<int>     *STUB_layer;
  vector<int>     *STUB_module;
  vector<int>     *STUB_ladder;
  vector<int>     *STUB_seg;
  vector<int>     *STUB_strip;
  vector<float>   *STUB_x;
  vector<float>   *STUB_y;
  vector<float>   *STUB_z;
  vector<int>     *STUB_clust1;
  vector<int>     *STUB_clust2;
  vector<float>   *STUB_deltas;
  vector<int>     *STUB_tp;
  vector<int>     *STUB_pdgID;
  vector<float>   *STUB_X0;
  vector<float>   *STUB_Y0;
  vector<float>   *STUB_Z0;
  vector<float>   *STUB_PHI0;

  // List of branches
  TBranch        *b_evt;   //!
  TBranch        *b_CLUS_n;   //!
  TBranch        *b_CLUS_x;   //!
  TBranch        *b_CLUS_y;   //!
  TBranch        *b_CLUS_z;   //!
  TBranch        *b_CLUS_xmc;   //!
  TBranch        *b_CLUS_ymc;   //!
  TBranch        *b_CLUS_charge;   //!
  TBranch        *b_CLUS_layer;   //!
  TBranch        *b_CLUS_module;   //!
  TBranch        *b_CLUS_ladder;   //!
  TBranch        *b_CLUS_seg;   //!
  TBranch        *b_CLUS_strip;   //!
  TBranch        *b_CLUS_nstrip;   //!
  TBranch        *b_CLUS_nsat;   //!
  TBranch        *b_CLUS_match;   //!
  TBranch        *b_CLUS_tp;   //!
  TBranch        *b_CLUS_hits;   //!
  TBranch        *b_CLUS_zmc;   //!
  TBranch        *b_STUB_n;   //!
  TBranch        *b_STUB_pt;   //!
  TBranch        *b_STUB_ptMC;   //!
  TBranch        *b_STUB_pxGEN;   //!
  TBranch        *b_STUB_pyGEN;   //!
  TBranch        *b_STUB_etaGEN;   //!
  TBranch        *b_STUB_layer;   //!
  TBranch        *b_STUB_module;   //!
  TBranch        *b_STUB_ladder;   //!
  TBranch        *b_STUB_seg;   //!
  TBranch        *b_STUB_strip;   //!
  TBranch        *b_STUB_x;   //!
  TBranch        *b_STUB_y;   //!
  TBranch        *b_STUB_z;   //!
  TBranch        *b_STUB_clust1;   //!
  TBranch        *b_STUB_clust2;   //!
  TBranch        *b_STUB_deltas;   //!
  TBranch        *b_STUB_tp;   //!
  TBranch        *b_STUB_pdgID;   //!
  TBranch        *b_STUB_X0;   //!
  TBranch        *b_STUB_Y0;   //!
  TBranch        *b_STUB_Z0;   //!
  TBranch        *b_STUB_PHI0;   //!

  L1TrackTrigger(TTree *tree);
  L1TrackTrigger(){;}
	 
  void do_ana(std::string filename,uint32_t nevmax=0);
  void event_hough(DCHistogramHandler* d);
  void fill_histos(DCHistogramHandler* d);
  void associate();
  virtual ~L1TrackTrigger();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop(TApplication* theApp=NULL);
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
	 
  DCHistogramHandler rootHandler_;
  HOUGHLOCAL* htl_;
  HOUGHLOCAL* htp_;
  uint32_t ngoodmc_,nmiss_,nfake_;
};

#endif

#ifdef L1TrackTrigger_cxx
L1TrackTrigger::L1TrackTrigger(TTree *tree) : fChain(0) 
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("C:/Users/laurent/Documents/SLHC_extr_POS_9.root");
    if (!f || !f->IsOpen()) {
      f = new TFile("C:/Users/laurent/Documents/SLHC_extr_POS_9.root");
    }
    f->GetObject("L1TrackTrigger",tree);

  }
  Init(tree);
}

L1TrackTrigger::~L1TrackTrigger()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t L1TrackTrigger::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t L1TrackTrigger::LoadTree(Long64_t entry)
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

void L1TrackTrigger::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set object pointer
  CLUS_x = 0;
  CLUS_y = 0;
  CLUS_z = 0;
  CLUS_xmc = 0;
  CLUS_ymc = 0;
  CLUS_charge = 0;
  CLUS_layer = 0;
  CLUS_module = 0;
  CLUS_ladder = 0;
  CLUS_seg = 0;
  CLUS_strip = 0;
  CLUS_nstrip = 0;
  CLUS_nsat = 0;
  CLUS_match = 0;
  CLUS_tp = 0;
  CLUS_hits = 0;
  CLUS_zmc = 0;
  STUB_pt = 0;
  STUB_ptMC = 0;
  STUB_pxGEN = 0;
  STUB_pyGEN = 0;
  STUB_etaGEN = 0;
  STUB_layer = 0;
  STUB_module = 0;
  STUB_ladder = 0;
  STUB_seg = 0;
  STUB_strip = 0;
  STUB_x = 0;
  STUB_y = 0;
  STUB_z = 0;
  STUB_clust1 = 0;
  STUB_clust2 = 0;
  STUB_deltas = 0;
  STUB_tp = 0;
  STUB_pdgID = 0;
  STUB_X0 = 0;
  STUB_Y0 = 0;
  STUB_Z0 = 0;
  STUB_PHI0 = 0;
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("evt", &evt, &b_evt);
  fChain->SetBranchAddress("CLUS_n", &CLUS_n, &b_CLUS_n);
  fChain->SetBranchAddress("CLUS_x", &CLUS_x, &b_CLUS_x);
  fChain->SetBranchAddress("CLUS_y", &CLUS_y, &b_CLUS_y);
  fChain->SetBranchAddress("CLUS_z", &CLUS_z, &b_CLUS_z);
  fChain->SetBranchAddress("CLUS_xmc", &CLUS_xmc, &b_CLUS_xmc);
  fChain->SetBranchAddress("CLUS_ymc", &CLUS_ymc, &b_CLUS_ymc);
  fChain->SetBranchAddress("CLUS_charge", &CLUS_charge, &b_CLUS_charge);
  fChain->SetBranchAddress("CLUS_layer", &CLUS_layer, &b_CLUS_layer);
  fChain->SetBranchAddress("CLUS_module", &CLUS_module, &b_CLUS_module);
  fChain->SetBranchAddress("CLUS_ladder", &CLUS_ladder, &b_CLUS_ladder);
  fChain->SetBranchAddress("CLUS_seg", &CLUS_seg, &b_CLUS_seg);
  fChain->SetBranchAddress("CLUS_strip", &CLUS_strip, &b_CLUS_strip);
  fChain->SetBranchAddress("CLUS_nstrip", &CLUS_nstrip, &b_CLUS_nstrip);
  fChain->SetBranchAddress("CLUS_nsat", &CLUS_nsat, &b_CLUS_nsat);
  fChain->SetBranchAddress("CLUS_match", &CLUS_match, &b_CLUS_match);
  fChain->SetBranchAddress("CLUS_tp", &CLUS_tp, &b_CLUS_tp);
  fChain->SetBranchAddress("CLUS_hits", &CLUS_hits, &b_CLUS_hits);
  fChain->SetBranchAddress("CLUS_zmc", &CLUS_zmc, &b_CLUS_zmc);
  fChain->SetBranchAddress("STUB_n", &STUB_n, &b_STUB_n);
  fChain->SetBranchAddress("STUB_pt", &STUB_pt, &b_STUB_pt);
  fChain->SetBranchAddress("STUB_ptMC", &STUB_ptMC, &b_STUB_ptMC);
  fChain->SetBranchAddress("STUB_pxGEN", &STUB_pxGEN, &b_STUB_pxGEN);
  fChain->SetBranchAddress("STUB_pyGEN", &STUB_pyGEN, &b_STUB_pyGEN);
  fChain->SetBranchAddress("STUB_etaGEN", &STUB_etaGEN, &b_STUB_etaGEN);
  fChain->SetBranchAddress("STUB_layer", &STUB_layer, &b_STUB_layer);
  fChain->SetBranchAddress("STUB_module", &STUB_module, &b_STUB_module);
  fChain->SetBranchAddress("STUB_ladder", &STUB_ladder, &b_STUB_ladder);
  fChain->SetBranchAddress("STUB_seg", &STUB_seg, &b_STUB_seg);
  fChain->SetBranchAddress("STUB_strip", &STUB_strip, &b_STUB_strip);
  fChain->SetBranchAddress("STUB_x", &STUB_x, &b_STUB_x);
  fChain->SetBranchAddress("STUB_y", &STUB_y, &b_STUB_y);
  fChain->SetBranchAddress("STUB_z", &STUB_z, &b_STUB_z);
  fChain->SetBranchAddress("STUB_clust1", &STUB_clust1, &b_STUB_clust1);
  fChain->SetBranchAddress("STUB_clust2", &STUB_clust2, &b_STUB_clust2);
  fChain->SetBranchAddress("STUB_deltas", &STUB_deltas, &b_STUB_deltas);
  fChain->SetBranchAddress("STUB_tp", &STUB_tp, &b_STUB_tp);
  fChain->SetBranchAddress("STUB_pdgID", &STUB_pdgID, &b_STUB_pdgID);
  fChain->SetBranchAddress("STUB_X0", &STUB_X0, &b_STUB_X0);
  fChain->SetBranchAddress("STUB_Y0", &STUB_Y0, &b_STUB_Y0);
  fChain->SetBranchAddress("STUB_Z0", &STUB_Z0, &b_STUB_Z0);
  fChain->SetBranchAddress("STUB_PHI0", &STUB_PHI0, &b_STUB_PHI0);
  Notify();
}

Bool_t L1TrackTrigger::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void L1TrackTrigger::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t L1TrackTrigger::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
#endif // #ifdef L1TrackTrigger_cxx
