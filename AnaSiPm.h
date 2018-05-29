//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon May  7 14:31:34 2018 by ROOT version 6.11/02
// from TTree B1/physics
// found on file: PtermcX0_Z175_NOCuD_Fil0_TBR10_ExtSr_60035.root
//////////////////////////////////////////////////////////

#ifndef AnaSiPm_h
#define AnaSiPm_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class AnaSiPm {
public :
	TFile* OutFile;
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        Eabs;
   vector<double>  *EabsComp;
   Double_t        PrePterTrackN;
   vector<double>  *PrePterPart;
   vector<double>  *PrePterEn;
   Double_t        InPterTrackN;
   vector<double>  *InPterPart;
   vector<double>  *InPterEn;
   vector<double>  *InPterEnPrim;
   vector<float>   *InPterTime;
   vector<double>  *InPterX;
   vector<double>  *InPterY;
   vector<double>  *InPterZ;
   vector<double>  *PixelID;
   vector<double>  *PixXPos;
   vector<double>  *PixYPos;
   Double_t        SourceX;
   Double_t        SourceY;
   Double_t        SourceZ;
   vector<double>  *SourceCosX;
   vector<double>  *SourceCosY;
   vector<double>  *SourceCosZ;
   vector<double>  *SourceEne;
   vector<double>  *SourceIsotope;
   Int_t           Nev;

   // List of branches
   TBranch        *b_Eabs;   //!
   TBranch        *b_EabsComp;   //!
   TBranch        *b_PrePterTrackN;   //!
   TBranch        *b_PrePterPart;   //!
   TBranch        *b_PrePterEn;   //!
   TBranch        *b_InPterTrackN;   //!
   TBranch        *b_InPterPart;   //!
   TBranch        *b_InPterEn;   //!
   TBranch        *b_InPterEnPrim;   //!
   TBranch        *b_InPterTime;   //!
   TBranch        *b_InPterX;   //!
   TBranch        *b_InPterY;   //!
   TBranch        *b_InPterZ;   //!
   TBranch        *b_PixelID;   //!
   TBranch        *b_PixXPos;   //!
   TBranch        *b_PixYPos;   //!
   TBranch        *b_SourceX;   //!
   TBranch        *b_SourceY;   //!
   TBranch        *b_SourceZ;   //!
   TBranch        *b_SourceCosX;   //!
   TBranch        *b_SourceCosY;   //!
   TBranch        *b_SourceCosZ;   //!
   TBranch        *b_SourceEne;   //!
   TBranch        *b_SourceIsotope;   //!
   TBranch        *b_Nev;   //!

   AnaSiPm(TString filenameMC, double Soglia=5);
   virtual ~AnaSiPm();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
	double fSoglia;
};

#endif

#ifdef AnaSiPm_cxx
AnaSiPm::AnaSiPm(TString filenameMC, double Soglia) : fChain(0)
{

	TFile *fileMC = new TFile(Form("%s.root",filenameMC.Data()));
	TTree *treeMC = (TTree*)gDirectory->Get("B1");
	
	fSoglia=Soglia;
	
	OutFile = new TFile(Form("%s_Thr%.0lf_Out.root",filenameMC.Data(), fSoglia),"RECREATE");
	
	Init(treeMC);
	
	
}

AnaSiPm::~AnaSiPm()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t AnaSiPm::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t AnaSiPm::LoadTree(Long64_t entry)
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

void AnaSiPm::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   EabsComp = 0;
   PrePterPart = 0;
   PrePterEn = 0;
   InPterPart = 0;
   InPterEn = 0;
   InPterEnPrim = 0;
   InPterTime = 0;
   InPterX = 0;
   InPterY = 0;
   InPterZ = 0;
   PixelID = 0;
   PixXPos = 0;
   PixYPos = 0;
   SourceCosX = 0;
   SourceCosY = 0;
   SourceCosZ = 0;
   SourceEne = 0;
   SourceIsotope = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Eabs", &Eabs, &b_Eabs);
   fChain->SetBranchAddress("EabsComp", &EabsComp, &b_EabsComp);
   fChain->SetBranchAddress("PrePterTrackN", &PrePterTrackN, &b_PrePterTrackN);
   fChain->SetBranchAddress("PrePterPart", &PrePterPart, &b_PrePterPart);
   fChain->SetBranchAddress("PrePterEn", &PrePterEn, &b_PrePterEn);
   fChain->SetBranchAddress("InPterTrackN", &InPterTrackN, &b_InPterTrackN);
   fChain->SetBranchAddress("InPterPart", &InPterPart, &b_InPterPart);
   fChain->SetBranchAddress("InPterEn", &InPterEn, &b_InPterEn);
   fChain->SetBranchAddress("InPterEnPrim", &InPterEnPrim, &b_InPterEnPrim);
   fChain->SetBranchAddress("InPterTime", &InPterTime, &b_InPterTime);
   fChain->SetBranchAddress("InPterX", &InPterX, &b_InPterX);
   fChain->SetBranchAddress("InPterY", &InPterY, &b_InPterY);
   fChain->SetBranchAddress("InPterZ", &InPterZ, &b_InPterZ);
   fChain->SetBranchAddress("PixelID", &PixelID, &b_PixelID);
   fChain->SetBranchAddress("PixXPos", &PixXPos, &b_PixXPos);
   fChain->SetBranchAddress("PixYPos", &PixYPos, &b_PixYPos);
   fChain->SetBranchAddress("SourceX", &SourceX, &b_SourceX);
   fChain->SetBranchAddress("SourceY", &SourceY, &b_SourceY);
   fChain->SetBranchAddress("SourceZ", &SourceZ, &b_SourceZ);
   fChain->SetBranchAddress("SourceCosX", &SourceCosX, &b_SourceCosX);
   fChain->SetBranchAddress("SourceCosY", &SourceCosY, &b_SourceCosY);
   fChain->SetBranchAddress("SourceCosZ", &SourceCosZ, &b_SourceCosZ);
   fChain->SetBranchAddress("SourceEne", &SourceEne, &b_SourceEne);
   fChain->SetBranchAddress("SourceIsotope", &SourceIsotope, &b_SourceIsotope);
   fChain->SetBranchAddress("Nev", &Nev, &b_Nev);
   Notify();
}

Bool_t AnaSiPm::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void AnaSiPm::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t AnaSiPm::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef AnaSiPm_cxx
