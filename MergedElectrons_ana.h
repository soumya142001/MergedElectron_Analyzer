//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Sep 14 15:26:55 2023 by ROOT version 6.24/06
// from TTree tree/tree
// found on file: /home/soumya/CMSSW_files/NTuples/root_files/ZtoeepT40to300_ntuple.root
//////////////////////////////////////////////////////////

#ifndef MergedElectrons_ana_h
#define MergedElectrons_ana_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
#include <vector>
#include <fstream>
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <iostream>
#include "TString.h"
#include <bitset>



class MergedElectrons_ana : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTreeReader     fReader_MC;
   TTreeReader     fReader_Data;
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

  //Adding a new branch to the tree
  // TBranch *b1 = fChain->GetTree()->Branch("dEraw_over_Ecorr","std::vector<float>",&dEraw_over_Ecorr_);
  
  
  // Readers to access the data (delete the ones you do not need).

  TTreeReaderValue<UInt_t> run = {fReader, "run"};
  TTreeReaderValue<UInt_t> event = {fReader, "event"};
  TTreeReaderValue<UInt_t> lumi = {fReader, "lumi"};

  //Gsf electron or reco::Electron variables
  TTreeReaderValue<UInt_t> ngsf = {fReader, "ngsf"};
  TTreeReaderArray<float> gsfPx = {fReader, "gsfPx"};
  TTreeReaderArray<float> gsfPy = {fReader, "gsfPy"};
  TTreeReaderArray<float> gsfPz = {fReader, "gsfPz"};
  TTreeReaderArray<float> gsfE = {fReader, "gsfE"};
  TTreeReaderArray<float> gsfPt = {fReader, "gsfPt"};
  TTreeReaderArray<float> gsfEta = {fReader, "gsfEta"};
  TTreeReaderArray<float> gsfPhi = {fReader, "gsfPhi"};
  TTreeReaderArray<int> gsfId = {fReader, "gsfId"};
  TTreeReaderArray<float> gsfscPixCharge = {fReader, "gsfscPixCharge"};
  TTreeReaderArray<int> gsfseedtrackcharge = {fReader, "gsfseedtrackcharge"};
  TTreeReaderArray<float> gsfseedtrackPtmode = {fReader, "gsfseedtrackPtmode"};
  TTreeReaderArray<float> gsfseedtrackEtamode = {fReader, "gsfseedtrackEtamode"};
  TTreeReaderArray<float> gsfseedtrackPhimode = {fReader, "gsfseedtrackPhimode"};
  TTreeReaderArray<float> gsfseedtrackPmode = {fReader, "gsfseedtrackPmode"};
  TTreeReaderArray<float> gsfseedtrackPxmode = {fReader, "gsfseedtrackPxmode"};
  TTreeReaderArray<float> gsfseedtrackPymode = {fReader, "gsfseedtrackPymode"};
  TTreeReaderArray<float> gsfseedtrackPzmode = {fReader, "gsfseedtrackPzmode"};
  TTreeReaderArray<float> gsfseedtrackPt = {fReader, "gsfseedtrackPt"};
  TTreeReaderArray<float> gsfseedtrackEta = {fReader, "gsfseedtrackEta"};
  TTreeReaderArray<float> gsfseedtrackPhi = {fReader, "gsfseedtrackPhi"};
  TTreeReaderArray<float> gsfseedtrackP = {fReader, "gsfseedtrackP"};
  TTreeReaderArray<float> gsfseedtrackPx = {fReader, "gsfseedtrackPx"};
  TTreeReaderArray<float> gsfseedtrackPy = {fReader, "gsfseedtrackPy"};
  TTreeReaderArray<float> gsfseedtrackPz = {fReader, "gsfseedtrackPz"};
  TTreeReaderArray<float> gsfsigmaEtaEta = {fReader, "gsfsigmaEtaEta"};
  TTreeReaderArray<float> gsfsigmaIetaIeta = {fReader, "gsfsigmaIetaIeta"};
  TTreeReaderArray<float> gsfsigmaphiIphi = {fReader, "gsfsigmaphiIphi"};
  TTreeReaderArray<float> gsfe1x5 = {fReader, "gsfe1x5"};
  TTreeReaderArray<float> gsfe5x5 = {fReader, "gsfe5x5"};
  TTreeReaderArray<float> gsfe2x5Max = {fReader, "gsfe2x5Max"};
  TTreeReaderArray<float> gsfr9 = {fReader, "gsfr9"};
  TTreeReaderArray<float> gsfeSuperClusterOverP = {fReader, "gsfeSuperClusterOverP"};
  TTreeReaderArray<float> gsfeSeedClusterOverP = {fReader, "gsfeSeedClusterOverP"};
  TTreeReaderArray<float> gsfeSeedClusterOverPout = {fReader, "gsfeSeedClusterOverPout"};
  TTreeReaderArray<float> gsfeEleClusterOverPout = {fReader, "gsfeEleClusterOverPout"};
  TTreeReaderArray<float> gsfdeltaEtaSuperClusterTrackAtVtx = {fReader, "gsfdeltaEtaSuperClusterTrackAtVtx"};
  TTreeReaderArray<float> gsfdeltaEtaSeedClusterTrackAtCalo = {fReader, "gsfdeltaEtaSeedClusterTrackAtCalo"};
  TTreeReaderArray<float> gsfdeltaEtaEleClusterTrackAtCalo = {fReader, "gsfdeltaEtaEleClusterTrackAtCalo"};
  TTreeReaderArray<float> gsfdeltaPhiEleClusterTrackAtCalo = {fReader, "gsfdeltaPhiEleClusterTrackAtCalo"};
  TTreeReaderArray<float> gsfdeltaPhiSuperClusterTrackAtVtx = {fReader, "gsfdeltaPhiSuperClusterTrackAtVtx"};
  TTreeReaderArray<float> gsfdeltaPhiSeedClusterTrackAtCalo = {fReader, "gsfdeltaPhiSeedClusterTrackAtCalo"};
  TTreeReaderArray<float> gsfdeltaEtaSeedClusterTrackAtVtx = {fReader, "gsfdeltaEtaSeedClusterTrackAtVtx"};
  TTreeReaderArray<float> gsfrawEnergy = {fReader, "gsfrawEnergy"};
  TTreeReaderArray<float> gsfcorrectedEcalEnergy = {fReader, "gsfcorrectedEcalEnergy"};
  TTreeReaderArray<int> gsfseedieta = {fReader, "gsfseedieta"};
  TTreeReaderArray<int> gsfseediphi = {fReader, "gsfseediphi"};
  TTreeReaderArray<float> gsfseedclusterenergy = {fReader, "gsfseedclusterenergy"};
  TTreeReaderArray<float> gsfseedclustereta = {fReader, "gsfseedclustereta"};
  TTreeReaderArray<float> gsfseedclusterphi = {fReader, "gsfseedclusterphi"};
  TTreeReaderValue<vector<bool>> gsfisEB = {fReader, "gsfisEB"};
  TTreeReaderValue<vector<bool>> gsfisEE = {fReader, "gsfisEE"};
  TTreeReaderArray<float> gsfdr03TkSumPt = {fReader, "gsfdr03TkSumPt"};
  TTreeReaderArray<float> gsfdr03TkSumPtHEEP = {fReader, "gsfdr03TkSumPtHEEP"};
  TTreeReaderArray<float> gsfdr03EcalRecHitSumEt = {fReader, "gsfdr03EcalRecHitSumEt"};
  TTreeReaderArray<float> gsfdr04TkSumPt = {fReader, "gsfdr04TkSumPt"};
  TTreeReaderArray<float> gsfdr04TkSumPtHEEP = {fReader, "gsfdr04TkSumPtHEEP"};
  TTreeReaderArray<float> gsfdr04EcalRecHitSumEt = {fReader, "gsfdr04EcalRecHitSumEt"};
  TTreeReaderArray<float> gsfdr04HcalTowerSumEt = {fReader, "gsfdr04HcalTowerSumEt"};
  TTreeReaderArray<float> gsfdr04HcalTowerSumEtBc = {fReader, "gsfdr04HcalTowerSumEtBc"};
  TTreeReaderArray<float> gsfhcalOverEcal = {fReader, "gsfhcalOverEcal"};
  TTreeReaderArray<float> gsfhcalOverEcalBc = {fReader, "gsfhcalOverEcalBc"};
  TTreeReaderArray<float> gsffull5x5_sigmaEtaEta = {fReader, "gsffull5x5_sigmaEtaEta"};
  TTreeReaderArray<float> gsffull5x5_sigmaIetaIeta = {fReader, "gsffull5x5_sigmaIetaIeta"};
  TTreeReaderArray<float> gsffull5x5_sigmaIphiIphi = {fReader, "gsffull5x5_sigmaIphiIphi"};
  TTreeReaderArray<float> gsffull5x5_e1x5 = {fReader, "gsffull5x5_e1x5"};
  TTreeReaderArray<float> gsffull5x5_e2x5Max = {fReader, "gsffull5x5_e2x5Max"};
  TTreeReaderArray<float> gsffull5x5_e5x5 = {fReader, "gsffull5x5_e5x5"};
  TTreeReaderArray<float> gsffull5x5_r9 = {fReader, "gsffull5x5_r9"};
  TTreeReaderArray<float> gsffull5x5_hcalOverEcal = {fReader, "gsffull5x5_hcalOverEcal"};
  TTreeReaderArray<float> gsffull5x5_hcalOverEcalBc = {fReader, "gsffull5x5_hcalOverEcalBc"};
  TTreeReaderValue<vector<bool>> gsfambiguous = {fReader, "gsfambiguous"};

  //General Tracks or KF track variables
  TTreeReaderValue<UInt_t> ntrack = {fReader, "ntrack"};
  TTreeReaderArray<float> trackchi2 = {fReader, "trackchi2"};
  TTreeReaderArray<float> trackndof = {fReader, "trackndof"};
  TTreeReaderArray<float> tracknormalizedChi2 = {fReader, "tracknormalizedChi2"};
  TTreeReaderArray<int> trackcharge = {fReader, "trackcharge"};
  TTreeReaderArray<float> trackqoverp = {fReader, "trackqoverp"};
  TTreeReaderArray<float> tracktheta = {fReader, "tracktheta"};
  TTreeReaderArray<float> tracklambda = {fReader, "tracklambda"};
  TTreeReaderArray<float> trackdxy = {fReader, "trackdxy"};
  TTreeReaderArray<float> trackdz = {fReader, "trackdz"};
  TTreeReaderArray<float> trackp2 = {fReader, "trackp2"};
  TTreeReaderArray<float> trackp = {fReader, "trackp"};
  TTreeReaderArray<float> trackpt2 = {fReader, "trackpt2"};
  TTreeReaderArray<float> trackpt = {fReader, "trackpt"};
  TTreeReaderArray<float> trackpx = {fReader, "trackpx"};
  TTreeReaderArray<float> trackpy = {fReader, "trackpy"};
  TTreeReaderArray<float> trackpz = {fReader, "trackpz"};
  TTreeReaderArray<float> trackphi = {fReader, "trackphi"};
  TTreeReaderArray<float> tracketa = {fReader, "tracketa"};
  TTreeReaderArray<float> trackvx = {fReader, "trackvx"};
  TTreeReaderArray<float> trackvy = {fReader, "trackvy"};
  TTreeReaderArray<float> trackvz = {fReader, "trackvz"};
  TTreeReaderArray<float> trackbeta = {fReader, "trackbeta"};
  TTreeReaderArray<float> trackpError = {fReader, "trackpError"};
  TTreeReaderArray<float> trackptError2 = {fReader, "trackptError2"};
  TTreeReaderArray<float> trackptError = {fReader, "trackptError"};
  TTreeReaderArray<float> trackthetaError = {fReader, "trackthetaError"};
  TTreeReaderArray<float> tracklambdaError = {fReader, "tracklambdaError"};
  TTreeReaderArray<float> tracketaError = {fReader, "tracketaError"};
  TTreeReaderArray<float> trackphiError = {fReader, "trackphiError"};
  TTreeReaderArray<float> trackdxyError = {fReader, "trackdxyError"};
  TTreeReaderArray<float> trackdzError = {fReader, "trackdzError"};
  TTreeReaderArray<float> trackbetaError = {fReader, "trackbetaError"};
  TTreeReaderArray<int> tracknumberofValidHits = {fReader, "tracknumberofValidHits"};
  TTreeReaderArray<int> tracknumberOfLostHits = {fReader, "tracknumberOfLostHits"};
  TTreeReaderArray<int> trackmissingInnerHits = {fReader, "trackmissingInnerHits"};
  TTreeReaderArray<int> trackmissingOuterHits = {fReader, "trackmissingOuterHits"};
  TTreeReaderArray<float> trackvalidFraction = {fReader, "trackvalidFraction"};
  TTreeReaderArray<int> trackrecHitSize = {fReader, "trackrecHitSize"};

  //GsfTrack variables
  TTreeReaderValue<UInt_t> ngsftrack = {fReader, "ngsftrack"};
  TTreeReaderArray<float> gsftrackchargeMode = {fReader, "gsftrackchargeMode"};
  TTreeReaderArray<float> gsftrackqoverpMode = {fReader, "gsftrackqoverpMode"};
  TTreeReaderArray<float> gsftrackthetaMode = {fReader, "gsftrackthetaMode"};
  TTreeReaderArray<float> gsftracklambdaMode = {fReader, "gsftracklambdaMode"};
  TTreeReaderArray<float> gsftrackpMode = {fReader, "gsftrackpMode"};
  TTreeReaderArray<float> gsftrackptMode = {fReader, "gsftrackptMode"};
  TTreeReaderArray<float> gsftrackpxMode = {fReader, "gsftrackpxMode"};
  TTreeReaderArray<float> gsftrackpyMode = {fReader, "gsftrackpyMode"};
  TTreeReaderArray<float> gsftrackpzMode = {fReader, "gsftrackpzMode"};
  TTreeReaderArray<float> gsftrackphiMode = {fReader, "gsftrackphiMode"};
  TTreeReaderArray<float> gsftracketaMode = {fReader, "gsftracketaMode"};
  TTreeReaderArray<float> gsftrackpModeError = {fReader, "gsftrackpModeError"};
  TTreeReaderArray<float> gsftrackptModeError = {fReader, "gsftrackptModeError"};
  TTreeReaderArray<float> gsftrackthetaModeError = {fReader, "gsftrackthetaModeError"};
  TTreeReaderArray<float> gsftracklambdaModeError = {fReader, "gsftracklambdaModeError"};
  TTreeReaderArray<float> gsftracketaModeError = {fReader, "gsftracketaModeError"};
  TTreeReaderArray<float> gsftrackphiModeError = {fReader, "gsftrackphiModeError"};
  TTreeReaderArray<float> gsftrackPt = {fReader, "gsftrackPt"};
  TTreeReaderArray<float> gsftrackEta = {fReader, "gsftrackEta"};
  TTreeReaderArray<float> gsftrackPhi = {fReader, "gsftrackPhi"};
  TTreeReaderArray<float> gsftrackP = {fReader, "gsftrackP"};
  TTreeReaderArray<float> gsftrackPx = {fReader, "gsftrackPx"};
  TTreeReaderArray<float> gsftrackPy = {fReader, "gsftrackPy"};
  TTreeReaderArray<float> gsftrackPz = {fReader, "gsftrackPz"};
  TTreeReaderArray<float> gsftrackdxy = {fReader, "gsftrackdxy"};
  TTreeReaderArray<float> gsftrackdz = {fReader, "gsftrackdz"};
  TTreeReaderArray<int> gsftrackcharge = {fReader, "gsftrackcharge"};
  TTreeReaderArray<float> gsftrackchi2 = {fReader, "gsftrackchi2"};
  TTreeReaderArray<int> gsftrackndof = {fReader, "gsftrackndof"};
  TTreeReaderArray<float> gsftracknormalizedChi2 = {fReader, "gsftracknormalizedChi2"};
  
  //GenParticle variables
  TTreeReaderValue<UInt_t> ngenparticle = {fReader_MC, "ngenparticle"};
  TTreeReaderArray<int> nmother = {fReader_MC, "nmother"};
  TTreeReaderArray<int> ndaughter = {fReader_MC, "ndaughter"};
  TTreeReaderArray<int> genstatus = {fReader_MC, "genstatus"};
  TTreeReaderArray<int> pdgid = {fReader_MC, "pdgid"};
  TTreeReaderArray<int> momid = {fReader_MC, "momid"};
  TTreeReaderArray<vector<int>> dauid = {fReader_MC, "dauid"};
  TTreeReaderArray<float> genpT = {fReader_MC, "genpT"};
  TTreeReaderArray<float> genEta = {fReader_MC, "genEta"};
  TTreeReaderArray<float> genPhi = {fReader_MC, "genPhi"};
  TTreeReaderArray<float> genM = {fReader_MC, "genM"};
  TTreeReaderArray<float> genE = {fReader_MC, "genE"};
  TTreeReaderArray<float> genvx = {fReader_MC, "genvx"};
  TTreeReaderArray<float> genvy = {fReader_MC, "genvy"};
  TTreeReaderArray<float> genvz = {fReader_MC, "genvz"};
  
   MergedElectrons_ana(TTree * /*tree*/ =0) { }
   virtual ~MergedElectrons_ana() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();
 
  //User defined functions here

  void SetHstFileName(const char *HstFileName){_HstFileName = HstFileName;}
  void SetSumFileName(const char *SumFileName){_SumFileName = SumFileName;}
  void SetTreeFileName(const char *TreeFileName){TreeFileName_ = TreeFileName;}
  void SetVerbose(int verbose){_verbosity = verbose;}
  void BookHistograms();
  void SetData(int data){_data=data;}
  void sort(int opt);
public:
  struct Hists{
    //Define your histograms here
    TH1F *nGsf;
    TH1F *pT;
    TH1F *dntrack_ngsftrack;
    TH1F *dngsf_ngsftrack;
    TH1F *ntrack_dR01;
    TH1F *ntrack_dR0p05;
    TH1F *Gsfprop_histo[50];
    TH1F *trackprop_histo[20];
    TH1F *misc_histo[100];
    TH1F *genElec_histo[100];
    TH1F *gsftrack_histo[100];
    TH1F *gsf_histo[100];
    TH1F *genPart_histo[100];
    TH1F *genRHN_histo[100];
    TH1F *mergedgsf_histo[100];
    //2D histograms
    TH2F *gen_histo[100];
  };

  struct Lepton{
    int id;
    TLorentzVector v;
    int momtag; //A tag use to identify mother (1-RHN, 2-W, 3-No mom)
    float seedtrkpx,seedtrkpy,seedtrkpz;
    float seedtrkpT,seedtrkEta,seedtrkPhi;
    float px,py,pz;
    float pT,eta,phi;
    float vx,vy,vz;
    float rawEnergy;
    float correctedEnergy;
  };
    

protected:
  Hists h;

private:
  TFile *_HstFile;
  const char *_HstFileName;
  const char *_SumFileName;
  const char *_MyFileName;
  int _verbosity;
  int _data;
  int nEvtTotal, nEvtpass;
  int nrecopass,ntrackpass,nonereco;
  int nrecodRpass,ntrackdRpass;
  int nrecodRpass2,ntrackdRpass2;
  int ngendRl0p1;
  int ngen_dRl0p1, ngen_dRg0p2;
  int n0reco_dRl0p1,n1reco_dRl0p1,n2reco_dRl0p1;
  int n0reco_dRg0p2,n1reco_dRg0p2,n2reco_dRg0p2;

  /*********************************

    Declare Vectors and arrays here
   
  **********************************/

  vector<Lepton> gsftrack;
  vector<Lepton> track;
  vector<Lepton> genElec;
  vector<Lepton> gsfElec;
  vector<Lepton> gsfElec_momRHN;
  vector<Lepton> gsfElec_momW;
  vector<Lepton> genElec_momRHN;
  vector<Lepton> genElec_momW;
  vector<Lepton> gsfseedtrack;
  vector<Lepton> genJPsi;
  vector<Lepton> gsftrack_dR0p1;
  vector<Lepton> gsftrack_dR0p05;
  vector<Lepton> genPart;
  vector<Lepton> genRHN;
  

  //For Tree making
  const char *TreeFileName_;
  TFile      *treefile_;
  TTree      *mytree;

  //Variables that goes into myTreeroot files

  //Integers
  unsigned int ngsf_, ngsf_momRHN_, ngsftrack_, ngenelectrons_, ngenJPsi_;

  //Event number filled for each electron if 2 electron then 2 event numbers (to identify them when dataframe is exploded in python NN)

  vector<int>  gsfeventno_;
  
  //Direct and derived Gsf Variables

  vector<float> gsfPt_;
  vector<float> gsfEta_;
  vector<float> gsfPhi_;
  vector<float> gsfseedtrackPt_;
  vector<float> gsfseedtrackEta_;
  vector<float> gsfseedtrackPhi_;
  vector<float> gsfseedtrackP_;
  vector<float> gsfsigmaEtaEta_;
  vector<float> gsfsigmaIetaIeta_;
  vector<float> gsfsigmaIphiIphi_;
  vector<float> gsfr9_;
  vector<float> gsfeSuperClusterOverP_;
  vector<float> gsfeSeedClusterOverP_;
  vector<float> gsfdeltaEtaSuperClusterTrackAtVtx_;
  vector<float> gsfdeltaEtaSeedClusterTrackAtCalo_;
  vector<float> gsfdeltaEtaEleClusterTrackAtCalo_;
  vector<float> gsfdeltaEtaSeedClusterTrackAtVtx_;
  vector<float> gsfdeltaPhiSuperClusterTrackAtVtx_;
  vector<float> gsfdeltaPhiSeedClusterTrackAtCalo_;
  vector<float> gsfdeltaPhiEleClusterTrackAtCalo_;
  vector<float> gsfdeltaErawEcorrOverEcorr_;
  vector<float> gsfhcalOverEcal_;
  vector<float> gsfhcalOverEcalBc_;
  vector<float> gsfrawEnergy_;
  vector<float> gsfcorrectedEcalEnergy_;
  vector<float> gsfErawOvertrackP_;
  vector<float> gsfEcorrOvertrackP_;
  vector<int> gsfmomRHNtag_;
  vector<int> gsfmomWtag_;
  vector<int> gsfmomtag_;

  //Gsf tracks within dR<0.1 of the seedgsftrack

  vector<float> trkdR0p1_Mtrktrk_;
  

  //Gen Variables

  float         gendR_;

  //Gsf Track variables

  unsigned int ngsftrack_dR0p1_;

  //Gsfelectron and Gsftrack hybrid variables

  vector<float> gsf_PtdiffgsfgsftrackOverPtgsf_;
  
  
  ClassDef(MergedElectrons_ana,0);
  
};

#endif

#ifdef MergedElectrons_ana_cxx
void MergedElectrons_ana::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

  fReader              .SetTree(tree); //fReader is used to read all the common branches
  if(_data == 0)                        
    fReader_MC         .SetTree(tree); //If Input file is MC
  else if(_data == 1) 
    fReader_Data       .SetTree(tree); //If Input file is Data

   fReader.SetTree(tree);

   // Tree maker to selected branches for NN
   
   treefile_ = new TFile(TreeFileName_,"RECREATE");
   mytree = new TTree("myVariables","myVariables");

   //Filling integer numbers
   mytree->Branch("ngsf",          &ngsf_,          "ngsf/I"          );
   mytree->Branch("ngsftrack",     &ngsftrack_,     "ngsftrack/I"     );
   mytree->Branch("ngenelectrons", &ngenelectrons_, "ngenelectrons/I" );
   mytree->Branch("ngsf_momRHN",   &ngsf_momRHN_,   "ngsf_momRHN/I"   );

   //Event number filled for each electron if 2 electron then 2 event numbers (to identify them when dataframe is exploded in python NN)

   mytree->Branch("gsfEventno",                          "vector<int>",        &gsfeventno_                         );

   //Filling Direct and Derived Gsf variables

   mytree->Branch("gsfPt",                               "vector<float>",      &gsfPt_                              );
   mytree->Branch("gsfEta",                              "vector<float>",      &gsfEta_                             );
   mytree->Branch("gsfPhi",                              "vector<float>",      &gsfPhi_                             );
   mytree->Branch("gsfseedtrackPt",                      "vector<float>",      &gsfseedtrackPt_                     );
   mytree->Branch("gsfseedtrackEta",                     "vector<float>",      &gsfseedtrackEta_                    );
   mytree->Branch("gsfseedtrackPhi",                     "vector<float>",      &gsfseedtrackPhi_                    );
   mytree->Branch("gsfseedtrackP",                       "vector<float>",      &gsfseedtrackP_                      );
   mytree->Branch("gsfsigmaEtaEta",                      "vector<float>",      &gsfsigmaEtaEta_                     );
   mytree->Branch("gsfsigmaIetaIeta",                    "vector<float>",      &gsfsigmaIetaIeta_                   );
   mytree->Branch("gsfsigmaIphiIphi",                    "vector<float>",      &gsfsigmaIphiIphi_                   );
   mytree->Branch("gsfr9",                               "vector<float>",      &gsfr9_                              );
   mytree->Branch("gsfeSuperClusterOverP",               "vector<float>",      &gsfeSuperClusterOverP_              );
   mytree->Branch("gsfeSeedClusterOverP",                "vector<float>",      &gsfeSeedClusterOverP_               );
   mytree->Branch("gsfdeltaEtaSuperClusterTrackAtVtx",   "vector<float>",      &gsfdeltaEtaSuperClusterTrackAtVtx_  );
   mytree->Branch("gsfdeltaEtaSeedClusterTrackAtCalo",   "vector<float>",      &gsfdeltaEtaSeedClusterTrackAtCalo_  );
   mytree->Branch("gsfdeltaEtaEleClusterTrackAtCalo",    "vector<float>",      &gsfdeltaEtaEleClusterTrackAtCalo_   );
   mytree->Branch("gsfdeltaEtaSeedClusterTrackAtVtx",    "vector<float>",      &gsfdeltaEtaSeedClusterTrackAtVtx_   );
   mytree->Branch("gsfdeltaPhiSuperClusterTrackAtVtx",   "vector<float>",      &gsfdeltaPhiSuperClusterTrackAtVtx_  );
   mytree->Branch("gsfdeltaPhiSeedClusterTrackAtCalo",   "vector<float>",      &gsfdeltaPhiSeedClusterTrackAtCalo_  );
   mytree->Branch("gsfdeltaPhiEleClusterTrackAtCalo",    "vector<float>",      &gsfdeltaPhiEleClusterTrackAtCalo_   );
   mytree->Branch("gsfdeltaErawEcorrOverEcorr",          "vector<float>",      &gsfdeltaErawEcorrOverEcorr_         );
   mytree->Branch("gsfhcalOverEcal",                     "vector<float>",      &gsfhcalOverEcal_                    );
   mytree->Branch("gsfhcalOverEcalBc",                   "vector<float>",      &gsfhcalOverEcalBc_                  );
   mytree->Branch("gsfrawEnergy",                        "vector<float>",      &gsfrawEnergy_                       );
   mytree->Branch("gsfcorrectedEcalEnergy",              "vector<float>",      &gsfcorrectedEcalEnergy_             );
   mytree->Branch("gsfErawOvertrackP",                   "vector<float>",      &gsfErawOvertrackP_                  );
   mytree->Branch("gsfEcorrOvertrackP",                  "vector<float>",      &gsfEcorrOvertrackP_                 );
   mytree->Branch("gsfmomRHNtag",                        "vector<int>",        &gsfmomRHNtag_                       );
   mytree->Branch("gsfmomWtag",                          "vector<int>",        &gsfmomWtag_                         );
   mytree->Branch("gsfmomtag",                           "vector<int>",        &gsfmomtag_                          );

   //Filling gen info

   mytree->Branch("gendR",&gendR_,"gendR/F");

   //Filling the Gsf track info

   mytree->Branch("ngsftrack_dR0p1", &ngsftrack_dR0p1_, "ngsftrack_dR0p1/I" );

   //Gsfelectron and Gsftrack hybrid variables

   mytree->Branch("gsf_PtdiffgsfgsftrackOverPtgsf",  "vector<float>",      &gsf_PtdiffgsfgsftrackOverPtgsf_    );

   //Gsf tracks within dR<0.1 of the seedgsftrack

   mytree->Branch("trkdR0p1_Mtrktrk",                 "vector<float>",      &trkdR0p1_Mtrktrk_                  );
   
     
}

Bool_t MergedElectrons_ana::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef MergedElectrons_ana_cxx
