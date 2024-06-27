#define MergedElectrons_ana_cxx
// The class definition in MergedElectrons_ana.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("MergedElectrons_ana.C")
// root> T->Process("MergedElectrons_ana.C","some options")
// root> T->Process("MergedElectrons_ana.C+")
//


#include "MergedElectrons_ana.h"
#include <TH2.h>
#include <TStyle.h>

void MergedElectrons_ana::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
}

void MergedElectrons_ana::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

   //Initialization of the counters

   nEvtTotal     = 0;
   nEvtpass      = 0;
   nrecopass     = 0;
   nonereco      = 0;
   ntrackpass    = 0;
   nrecodRpass   = 0;
   ntrackdRpass  = 0;
   nrecodRpass2  = 0;
   ntrackdRpass2 = 0;
   ngendRl0p1    = 0;
   ngen_dRl0p1   = 0;
   ngen_dRg0p2   = 0;
   n0reco_dRl0p1 = 0;
   n1reco_dRl0p1 = 0;
   n2reco_dRl0p1 = 0;
   n0reco_dRg0p2 = 0;
   n1reco_dRg0p2 = 0;
   n2reco_dRg0p2 = 0;

   _HstFile = new TFile(_HstFileName,"recreate");
   BookHistograms();

}

Bool_t MergedElectrons_ana::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // When processing keyed objects with PROOF, the object is already loaded
   // and is available via the fObject pointer.
   //
   // This function should contain the \"body\" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

   fReader.SetLocalEntry(entry);
   if(_data == 0)
     fReader_MC.SetLocalEntry(entry);
   if(_data == 1)
     fReader_Data.SetLocalEntry(entry);
   
   nEvtTotal++;

   /**************************************************************************************************************************************************************************************************
                 
                                     Weights should be applied on Z sample(1) as we reweight Z pT to JPsi pT. So remove it when running JPsi sample(0)
				     
   *************************************************************************************************************************************************************************************************/   

   track.clear();
   gsftrack.clear();
   gsfseedtrack.clear();
   gsftrack_dR0p1.clear();
   genElec.clear();
   genElec_momRHN.clear();
   genElec_momW.clear();
   gsfElec.clear();
   gsfElec_momRHN.clear();
   gsfElec_momW.clear();
   genJPsi.clear();
   genPart.clear();
   genRHN.clear();

   
   //Weights for pT 0 to 50 GeV
   float w1[50] = {2.8,0.875,4.63158,3.78409,2.95,2.86667,2.14407,1.75185,1.56418,1.36508,1.33333,1.27696,1.28087,1.0552,1.12391,1.0383,1.07835,0.865832,0.869888,0.83213,0.84657,0.787879,0.77836,0.847518,0.766779,0.729055,0.692557,0.698305,0.755814,0.833333,0.700508,0.718802,0.711111,0.700173,0.801802,0.769084,0.793169,0.863905,0.800377,0.715044,0.697674,0.877016,0.68097,0.702555,0.653992,0.655678,0.742393,0.630597,0.57485,0.66666};
   //Weights for pT 50 to 250 GeV
   float w2[40] = {0.597631,0.583301,0.507107,0.4998,0.467422,0.467314,0.483051,0.521042,0.496489,0.542565,0.531975,0.527693,0.52752,0.549677,0.610253,0.581148,0.588598,0.589506,0.652529,0.620631,0.685506,0.769307,0.688689,0.786243,0.752464,0.829664,0.818849,0.886889,0.957004,1.0956,1.03566,1.12727,1.21964,1.40113,1.3457,1.36454,1.51209,1.677,1.89017,2.50362};
   
   h.nGsf->Fill(*ngsf);

       //Condition only for JPsi sample as many has only one electron events whereas in reality they are two electron events( No more used)
       for(unsigned int i=0;i<(*ngsf);i++){
	 Lepton temp;
	 temp.v.SetPtEtaPhiM(gsfPt[i],gsfEta[i],gsfPhi[i],0.000511);
	 temp.rawEnergy = gsfrawEnergy[i];
	 temp.correctedEnergy = gsfcorrectedEcalEnergy[i];
	 temp.seedtrkpx = gsfseedtrackPx[i];
	 temp.seedtrkpy = gsfseedtrackPy[i];
	 temp.seedtrkpz = gsfseedtrackPz[i];
	 temp.seedtrkpT = gsfseedtrackPt[i];
	 temp.seedtrkEta = gsfseedtrackEta[i];
	 temp.seedtrkPhi = gsfseedtrackPhi[i];
	 temp.pT =  gsfPt[i];
	 temp.eta = gsfEta[i];
	 temp.phi = gsfPhi[i];
	 temp.px = gsfPx[i];
	 temp.py = gsfPy[i];
	 temp.pz = gsfPz[i];
	 gsfElec.push_back(temp);
	 temp.v.SetPtEtaPhiM(gsfseedtrackPt[i],gsfseedtrackEta[i],gsfseedtrackPhi[i],0.000511);
	 gsfseedtrack.push_back(temp);
	 float px = gsfPx[i];
	 float py = gsfPy[i];
	 float pz = gsfPz[i];
	 float p = TMath::Sqrt(px*px+py*py+pz*pz);
	 float pT = gsfPt[i];
	 float Eta = gsfEta[i];
	 float sigmaEtaEta = gsfsigmaEtaEta[i];
	 float sigmaIetaIeta = gsfsigmaIetaIeta[i];
	 float sigmaIphiIphi = gsfsigmaphiIphi[i];
	 float e1x5 = gsfe1x5[i];
	 float e2x5Max = gsfe2x5Max[i];
	 float e5x5 = gsfe5x5[i];
	 float r9 = gsfr9[i];
	 float eSuperClusterOverP = gsfeSuperClusterOverP[i];
	 float E_raw = gsfrawEnergy[i];
	 float E_corr = gsfcorrectedEcalEnergy[i];
	 float dE_raw_overE_corr = fabs(E_raw-E_corr)/E_corr;
	 float Eraw_overP = E_raw/p;
	 float Ecorr_overP = E_corr/p;
	 float dr03TkSumPt = gsfdr03TkSumPt[i];
	 float dr03TkSumPtHEEP = gsfdr03TkSumPtHEEP[i];
	 float dr04TkSumPt = gsfdr04TkSumPt[i];
	 float dr04TkSumPtHEEP = gsfdr04TkSumPtHEEP[i];
	 float detascTkatVtx = gsfdeltaEtaSuperClusterTrackAtVtx[i];
	 float detaseedTkatCalo = gsfdeltaEtaSeedClusterTrackAtCalo[i];
	 float dphiscTkatVtx = gsfdeltaPhiSuperClusterTrackAtVtx[i];
	 float dphiseedTkatCalo = gsfdeltaPhiSeedClusterTrackAtCalo[i];
	 float full5x5_sigmaEtaEta = gsffull5x5_sigmaEtaEta[i];
	 float full5x5_sigmaIetaIeta = gsffull5x5_sigmaIetaIeta[i];
	 float full5x5_sigmaIphiIphi = gsffull5x5_sigmaIphiIphi[i];
	 float full5x5_r9 = gsffull5x5_r9[i];
	 float full5x5_hcaloverEcal = gsffull5x5_hcalOverEcal[i];
	 float hcaloverEcal = gsfhcalOverEcal[i];
	 float Eraw_over_gsftrackp = E_raw/gsftrackpMode[i];
	 float var[28] = {pT,Eta,sigmaEtaEta,sigmaIetaIeta,sigmaIphiIphi,e1x5,e2x5Max,e5x5,r9,eSuperClusterOverP,dE_raw_overE_corr,Eraw_overP,Ecorr_overP,dr03TkSumPt,dr03TkSumPtHEEP,dr04TkSumPt,dr04TkSumPtHEEP,detascTkatVtx,detaseedTkatCalo,dphiscTkatVtx,dphiseedTkatCalo,full5x5_sigmaEtaEta,full5x5_sigmaIetaIeta,full5x5_sigmaIphiIphi,full5x5_r9,full5x5_hcaloverEcal,hcaloverEcal, Eraw_over_gsftrackp};

	 int varsize = sizeof(var)/sizeof(float);
     
	 for(int ivar=0;ivar<varsize;ivar++){
	   for(int ipt=0;ipt<50;ipt++){
	     if(pT>=ipt && pT<(ipt+1)) h.Gsfprop_histo[ivar]->Fill(var[ivar],w1[ipt]);
	   }
       
	   for(int ipt=10;ipt<50;ipt++){
	     if(pT>=5*(ipt) && pT<5*(ipt+1)) h.Gsfprop_histo[ivar]->Fill(var[ivar],w2[ipt-10]);
	   }
	 }

	 if(fabs(Eta)<1.4) h.misc_histo[0]->Fill(sigmaEtaEta);
	 if(fabs(Eta)>1.4 && fabs(Eta)<2.4) h.misc_histo[1]->Fill(sigmaEtaEta);

	 float ErawOvertrackP = gsfrawEnergy[i]/gsfseedtrackP[i];
	 float EcorrOvertrackP = gsfcorrectedEcalEnergy[i]/gsfseedtrackP[i];
	 float EcorrOverEraw = gsfcorrectedEcalEnergy[i]/gsfrawEnergy[i];
	 float dErawEcorrOverEraw = abs(gsfcorrectedEcalEnergy[i]-gsfrawEnergy[i])/(gsfrawEnergy[i]);
	 float HoverE = gsfhcalOverEcal[i];
	 float HoverEbc = gsfhcalOverEcalBc[i];
	 float pTdiffgsfgtOvergsf = (gsfPt[i]-gsfseedtrackPt[i])/gsfPt[i];
	 h.gsf_histo[11]->Fill(ErawOvertrackP);
	 h.gsf_histo[12]->Fill(EcorrOvertrackP);
	 h.gsf_histo[13]->Fill(EcorrOverEraw);
	 h.gsf_histo[14]->Fill(dErawEcorrOverEraw);
	 h.gsf_histo[15]->Fill(HoverE);
	 h.gsf_histo[16]->Fill(HoverEbc);
	 h.gsf_histo[44]->Fill(gsfdeltaEtaSuperClusterTrackAtVtx[i]);
	 h.gsf_histo[45]->Fill(gsfdeltaPhiSuperClusterTrackAtVtx[i]);
	 h.gsf_histo[46]->Fill(gsfr9[i]);
	 
	 if((*ngsf)==2){
	   h.gsf_histo[17]->Fill(gsfseedtrackPt[i]);
	   h.gsf_histo[18]->Fill(gsfPt[i]);
	 }
	 if((*ngsf)==1){
	   h.gsf_histo[19]->Fill(gsfseedtrackPt[i]);
	   h.gsf_histo[20]->Fill(gsfPt[i]);
	   h.gsf_histo[25]->Fill(ErawOvertrackP);
	   h.gsf_histo[26]->Fill(EcorrOvertrackP);
	   h.gsf_histo[27]->Fill(pTdiffgsfgtOvergsf);

	   //Plotting Some properties for MS thesis plots
	   
	   h.mergedgsf_histo[0]->Fill(gsfr9[i]);
	   h.mergedgsf_histo[1]->Fill(gsfsigmaEtaEta[i]);
	   h.mergedgsf_histo[3]->Fill(gsfsigmaIetaIeta[i]);
	   h.mergedgsf_histo[4]->Fill(gsfsigmaphiIphi[i]);
	   h.mergedgsf_histo[5]->Fill(gsfcorrectedEcalEnergy[i]/gsfrawEnergy[i]);
	   h.mergedgsf_histo[6]->Fill(gsfdeltaEtaSuperClusterTrackAtVtx[i]);
	   h.mergedgsf_histo[7]->Fill(gsfdeltaPhiSuperClusterTrackAtVtx[i]);
	   h.mergedgsf_histo[8]->Fill(ErawOvertrackP);
	   h.mergedgsf_histo[9]->Fill(EcorrOvertrackP);
	   
	 }
	 
       }
       h.dntrack_ngsftrack->Fill(fabs(*ntrack-*ngsftrack));
       h.dngsf_ngsftrack->Fill(fabs(*ngsf-*ngsftrack));

       //Track 
       
       for(unsigned int i=0;i<(*ntrack);i++){
	 Lepton temp;
	 float normalizedChi2 = tracknormalizedChi2[i];
	 temp.v.SetPtEtaPhiM(trackpt[i],tracketa[i],trackphi[i],0.00511);
	 track.push_back(temp);
	 h.trackprop_histo[0]->Fill(normalizedChi2);
       }

       for(unsigned int i=0;i<(*ngsftrack);i++){
	 Lepton temp;
	 temp.v.SetPtEtaPhiM(gsftrackptMode[i],gsftracketaMode[i],gsftrackphiMode[i],0.00511);
	 gsftrack.push_back(temp);
       }

       if(_data == 0){
       //GenParticle
	 for(unsigned int i=0;i<(*ngenparticle);i++){
	   Lepton temp;
	   int id = pdgid[i];
	   int status = genstatus[i];
	   int nmom = nmother[i];
	   float pT = genpT[i];
	   temp.id = pdgid[i];
	   temp.v.SetPtEtaPhiM(genpT[i],genEta[i],genPhi[i],genM[i]);
	   temp.vx = genvx[i];
	   temp.vy = genvy[i];
	   temp.vz = genvz[i];
	   if(fabs(pdgid[i])==11 && genstatus[i]==1){
	     genElec.push_back(temp);
	   }
	   if(fabs(pdgid[i])==443){
	     genJPsi.push_back(temp);
	   }
	   genPart.push_back(temp);
	   if(fabs(pdgid[i])==9900012) genRHN.push_back(temp);

	   //Extracting Mother information
	    if(nmom>0 && status==1 && fabs(id)==11 && fabs(momid[i])==9900012) genElec_momRHN.push_back(temp);
	    if(nmom>0 && status==1 && fabs(id)==11 && (fabs(momid[i]) == 24 || fabs(momid[i])<=6)) genElec_momW.push_back(temp);

	    if(status==1 && fabs(id)==11) h.genElec_histo[18]->Fill(momid[i]);
	 }
       }
       
       h.genRHN_histo[0]->Fill(genRHN.size());
       
       h.genElec_histo[0]->Fill(genElec.size());
       h.genElec_histo[2]->Fill(genJPsi.size());
       h.genElec_histo[16]->Fill(genElec_momRHN.size());
       h.genElec_histo[17]->Fill(genElec_momW.size());


       /*******************************************************************************************************************************************************************************************

                                                                                         Analysis Part of the Code

       ******************************************************************************************************************************************************************************************/      

       sort(0);
       
       for(unsigned int i=0;i<gsftrack.size();i++){
	 int n1 = 0;
	 int n2 = 0;
	 for(unsigned int j=0;j<track.size();j++){
	   float dR = gsftrack.at(i).v.DeltaR(track.at(j).v);
	   if(dR<0.1) n1++;
	   if(dR<0.05) n2++;
	 }
	 h.ntrack_dR01->Fill(n1);
	 h.ntrack_dR0p05->Fill(n2);
       }

       //Gsftracks and KF tracks within dR<0.1 of a Gsf electron 
       
       for(unsigned int i=0;i<gsfseedtrack.size();i++){
	 gsftrack_dR0p1.clear();
	 gsftrack_dR0p05.clear();
	 float dRmin0 = 999;
	 float dRmin1 = 999;
	 int ind0 = -1;
	 int ind1 = -1;
	 for(unsigned int j=0;j<gsftrack.size();j++){
	   float dR = gsfseedtrack.at(i).v.DeltaR(gsftrack.at(j).v);
	   if(dR < dRmin0){
	     dRmin1 = dRmin0;
	     ind1 = ind0;
	     dRmin0 = dR;
	     ind0 = j;
	   }
	   if(dR<dRmin1 && dR>dRmin0){
	     dRmin1 = dR;
	     ind1 = j;
	   }
	 }

	 if(dRmin0<0.1 && dRmin1<0.1){
	   gsftrack_dR0p1.push_back(gsftrack.at(ind0));
	   gsftrack_dR0p1.push_back(gsftrack.at(ind1));
	 }

	 if(dRmin0<0.05 && dRmin1<0.05){
	   gsftrack_dR0p05.push_back(gsftrack.at(ind0));
	   gsftrack_dR0p05.push_back(gsftrack.at(ind1));
	 }

	 if(dRmin0<0.1 && dRmin1>0.1) gsftrack_dR0p1.push_back(gsftrack.at(ind0));
	 if(dRmin0<0.05 && dRmin1>0.05) gsftrack_dR0p05.push_back(gsftrack.at(ind0));
  
	 h.misc_histo[3]->Fill(gsftrack_dR0p1.size());
	 h.misc_histo[7]->Fill(gsftrack_dR0p05.size());
	 
	 //Plotting the invariant mass of the two Gsf tracks

	 if(gsftrack_dR0p1.size()>1){
	   float Mgtgt = (gsftrack_dR0p1.at(0).v+gsftrack_dR0p1.at(1).v).M();
	   float dRgtgt = gsftrack_dR0p1.at(0).v.DeltaR(gsftrack_dR0p1.at(1).v);
	   float pTgtgt = (gsftrack_dR0p1.at(0).v+gsftrack_dR0p1.at(1).v).Pt();
	   float pTgt0  = gsftrack_dR0p1.at(0).v.Pt();
	   float pTgt1  = gsftrack_dR0p1.at(1).v.Pt();
	   h.gsftrack_histo[0]->Fill(Mgtgt);
	   h.gsftrack_histo[1]->Fill(dRgtgt);
	   h.gsftrack_histo[3]->Fill(pTgtgt);
	   if((*ngsf)==1){
	     float pTgsf = gsfElec.at(0).v.Pt();
	     float pTgsfovergt = pTgsf/pTgtgt;
	     float pTasym      = fabs(pTgt0-pTgt1)/pTgt0;
	     h.gsftrack_histo[4]->Fill(pTgtgt);
	     h.gsftrack_histo[5]->Fill(pTgsfovergt);
	     if(pTgsfovergt>0.9 && pTgsfovergt<1.1){
	       
	       h.gsftrack_histo[6]->Fill(pTgtgt);
	       h.gsftrack_histo[7]->Fill(pTgsf);
	       h.gsftrack_histo[8]->Fill((pTgsf-pTgtgt));
	       h.gsftrack_histo[9]->Fill((pTgt0-pTgt1));
	       h.gsftrack_histo[10]->Fill(pTasym);
	     }
	     if(genElec.size()>1){
	       float pTJPsi = (genElec.at(0).v+genElec.at(1).v).Pt();
	       float pTdiff = (pTJPsi-pTgtgt);
	       float pTdiffoverpTJPsi = (pTdiff/pTJPsi);

	       h.gsftrack_histo[11]->Fill(pTdiff);
	       h.gsftrack_histo[12]->Fill(pTdiffoverpTJPsi);
	     }
	     
	   }
	   if(genElec.size()>1){
	     float Mee = (genElec.at(0).v+genElec.at(1).v).M();
	     float dRgen = genElec.at(0).v.DeltaR(genElec.at(1).v);

	     if(dRgen<0.03 && dRgtgt<0.03){
	       float Mreso = (Mee-Mgtgt)/Mee;
	       h.misc_histo[4]->Fill(Mreso);
	     }

	     if(dRgen<0.06 && dRgen>0.03 && dRgtgt<0.06 && dRgtgt>0.03){
	       float Mreso = (Mee-Mgtgt)/Mee;
	       h.misc_histo[5]->Fill(Mreso);
	     }

	     if(dRgen<0.1 && dRgen>0.06 && dRgtgt<0.1 && dRgtgt>0.06){
	       float Mreso = (Mee-Mgtgt)/Mee;
	       h.misc_histo[6]->Fill(Mreso);
	     }
	   }  
	 }

	 //***************** One reco electron + two gsftracks **************************//
	 
	 if(gsfElec.size()==1 && gsftrack_dR0p1.size()>1){
	   ntrackpass++;
	   TLorentzVector vec1,vec2;
	   float rawE = gsfElec.at(0).rawEnergy;
	   float corrE = gsfElec.at(0).v.E();
	   float px = gsftrack_dR0p1.at(0).v.Px()+gsftrack_dR0p1.at(1).v.Px();
	   float py = gsftrack_dR0p1.at(0).v.Py()+gsftrack_dR0p1.at(1).v.Py();
	   float pz = gsftrack_dR0p1.at(0).v.Pz()+gsftrack_dR0p1.at(1).v.Pz();
	   vec1.SetPxPyPzE(px,py,pz,rawE);
	   vec2.SetPxPyPzE(px,py,pz,corrE);
	   float Mcorr = TMath::Sqrt((corrE)*(corrE)-px*px-py*py-pz*pz);
	   float Mraw =  TMath::Sqrt((rawE)*(rawE)-px*px-py*py-pz*pz);
	   float Mcorr2 = (corrE)*(corrE)-px*px-py*py-pz*pz;
	   float Mraw2 =  (rawE)*(rawE)-px*px-py*py-pz*pz;
	   
	   h.gsf_histo[4]->Fill(Mraw);
	   h.gsf_histo[5]->Fill(Mcorr);
	   h.gsf_histo[6]->Fill(Mraw2);
	   h.gsf_histo[7]->Fill(Mcorr2);
	   
	 }
	 if(gsfElec.size()==1 && gsftrack_dR0p1.size()>1 && genElec.size()>1){
	   float dR = genElec.at(0).v.DeltaR(genElec.at(1).v);
	   if(dR<0.1) ntrackdRpass++;
	   if(dR>0.2) ntrackdRpass2++;
	 }
	 
       }

       
       if(gsfElec.size()>1){
	 float dR = gsfElec.at(0).v.DeltaR(gsfElec.at(1).v);
	 h.gsf_histo[0]->Fill(dR);
       }
       if(gsftrack.size()>1){
	 float dR = gsftrack.at(0).v.DeltaR(gsftrack.at(1).v);
	 h.gsftrack_histo[2]->Fill(dR);
       }

       
       for(unsigned int i=0;i<gsfElec.size();i++){
	 float dRmin = 999;
	 for(unsigned int j=0;j<genElec.size();j++){
	   float dR = genElec.at(j).v.DeltaR(gsfElec.at(i).v);
	   if(dR<dRmin) dRmin = dR;
	 }
       }

       //Genparticle analysis

       if(genElec.size()>1){
	 float pT_JPsi = (genElec.at(0).v+genElec.at(1).v).Pt();
	 float dR = genElec.at(0).v.DeltaR(genElec.at(1).v);
	 float Mee = (genElec.at(0).v+genElec.at(1).v).M();
	 h.gen_histo[0]->Fill(dR,pT_JPsi);
	 h.genElec_histo[1]->Fill(Mee);
	 h.genElec_histo[3]->Fill(dR);
	 if(dR<0.1) ngendRl0p1++;
       }
       

       //Two reco electrons
       if(gsfElec.size()>1 && genElec.size()==2){
	 nrecopass++;
	 float dR = genElec.at(0).v.DeltaR(genElec.at(1).v);
	 if(dR<0.1) nrecodRpass++;
	 if(dR>0.2) nrecodRpass2++;
       }
       if(genElec.size()==2) nEvtpass++;

       //Finding the ratios alpha and beta
       
       if(genElec.size()>1){
         float dR = genElec.at(0).v.DeltaR(genElec.at(1).v);
	 float E0 = genElec.at(0).v.E();
	 float E1 = genElec.at(1).v.E();
	 float pT0 = genElec.at(0).v.Pt();
	 float pT1 = genElec.at(1).v.Pt();
	 float pTJPsi = (genElec.at(0).v+genElec.at(1).v).Pt();
	 float Etot = E0+E1;
	 float alpha = E0/(E0+E1);
	 float beta = E1/(E0+E1);
	 float alpha_pT = pT0/(pT0+pT1);
	 float beta_pT = pT1/(pT0+pT1);
	 float pT_reso = (pT0-pT1)/pT0;

	 h.genElec_histo[4]->Fill(alpha);
	 h.genElec_histo[5]->Fill(beta);
	 h.genElec_histo[6]->Fill(alpha+beta);
	 h.genElec_histo[7]->Fill(alpha_pT);
	 h.genElec_histo[8]->Fill(beta_pT);
	 h.genElec_histo[9]->Fill(pT_reso);
	 h.genElec_histo[14]->Fill(pTJPsi);
       }

       //Finding the ratios of electrons to raw energies

       for(unsigned int i=0;i<gsfElec.size();i++){
	 float alpha = (gsfElec.at(i).v.E())/(gsfElec.at(i).rawEnergy);
	 h.gsf_histo[1]->Fill(alpha);
	 h.gsf_histo[9]->Fill(gsfElec.at(i).rawEnergy);
	 h.gsf_histo[10]->Fill(gsfElec.at(i).correctedEnergy);
       }

       h.gsf_histo[2]->Fill(gsfElec.size());
       
       if(gsfElec.size()>1){
	 float Mee = (gsfElec.at(0).v+gsfElec.at(1).v).M();
	 h.gsf_histo[3]->Fill(Mee);
       }

       if(gsfElec.size()==1 && genElec.size()>1){
	 float E = genElec.at(0).v.E()+genElec.at(1).v.E();
	 float rawE = gsfElec.at(0).rawEnergy;
	 float corrE = gsfElec.at(0).v.E();
	 h.genElec_histo[10]->Fill(rawE/E);
	 h.genElec_histo[11]->Fill(corrE/E);
       }

       for(unsigned int i=0;i<gsfElec.size();i++){
	 float dRmin = 999;
	 for(unsigned int j=0;j<genElec.size();j++){
	   float dR = gsfElec.at(i).v.DeltaR(genElec.at(j).v);
	   if(dR<dRmin) dRmin = dR;
	 }
	 h.gsf_histo[8]->Fill(dRmin);
       }

       for(unsigned int i=0;i<genPart.size();i++){
	 if(abs(genPart.at(i).id)>1000) h.genPart_histo[0]->Fill(abs(genPart.at(i).id));
       }
       
       int nhadrons = 0;
       int nphotons = 0;
       for(unsigned int i=0;i<genElec.size();i++){
	 float dxy = TMath::Sqrt(genElec.at(i).vx*genElec.at(i).vx + genElec.at(i).vy*genElec.at(i).vy);
	 float dz = genElec.at(i).vz;
	 float dRmin = 999;
	 int ind = -1;
	 for(unsigned int j=0;j<genPart.size();j++){
	   if(fabs(genPart.at(j).id != 11)){
	     float dR = genElec.at(i).v.DeltaR(genPart.at(j).v);
	     if(dR<dRmin){
	       dRmin = dR;
	       ind = j;
	     }
	   }
	 }
	 h.genElec_histo[12]->Fill(dxy);
	 h.genElec_histo[13]->Fill(dz);
	 
	 h.genPart_histo[1]->Fill(dRmin);
	 if(dRmin<0.1 && fabs(genPart.at(ind).id)>30) nhadrons++;
	 if(dRmin<0.1 && fabs(genPart.at(ind).id)==22) nphotons++;
       }

       if(gsfElec.size()==1) nonereco++;
       
       h.genPart_histo[2]->Fill(nhadrons);
       h.genPart_histo[3]->Fill(nphotons);

       //Analysis with Gsf electron and closest Gsftrack

       for(unsigned int i=0;i<gsfElec.size();i++){
	 float dRmin = 999;
	 int ind = -1;
	 for(unsigned int j=0;j<gsftrack.size();j++){
	   float dR = gsfElec.at(i).v.DeltaR(gsftrack.at(j).v);
	   if(dR<dRmin){
	     dRmin = dR;
	     ind = j;
	   }	   
	 }
	 float pTgsf = gsfElec.at(i).v.Pt();
	 float pTgt  = gsftrack.at(ind).v.Pt();
	 if(dRmin<0.1 && gsfElec.size()==1){
	   h.gsf_histo[21]->Fill(pTgsf-pTgt);
	   h.gsf_histo[23]->Fill((pTgsf-pTgt)/pTgsf);
	 }
	 if(dRmin<0.1 && gsfElec.size()>1){
	   h.gsf_histo[22]->Fill(pTgsf-pTgt);
	   h.gsf_histo[24]->Fill((pTgsf-pTgt)/pTgsf);
	 }
       }


       /*******************************************************************************************************************************************************************************************

                                                                     Some more analysis using Gsf electrons Calo and Gsftrack properties

       ********************************************************************************************************************************************************************************************/

       
       for(unsigned int i=0;i<gsfElec.size();i++){
	 TLorentzVector v;
	 float trkpt = gsfElec.at(i).seedtrkpT;
	 float trketa = gsfElec.at(i).seedtrkEta;
	 float trkphi = gsfElec.at(i).seedtrkPhi;
	 float pT     = gsfElec.at(i).pT;
	 float eta    = gsfElec.at(i).eta;
	 float phi    = gsfElec.at(i).phi;
	 float Eraw  = gsfElec.at(i).rawEnergy;          
	 float Ecorr = gsfElec.at(i).correctedEnergy;

	 TLorentzVector vtrkraw,vtrkcorr,vraw,vcorr;
	 vtrkraw.SetPtEtaPhiE(trkpt,trketa,trkphi,Eraw);
	 vtrkcorr.SetPtEtaPhiE(trkpt,trketa,trkphi,Ecorr);
	 vraw.SetPtEtaPhiE(pT,eta,phi,Eraw);
	 vcorr.SetPtEtaPhiE(pT,eta,phi,Ecorr);

	 h.gsf_histo[28]->Fill(vtrkraw.M());
	 h.gsf_histo[29]->Fill(vtrkcorr.M());
	 h.gsf_histo[30]->Fill(vraw.M());
	 h.gsf_histo[31]->Fill(vcorr.M());

	 if(gsfElec.size()==1){
	   h.gsf_histo[32]->Fill(vtrkraw.M());
	   h.gsf_histo[33]->Fill(vtrkcorr.M());
	   h.gsf_histo[34]->Fill(vraw.M());
	   h.gsf_histo[35]->Fill(vcorr.M());
	 }

	 //Now let us do the same exercise above with gen Matching
	 float dRmin = 999;
	 for(unsigned int j=0;j<genElec.size();j++){
	   h.genElec_histo[15]->Fill(genElec.at(j).v.M());
	   float dR = genElec.at(j).v.DeltaR(gsfElec.at(i).v);
	   if(dR<dRmin) dRmin = dR;
	 }

	 if(dRmin<0.05){
	   h.gsf_histo[36]->Fill(vtrkraw.M());
	   h.gsf_histo[37]->Fill(vtrkcorr.M());
	   h.gsf_histo[38]->Fill(vraw.M());
	   h.gsf_histo[39]->Fill(vcorr.M());
	   
	   if(gsfElec.size()==1){
	     h.gsf_histo[40]->Fill(vtrkraw.M());
	     h.gsf_histo[41]->Fill(vtrkcorr.M());
	     h.gsf_histo[42]->Fill(vraw.M());
	     h.gsf_histo[43]->Fill(vcorr.M());
	   }

	 }
	 
       }

       /********************************************************************************************************************************************************************************************

                                                                      Calculating the Numbers of reconstrcuted electrons for gendR<0.1 and gendR>0.2
       
       *********************************************************************************************************************************************************************************************/

       if(genElec.size()>1){
	 float dR = genElec.at(0).v.DeltaR(genElec.at(1).v);
	 
	 if(dR<0.1){
	   ngen_dRl0p1++;
	   if(gsfElec.size() == 0) n0reco_dRl0p1++;
	   if(gsfElec.size() == 1) n1reco_dRl0p1++;
	   if(gsfElec.size() == 2) n2reco_dRl0p1++;
	 }
	 if(dR>0.1 && dR<0.2){
	   ngen_dRg0p2++;
	   if(gsfElec.size() == 0) n0reco_dRg0p2++;
	   if(gsfElec.size() == 1) n1reco_dRg0p2++;
	   if(gsfElec.size() == 2) n2reco_dRg0p2++;
	 }
	 
       }

       /******************************************************************************************************************************************************************************************

                                                                      RHN specific analysis Here. Both Gen electrons with mother RHN and reco electrons matched and analysed

       *******************************************************************************************************************************************************************************************/

       //Getting dR between gen electrons from W and RHN

       for(unsigned int i=0;i<genElec_momW.size();i++){
	 for(unsigned int j=0;j<genElec_momRHN.size();j++){
	   float dR = genElec_momW.at(i).v.DeltaR(genElec_momRHN.at(j).v);
	   h.genElec_histo[19]->Fill(dR);
	 }
       }
       
       for(unsigned int i=0;i<gsfElec.size();i++){
	 float dRminRHN = 999;
	 float dRminW   = 999;
	 for(unsigned int j=0;j<genElec_momRHN.size();j++){
	   float dR = gsfElec.at(i).v.DeltaR(genElec_momRHN.at(j).v);
	   if(dR<dRminRHN) dRminRHN = dR;	   
	 }

	 for(unsigned int j=0;j<genElec_momW.size();j++){
	   float dR = gsfElec.at(i).v.DeltaR(genElec_momW.at(j).v);
	   if(dR<dRminW) dRminW = dR;
	 }

	 h.misc_histo[8]->Fill(dRminRHN);
	 h.misc_histo[9]->Fill(dRminW);
	 
	 if(dRminRHN<0.1 && genElec_momRHN.size()==2){
	   gsfElec.at(i).momtag = 1;
	   gsfElec_momRHN.push_back(gsfElec.at(i));
	 }
	 else if(dRminW<0.1){
	   gsfElec.at(i).momtag = 2;
	   gsfElec_momW.push_back(gsfElec.at(i));
	 }

	 else gsfElec.at(i).momtag = 3;
    
       }

       h.gsf_histo[47]->Fill(gsfElec_momRHN.size());
       h.gsf_histo[48]->Fill(gsfElec_momW.size());
              
       
       /********************************************************************************************************************************************************************************************

                                                                      Declaring and selecting the variables for Filling the trees
       
       *********************************************************************************************************************************************************************************************/

       //Fill The trees and all only if we have atleast one recostructed Gsf electron in an event
   
       if(gsfElec.size()>0){
       
	 //Integer Variables
	 ngsf_          = gsfElec.size();
	 ngsftrack_     = gsftrack.size();
	 ngenelectrons_ = genElec.size();
	 ngsf_momRHN_   = gsfElec_momRHN.size();

	 //Gsf variables
	 gsfeventno_.clear();
	 gsfPt_.clear();
	 gsfEta_.clear();
	 gsfPhi_.clear();
	 gsfseedtrackPt_.clear();
	 gsfseedtrackEta_.clear();
	 gsfseedtrackPhi_.clear();
	 gsfseedtrackP_.clear();
	 gsfsigmaEtaEta_.clear();
	 gsfsigmaIetaIeta_.clear();
	 gsfsigmaIphiIphi_.clear();
	 gsfr9_.clear();
	 gsfeSuperClusterOverP_.clear();
	 gsfeSeedClusterOverP_.clear();
	 gsfdeltaEtaSuperClusterTrackAtVtx_.clear();
	 gsfdeltaEtaSeedClusterTrackAtCalo_.clear();
	 gsfdeltaEtaEleClusterTrackAtCalo_.clear();
	 gsfdeltaEtaSeedClusterTrackAtVtx_.clear();
	 gsfdeltaPhiSuperClusterTrackAtVtx_.clear();
	 gsfdeltaPhiSeedClusterTrackAtCalo_.clear();
	 gsfdeltaPhiEleClusterTrackAtCalo_.clear();
	 gsfdeltaErawEcorrOverEcorr_.clear();
	 gsfhcalOverEcal_.clear();
	 gsfhcalOverEcalBc_.clear();
	 gsfrawEnergy_.clear();
	 gsfcorrectedEcalEnergy_.clear();
	 gsfErawOvertrackP_.clear();
	 gsfEcorrOvertrackP_.clear();
	 gsf_PtdiffgsfgsftrackOverPtgsf_.clear();
	 gsfmomRHNtag_.clear();
	 gsfmomWtag_.clear();
	 gsfmomtag_.clear();

	 //Gsf tracks within dR<0.1 of the seedgsftrack
	 
	 trkdR0p1_Mtrktrk_.clear();
	 
	 for(unsigned int i=0; i<(*ngsf); i++){

	   //Storing the size of gsfElec_momRHN for each gsf electrons
	   gsfmomRHNtag_.push_back(gsfElec_momRHN.size());
	   gsfmomWtag_.push_back(gsfElec_momW.size());
	   gsfmomtag_.push_back(gsfElec.at(i).momtag);
	   
	   gsfeventno_.push_back(nEvtTotal); //Storing Event Numbers
	   gsfPt_.push_back(gsfPt[i]);
	   gsfEta_.push_back(gsfEta[i]);
	   gsfPhi_.push_back(gsfPhi[i]);
	   gsfseedtrackPt_.push_back(gsfseedtrackPt[i]);
	   gsfseedtrackEta_.push_back(gsfseedtrackEta[i]);
	   gsfseedtrackPhi_.push_back(gsfseedtrackPhi[i]);
	   gsfseedtrackP_.push_back(gsfseedtrackP[i]);
	   gsfsigmaEtaEta_.push_back(gsfsigmaEtaEta[i]);
	   gsfsigmaIetaIeta_.push_back(gsfsigmaIetaIeta[i]);
	   gsfsigmaIphiIphi_.push_back(gsfsigmaphiIphi[i]);
	   gsfr9_.push_back(gsfr9[i]);
	   gsfeSuperClusterOverP_.push_back(gsfeSuperClusterOverP[i]);
	   gsfeSeedClusterOverP_.push_back(gsfeSeedClusterOverP[i]);
	   gsfdeltaEtaSuperClusterTrackAtVtx_.push_back(gsfdeltaEtaSuperClusterTrackAtVtx[i]);
	   gsfdeltaEtaSeedClusterTrackAtCalo_.push_back(gsfdeltaEtaSeedClusterTrackAtCalo[i]);
	   gsfdeltaEtaEleClusterTrackAtCalo_.push_back(gsfdeltaEtaEleClusterTrackAtCalo[i]);
	   gsfdeltaEtaSeedClusterTrackAtVtx_.push_back(gsfdeltaEtaSeedClusterTrackAtVtx[i]);
	   gsfdeltaPhiSuperClusterTrackAtVtx_.push_back(gsfdeltaPhiSuperClusterTrackAtVtx[i]);
	   gsfdeltaPhiSeedClusterTrackAtCalo_.push_back(gsfdeltaPhiSeedClusterTrackAtCalo[i]);
	   gsfdeltaPhiEleClusterTrackAtCalo_.push_back(gsfdeltaPhiEleClusterTrackAtCalo[i]);
	   gsfdeltaErawEcorrOverEcorr_.push_back(fabs(gsfrawEnergy[i]-gsfcorrectedEcalEnergy[i])/(gsfcorrectedEcalEnergy[i]));
	   gsfhcalOverEcal_.push_back(gsfhcalOverEcal[i]);
	   gsfhcalOverEcalBc_.push_back(gsfhcalOverEcalBc[i]);
	   gsfrawEnergy_.push_back(gsfrawEnergy[i]);
	   gsfcorrectedEcalEnergy_.push_back(gsfcorrectedEcalEnergy[i]);
	   gsfErawOvertrackP_.push_back(gsfrawEnergy[i]/gsfseedtrackP[i]);
	   gsfEcorrOvertrackP_.push_back(gsfcorrectedEcalEnergy[i]/gsfseedtrackP[i]);
	   gsf_PtdiffgsfgsftrackOverPtgsf_.push_back((gsfPt[i]-gsfseedtrackPt[i])/gsfPt[i]);

	   //Calculation for tracks within dR<0.1 of seedgsftrack (Doing this in electron loop for consistency, for each gsf electron there should be 1 object)

	   TLorentzVector vseedtrk;
	   vseedtrk.SetPtEtaPhiM(gsfseedtrackPt[i],gsfseedtrackEta[i],gsfseedtrackPhi[i],0.000511);
	   float dRmin0 = 999;
	   float dRmin1 = 999;
	   int ind0 = -1;
	   int ind1 = -1;
	   float M = -999;
	   for(unsigned int j=0;j<gsftrack.size();j++){
	     float dR = vseedtrk.DeltaR(gsftrack.at(j).v);
	     if(dR<dRmin0){
	       dRmin1 = dRmin0;
	       ind1 = ind0;
	       dRmin0 = dR;
	       ind0 = j;
	     }
	     if(dR<dRmin1 && dR>dRmin0){
	       dRmin1 = dR;
	       ind1 = j;
	     }
	   }
	   if(dRmin0<0.1 && dRmin1<0.1){
	     M = (gsftrack.at(ind0).v+gsftrack.at(ind1).v).M();
	   }
	   
	   trkdR0p1_Mtrktrk_.push_back(M);	   
	 }
	 
   	 //Gen Variables
	 if(genElec.size()>1){
	   gendR_ = genElec.at(0).v.DeltaR(genElec.at(1).v);
	 }
	 
	 //Gsf track variables
	 
	 int ngt = 0;
	 for(unsigned int i=0;i<gsfseedtrack.size();i++){
	   float dRmin0 = 999;
	   float dRmin1 = 999;
	   for(unsigned int j=0;j<gsftrack.size();j++){
	     float dR = gsfseedtrack.at(i).v.DeltaR(gsftrack.at(j).v);
	     if(dR<dRmin0){
	       dRmin1 = dRmin0;
	       dRmin0 = dR;
	     }
	     if(dR<dRmin1 && dR>dRmin0) dRmin1 = dR;
	   }

	   if(dRmin1<0.1 && dRmin0<0.1) ngt++;
	 }
	 ngsftrack_dR0p1_ = ngt;

	 mytree->Fill();
       }    
       
       return kTRUE;
}

void MergedElectrons_ana::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

  _HstFile->Write();
  _HstFile->Close();

  //Writing and closing the Tree and root file

  mytree->Write();
  treefile_->Close();

  //The following lines are displayed at the end of the run in root prompt
  cout<<"Total good events: "<<nEvtTotal<<endl;
  cout<<"Events with two reco electrons: "<<nrecopass<<endl;
  cout<<"Events with one reco electrons: "<<nonereco<<endl;
  cout<<"Events with gen dR<0.1:         "<<ngendRl0p1<<endl;
  cout<<"Events with one reco electron but two gsf tracks (dR_reco_track<0.1): "<<ntrackpass<<endl;
  cout<<"Events with two reco electron and gendR<0.1 : "<<nrecodRpass<<endl;
  cout<<"Events with 1 reco electron + two gsf tracks and gendR<0.1 : "<<ntrackdRpass<<endl;
  cout<<"Events with two reco electron and gendR>0.2 : "<<nrecodRpass2<<endl;
  cout<<"Events with 1 reco electron + two gsf tracks and gendR>0.2 : "<<ntrackdRpass2<<endl;
  cout<<"Total events with two gen electron: "<<nEvtpass<<endl;
  cout<<"Fraction of events with two reco electrons: "<<(float)nrecopass/(float)nEvtpass<<endl;
  cout<<"Fraction of events with one reco electron and two gsf track (dR_reco_track<0.1): "<<(float)ntrackpass/(float)nEvtpass<<endl;
  cout<<"Fraction of events with two reco electrons and gendR<0.1: "<<(float)nrecodRpass/(float)nEvtpass<<endl;
  cout<<"Fraction of events with 1 reco + 2 gsf tracks and gendR<0.1: "<<(float)ntrackdRpass/(float)nEvtpass<<endl;
  cout<<"Fraction of events with two reco electrons and gendR>0.2: "<<(float)nrecodRpass2/(float)nEvtpass<<endl;
  cout<<"Fraction of events with 1 reco + 2 gsf tracks and gendR>0.2: "<<(float)ntrackdRpass2/(float)nEvtpass<<endl;
  cout<<"------------------------------------------------------------------------------------------------------------"<<endl;
  cout<<"------------------------------------------------------------------------------------------------------------"<<endl;
  cout<<"Number of events with gendR<0.1:                               "<<ngen_dRl0p1<<endl;
  cout<<"Number of events with 0.1<gendR<0.2:                           "<<ngen_dRg0p2<<endl;
  cout<<"Number of events with 0 reco+gendR<0.1:                        "<<n0reco_dRl0p1<<endl;
  cout<<"Number of events with 1 reco+gendR<0.1:                        "<<n1reco_dRl0p1<<endl;
  cout<<"Number of events with 2 reco+gendR<0.1:                        "<<n2reco_dRl0p1<<endl;
  cout<<"Number of events with 0 reco+0.1<gendR<0.2:                    "<<n0reco_dRg0p2<<endl;
  cout<<"Number of events with 1 reco+0.1<gendR<0.2:                    "<<n1reco_dRg0p2<<endl;
  cout<<"Number of events with 2 reco+0.1<gendR<0.2:                    "<<n2reco_dRg0p2<<endl;

  
  //The following lines are written on the sum_<sample>.txt file
  ofstream fout(_SumFileName);
  fout<<"Total events ran = "<<nEvtTotal<<endl;
  fout<<"Total good events: "<<nEvtTotal<<endl;
  fout<<"Events with two reco electrons: "<<nrecopass<<endl;
  fout<<"Events with one reco electron but two gsf tracks (dR_reco_track<0.1): "<<ntrackpass<<endl;
  fout<<"Events with two reco electron and gendR<0.1 : "<<nrecodRpass<<endl;
  fout<<"Events with 1 reco electron + two gsf tracks and gendR<0.1 : "<<ntrackdRpass<<endl;
  fout<<"Events with two reco electron and gendR>0.2 : "<<nrecodRpass2<<endl;
  fout<<"Events with 1 reco electron + two gsf tracks and gendR>0.2 : "<<ntrackdRpass2<<endl;
  fout<<"Total events with two gen electron: "<<nEvtpass<<endl;
  fout<<"Fraction of events with two reco electrons: "<<(float)nrecopass/(float)nEvtpass<<endl;
  fout<<"Fraction of events with one reco electron and two gsf track (dR_reco_track<0.1): "<<(float)ntrackpass/(float)nEvtpass<<endl;
  fout<<"Fraction of events with two reco electrons and gendR<0.1: "<<(float)nrecodRpass/(float)nEvtpass<<endl;
  fout<<"Fraction of events with 1 reco + 2 gsf tracks and gendR<0.1: "<<(float)ntrackdRpass/(float)nEvtpass<<endl;
  fout<<"Fraction of events with two reco electrons and gendR>0.2: "<<(float)nrecodRpass2/(float)nEvtpass<<endl;
  fout<<"Fraction of events with 1 reco + 2 gsf tracks and gendR>0.2: "<<(float)ntrackdRpass2/(float)nEvtpass<<endl;
  fout<<"------------------------------------------------------------------------------------------------------------"<<endl;
  fout<<"------------------------------------------------------------------------------------------------------------"<<endl;
  fout<<"Number of events with gendR<0.1:                               "<<ngen_dRl0p1<<endl;
  fout<<"Number of events with 0.1<gendR<0.2:                           "<<ngen_dRg0p2<<endl;
  fout<<"Number of events with 0 reco+gendR<0.1:                        "<<n0reco_dRl0p1<<endl;
  fout<<"Number of events with 1 reco+gendR<0.1:                        "<<n1reco_dRl0p1<<endl;
  fout<<"Number of events with 2 reco+gendR<0.1:                        "<<n2reco_dRl0p1<<endl;
  fout<<"Number of events with 0 reco+0.1<gendR<0.2:                    "<<n0reco_dRg0p2<<endl;
  fout<<"Number of events with 1 reco+0.1<gendR<0.2:                    "<<n1reco_dRg0p2<<endl;
  fout<<"Number of events with 2 reco+0.1<gendR<0.2:                    "<<n2reco_dRg0p2<<endl;
  
}

void MergedElectrons_ana::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

}

void MergedElectrons_ana::BookHistograms()
{
  
  //Some initial histograms 
  h.nGsf = new TH1F("Number of Gsf electrons","Number of Gsf electrons",10,0,10);
  h.pT = new TH1F("pT_gsf","pT_gsf",500,0,500);
  h.dntrack_ngsftrack = new TH1F("ntrack_minus_ngsftrack","ntrack_minus_ngsftrack",10,0,10);
  h.dngsf_ngsftrack = new TH1F("ngsf_minus_ngsftrack","ngsf_minus_ngsftrack",10,0,10);
  h.ntrack_dR01 = new TH1F("ntrack_dR01_ofgsftrack","ntrack_dR01_ofgsftrack",10,0,10);
  h.ntrack_dR0p05 = new TH1F("ntrack_dR0p05_ofgsftrack","ntrack_dR0p05_ofgsftrack",10,0,10);
  
  //Misc Histograms defined here
  h.misc_histo[0] = new TH1F("sigmaEtaEta_EB","sigmaEtaEta_EB",100,0,0.1);
  h.misc_histo[1] = new TH1F("sigmaEtaEta_EE","sigmaEtaEta_EE",100,0,0.1);
  h.misc_histo[2] = new TH1F("ntrack_dR0p1 of gsf electron","ntrack_dR0p1 of gsf electron",10,0,10);
  h.misc_histo[3] = new TH1F("ngsftrack_dR0p1 of gsf seedtrack","ngsftrack_dR0p1 of gsf seed track",10,0,10);
  h.misc_histo[4] = new TH1F("Mass reso if dR<0.03","Mass reso if dR<0.03",100,-5,5);
  h.misc_histo[5] = new TH1F("Mass reso if 0.03<dR<0.06","Mass reso if 0.03<dR<0.06",100,-5,5);
  h.misc_histo[6] = new TH1F("Mass reso if 0.06<dR<0.1","Mass reso if 0.06<dR<0.1",100,-5,5);
  h.misc_histo[7] = new TH1F("ngsftrack_dR0p05 of gsf seedtrack","ngsftrack_dR0p05 of gsf seed track",10,0,10);
  h.misc_histo[8] = new TH1F("dRminRHN_gsfElec_genElec","dRminRHN_gsfElec_genElec",5000,0,5);
  h.misc_histo[9] = new TH1F("dRminW_gsfElec_genElec","dRminW_gsfElec_genElec",5000,0,5);
  
  TString varname[28] = {"pT","Eta","sigmaEtaEta","sigmaIetaIeta","sigmaIphiIphi","e1x5","e5x5","e2x5Max","r9","eSuperClusterOverP","dE_raw_overE_corr","Eraw_overP","Ecorr_overP","dr03TkSumPt","dr03TkSumPtHEEP","dr04TkSumPt","dr04TkSumPtHEEP","deltaEtaSuperClusterTrackAtVtx","deltaEtaSeedClusterTrackAtCalo","deltaPhiSuperClusterTrackatVtx","deltaPhiSeedClusterTrackatCalo","full5x5_sigmaEtaEta","full5x5_sigmaIetaIeta","full5x5_sigmaIphiIphi","full5x5_r9","full5x5_hcaloverEcal","hcaloverEcal","Eraw_over_gsftrackp"};
  
  float bin_min[28] = {0,-5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,-1,-1,0,0,0,0,0,0,0};
  float bin_max[28] = {250,5,0.1,0.1,0.1,3000,3000,3000,1.2,100,1,5,5,100,100,100,100,1,1,1,1,0.1,0.1,0.1,1.1,1.5,1,10};
  int nbin[28] = {250,100,100,100,100,3000,3000,3000,1200,10000,1000,1000,1000,1000,1000,1000,1000,2000,2000,2000,2000,100,100,100,1100,1500,1000,1000};
  int varsize = sizeof(varname)/sizeof(TString);
  for(int ivar=0;ivar<varsize;ivar++){
    TString name = varname[ivar];
    h.Gsfprop_histo[ivar] = new TH1F(name,name,nbin[ivar],bin_min[ivar],bin_max[ivar]);
  }
  
  //Track properties histogram 
  TString trackvarname[2] = {"normalizedChi2","gsftrackp"};
  h.trackprop_histo[0] = new TH1F(trackvarname[0],trackvarname[0],5000,0,50);
  
  //GenElectron histograms
  h.genElec_histo[0] = new TH1F("number of gen electron","number of genelectron",10,0,10);
  h.genElec_histo[1] = new TH1F("genMee","genMee",2000,3.08,3.1);
  h.genElec_histo[2] = new TH1F("genJpsi size","genJpsi size",10,0,10);
  h.genElec_histo[3] = new TH1F("gendRee","gendRee",5000,0,5);
  h.genElec_histo[4] = new TH1F("Alpha_dR0p1","Alpha_dR0p1",1100,0,1.1);
  h.genElec_histo[5] = new TH1F("beta_dR0p1","beta_dR0p1",1100,0,1.1);
  h.genElec_histo[6] = new TH1F("sum_Alpha_beta_dR0p1","sum_Alpha_beta_dR0p1",1100,0,1.1);
  h.genElec_histo[7] = new TH1F("alpha_pT","alpha_pT",1100,0,1.1);
  h.genElec_histo[8] = new TH1F("beta_pT","beta_pT",1100,0,1.1);
  h.genElec_histo[9] = new TH1F("pT_reso","pT_reso",2100,-1,1.1);
  h.genElec_histo[10] = new TH1F("rawE_o_genE","rawE_o_genE",200,0,2);
  h.genElec_histo[11] = new TH1F("corrE_o_genE","corrE_o_genE",200,0,2);
  h.genElec_histo[12] = new TH1F("genElec_dxy","genElec_dxy",10000,0,1);
  h.genElec_histo[13] = new TH1F("genElec_dz","genElec_dz",200,-10,10);
  h.genElec_histo[14] = new TH1F("genpTee","genpTee",600,0,600);
  h.genElec_histo[15] = new TH1F("genElec_mass","genElec_mass",20000,-1,1);
  h.genElec_histo[16] = new TH1F("genElec_momRHN_size","genElec_momRHN_size",10,0,10);
  h.genElec_histo[17] = new TH1F("genElec_momW_size","genElec_momW_size",10,0,10);
  h.genElec_histo[18] = new TH1F("genElec_mom_pdgid","genElec_momW_pdgid",1200,-600,600);
  h.genElec_histo[19] = new TH1F("dRee_momW_momRHN","dRee_momW_momRHN",1000,0,10);
  
  //Gsftrack histograms
  h.gsftrack_histo[0] = new TH1F("M_gsftrack_gsftrack","M_gsftrack_gsftrack",10000,0,100);
  h.gsftrack_histo[1] = new TH1F("dR_gsftrack_gsftrack","dR_gsftrack_gsftrack",1000,0,1);
  h.gsftrack_histo[2] = new TH1F("dR_all_gsftrack_gsftrack","dR_all_gsftrack_gsftrack",5000,0,5);
  h.gsftrack_histo[3] = new TH1F("pT_gsftrack_gsftrack","pT_gsftrack_gsftrack",600,0,600);
  h.gsftrack_histo[4] = new TH1F("pT_gsftrack_gsftrack_merged","pT_gsftrack_gsftrack_merged",600,0,600);
  h.gsftrack_histo[5] = new TH1F("pT_gtgtOvermergedpT","pT_gtgtOvermergedpT",10000,0,100);
  h.gsftrack_histo[6] = new TH1F("pTgtgt_pTgtgtOvermergedpT_0p9to1p1","pTgtgt_pTgtgtOvermergedpT_0p9to1p1",500,0,500);
  h.gsftrack_histo[7] = new TH1F("pTmergedgsf_pTgtgtOvermergedpT_0p9to1p1","pTmergedgsf_pTgtgtOvermergedpT_0p9to1p1",500,0,500);
  h.gsftrack_histo[8] = new TH1F("pTdiffmergedgsfgtgt_pTgtgtOvermergedpT_0p9to1p1","pTdiffmergedgsfgtgt_pTgtgtOvermergedpT_0p9to1p1",200,-100,100);
  h.gsftrack_histo[9] = new TH1F("pTdiffgt0gt1_pTgtgtOvermergedpT_0p9to1p1","pTdiffgt0gt1_pTgtgtOvermergedpT_0p9to1p1",1000,-500,500);
  h.gsftrack_histo[10] = new TH1F("pTassymetrygt0gt1_pTgtgtOvermergedpT_0p9to1p1","pTassymetrygt0gt1_pTgtgtOvermergedpT_0p9to1p1",1000,0,10);
  h.gsftrack_histo[11] = new TH1F("mergedpTdiffJPsigsftrack","mergedpTdiffJPsigsftrack",1000,-500,500);
  h.gsftrack_histo[12] = new TH1F("mergedpTdiffJPsigsftrackOverpTJPsi","mergedpTdiffJPsigsftrackOverpTJPsi",2000,-10,10);
  
  
  //Gsfelectron histograms
  h.gsf_histo[0] = new TH1F("dR_gsf_gsf","dR_gsf_gsf",10000,0,10);
  h.gsf_histo[1] = new TH1F("alpha_gsfenergy_over_rawenergy","alpha_gsfenergy_over_rawenergy",1000,0,10);
  h.gsf_histo[2] = new TH1F("GsfElec_size","GsfElec_size",10,0,10);
  h.gsf_histo[3] = new TH1F("reco_Mee","reco_Mee",10000,0,100);
  h.gsf_histo[4] = new TH1F("M_gsftrack_rawEnergy","M_gsftrack_rawEnergy",20000,-1000,1000);
  h.gsf_histo[5] = new TH1F("M_gsftrack_Energy","M_gsftrack_Energy",20000,-1000,1000);
  h.gsf_histo[6] = new TH1F("M2_gsftrack_rawEnergy","M2_gsftrack_rawEnergy",200000,-10000,10000);
  h.gsf_histo[7] = new TH1F("M2_gsftrack_Energy","M2_gsftrack_Energy",200000,-10000,10000);
  h.gsf_histo[8] = new TH1F("dRmin_genElec_gsfElec","dRmin_genElec_gsfelec",500,0,5);
  h.gsf_histo[9] = new TH1F("gsf_rawEnergy","gsf_rawEnergy",5000,0,500);
  h.gsf_histo[10] = new TH1F("gsf_correctedEnergy","gsf_correctedEnergy",5000,0,500);
  h.gsf_histo[11] = new TH1F("gsf_ErawOvertrackP","gsf_ErawOvertrackP",1000,0,10);
  h.gsf_histo[12] = new TH1F("gsf_EcorrOvertrackP","gsf_EcorrOvertrackP",1000,0,10);
  h.gsf_histo[13] = new TH1F("gsf_EcorrOverEraw","gsf_EcorrOverEraw",10000,0,10);
  h.gsf_histo[14] = new TH1F("gsf_deltaEcorrErawOverEraw","gsf_deltaEcorrErawOverEraw",1000,0,1);
  h.gsf_histo[15] = new TH1F("gsf_hcalOverEcal","gsf_hcalOverEcal",1000,0,10);
  h.gsf_histo[16] = new TH1F("gsf_hcalOverEcalBc","gsf_hcalOverEcalBc",1000,0,10);
  h.gsf_histo[17] = new TH1F("gsf_seedtrackPt","gsf_seedtrackPt",1000,0,1000);
  h.gsf_histo[18] = new TH1F("gsf_Pt","gsf_Pt",1000,0,1000);
  h.gsf_histo[19] = new TH1F("mergedgsf_seedtrackPt","mergedgsf_seedtrackPt",1000,0,1000);
  h.gsf_histo[20] = new TH1F("mergedgsf_Pt","mergedgsf_Pt",1000,0,1000);
  h.gsf_histo[21] = new TH1F("mergedgsf_Ptdiff_gsfgt","mergedgsf_Ptdiff_gsfgt",1000,-500,500);
  h.gsf_histo[22] = new TH1F("gsf_Ptdiff_gsfgt","gsf_Ptdiff_gsfgt",1000,-500,500);
  h.gsf_histo[23] = new TH1F("mergedgsf_PtdiffgsfgtOverPtgsf","mergedgsf_PtdiffgsfgtOverPtgsf",1000,-10,10);
  h.gsf_histo[24] = new TH1F("gsf_PtdiffgsfgtOverPtgsf","gsf_PtdiffgsfgtOverPtgsf",1000,-10,10);
  h.gsf_histo[25] = new TH1F("mergedgsf_ErawOvertrackP","mergedgsf_ErawOvertrackP",1000,0,10);
  h.gsf_histo[26] = new TH1F("mergedgsf_EcorrOvertrackP","mergedgsf_EcorrOvertrackP",1000,0,10);
  h.gsf_histo[27] = new TH1F("mergedgsf_pTdiffgsfgtOverPtgsf","mergedgsf_pTdiffgsfgtOverPtgsf",1000,-10,10);
  h.gsf_histo[28] = new TH1F("gsf_M_Eraw_seedgsftrk","gsf_M_Eraw_seedgsftrk",1100000,-10,100);
  h.gsf_histo[29] = new TH1F("gsf_M_Ecorr_seedgsftrk","gsf_M_Ecorr_seedgsftrk",1100000,-10,100);
  h.gsf_histo[30] = new TH1F("gsf_M_Eraw_gsf","gsf_M_Eraw_gsf",1100000,-10,100);
  h.gsf_histo[31] = new TH1F("gsf_M_Ecorr_gsf","gsf_M_Ecorr_gsf",1100000,-10,100);
  h.gsf_histo[32] = new TH1F("mergedgsf_M_Eraw_seedgsftrk","gsf_M_Eraw_seedgsftrk",1100000,-10,100);
  h.gsf_histo[33] = new TH1F("mergedgsf_M_Ecorr_seedgsftrk","gsf_M_Ecorr_seedgsftrk",1100000,-10,100);
  h.gsf_histo[34] = new TH1F("mergedgsf_M_Eraw_gsf","gsf_M_Eraw_gsf",1100000,-10,100);
  h.gsf_histo[35] = new TH1F("mergedgsf_M_Ecorr_gsf","gsf_M_Ecorr_gsf",1100000,-10,100);
  h.gsf_histo[36] = new TH1F("genmatched_gsf_M_Eraw_seedgsftrk","genmatched_gsf_M_Eraw_seedgsftrk",1100000,-10,100);
  h.gsf_histo[37] = new TH1F("genmatched_gsf_M_Ecorr_seedgsftrk","genmatched_gsf_M_Ecorr_seedgsftrk",1100000,-10,100);
  h.gsf_histo[38] = new TH1F("genmatched_gsf_M_Eraw_gsf","genmatched_gsf_M_Eraw_gsf",1100000,-10,100);
  h.gsf_histo[39] = new TH1F("genmatched_gsf_M_Ecorr_gsf","genmatched_gsf_M_Ecorr_gsf",1100000,-10,100);
  h.gsf_histo[40] = new TH1F("genmatched_mergedgsf_M_Eraw_seedgsftrk","genmatched_gsf_M_Eraw_seedgsftrk",1100000,-10,100);
  h.gsf_histo[41] = new TH1F("genmatched_mergedgsf_M_Ecorr_seedgsftrk","genmatched_gsf_M_Ecorr_seedgsftrk",1100000,-10,100);
  h.gsf_histo[42] = new TH1F("genmatched_mergedgsf_M_Eraw_gsf","genmatched_gsf_M_Eraw_gsf",1100000,-10,100);
  h.gsf_histo[43] = new TH1F("genmatched_mergedgsf_M_Ecorr_gsf","genmatched_gsf_M_Ecorr_gsf",1100000,-10,100);
  h.gsf_histo[44] = new TH1F("gsf_deltaEtaSCtrkatVtx","gsf_deltaEtaSCtrkatVtx",2000,-0.1,0.1);
  h.gsf_histo[45] = new TH1F("gsf_deltaPhiSCtrkatVtx","gsf_deltaPhiSCtrkatVtx",2000,-0.1,0.1);
  h.gsf_histo[46] = new TH1F("gsf_r9","gsf_r9",1200,0,1.2);
  h.gsf_histo[47] = new TH1F("gsfElec_momRHN_size","gsfElec_momRHN_size",10,0,10);
  h.gsf_histo[48] = new TH1F("gsfElec_momW_size","gsfElec_momW_size",10,0,10);


  //Merged gsf electron Histograms
  h.mergedgsf_histo[0] = new TH1F("mergedgsf_r9","mergedgsf_r9",1200,0,1.2);
  h.mergedgsf_histo[1] = new TH1F("mergedgsf_sigmaEtaEta","mergedgsf_sigmaEtaEta",100,0,0.1);
  h.mergedgsf_histo[2] = new TH1F("mergedgsf_sigmaPhiPhi","mergedgsf_sigmaPhiPhi",100,0,0.1);
  h.mergedgsf_histo[3] = new TH1F("mergedgsf_sigmaIetaIeta","mergedgsf_sigmaIetaIeta",100,0,0.1);
  h.mergedgsf_histo[4] = new TH1F("mergedgsf_sigmaIphiIphi","mergedgsf_sigmaIphiIphi",100,0,0.1);
  h.mergedgsf_histo[5] = new TH1F("mergedgsf_EcorrOverEraw","mergedgsf_EcorrOverEraw",10000,0,10);
  h.mergedgsf_histo[6] = new TH1F("mergedgsf_deltaEtaSCtrkatVtx","mergedgsf_deltaEtaSCtrkatVtx",2000,-0.1,0.1);
  h.mergedgsf_histo[7] = new TH1F("mergedgsf_deltaPhiSCtrkatVtx","mergedgsf_deltaPhiSCtrkatVtx",2000,-0.1,0.1);
  h.mergedgsf_histo[8] = new TH1F("mergedgsf_ErawOvertrkP","mergedgsf_ErawOvertrkP",1000,0,10);
  h.mergedgsf_histo[9] = new TH1F("mergedgsf_EcorrOvertrkP","mergedgsf_EcorrOvertrkP",1000,0,10);

  //GenParticle Histograms
  h.genPart_histo[0] = new TH1F("pdgid_genparticle_g1000","pdgid_genparticle",10000,1000,11000);
  h.genPart_histo[1] = new TH1F("dRmin_genElec_genPart","dRmin_genElec_genPart",500,0,5);
  h.genPart_histo[2] = new TH1F("nhadrons_dR0p1_gsfelectron","nhadrons_dR0p1_gsfelectron",20,0,20);
  h.genPart_histo[3] = new TH1F("nphotons_dR0p1_gsfelectron","nphotons_dR0p1_gsfelectron",20,0,20);

  //Gen RHN Histograms

  h.genRHN_histo[0] = new TH1F("size_genRHN","size_genRHN",10,0,10);

   
  //2D Histograms
  h.gen_histo[0] = new TH2F("dRee_vs_pTJPsi","dRee_vs_pTJPsi",500,0,5,1000,0,1000);

  
}

//Sorting the relevant arrays by their pT

void MergedElectrons_ana::sort(int opt)
{
  if(opt==0){
    
    for(int i=0;i<(int)gsfElec.size()-1;i++){
      for(int j=i+1;j<(int)gsfElec.size();j++){
	if(gsfElec[i].v.Pt() < gsfElec[j].v.Pt()) swap(gsfElec.at(i),gsfElec.at(j));
      }
    }
    
    for(int i=0;i<(int)genElec.size()-1;i++){
      for(int j=i+1;j<(int)genElec.size();j++){
	if(genElec[i].v.Pt() < genElec[j].v.Pt()) swap(genElec.at(i),genElec.at(j));
      }
    }
    
    for(int i=0;i<(int)gsftrack.size()-1;i++){
      for(int j=i+1;j<(int)gsftrack.size();j++){
	if(gsftrack[i].v.Pt() < gsftrack[j].v.Pt()) swap(gsftrack.at(i),gsftrack.at(j));
      }
    }
    
    for(int i=0;i<(int)gsftrack_dR0p1.size()-1;i++){
      for(int j=i+1;j<(int)gsftrack_dR0p1.size();j++){
	if(gsftrack_dR0p1[i].v.Pt() < gsftrack_dR0p1[j].v.Pt()) swap(gsftrack_dR0p1.at(i),gsftrack_dR0p1.at(j));
      }
    }
    
    
  }
}
