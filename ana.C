#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
/*
This is a driver script.
It decides which code to run over which sample, the names
of output files and so on.
*/


void ana(int sample=0){
  const char *hstfilename, *sumfilename, *myfilename, *treefilename;
  //Declare a chain for input files.
  TChain *chain = new TChain("demo/tree"); 
  //Declare an instance of our code class
  MergedElectrons_ana m_selec;
  
  if(sample==0){
    //Add one file to chain. This is the input file.
    chain->Add("/home/soumya/CMSSW_files/NTuples/root_files/JPsitoeepT40to500_ntuple.root");
    //Set Names of outputfiles
    hstfilename = "hst_JPsi.root";
    sumfilename = "sum_JPsi.txt";
    treefilename = "myTree_rootfiles/myTree_JPsi.root";
    //Set some options
    m_selec.SetData(0); //MC=0, data=1
    //m_selec.SetYear(2016);
  }
  if(sample==1){
    chain->Add("/home/soumya/CMSSW_files/NTuples/root_files/ZtoeepT40to300_ntuple.root");
    hstfilename = "hst_Z.root";
    sumfilename = "sum_Z.txt";
    treefilename = "myTree_rootfiles/myTree_Z.root";
    m_selec.SetData(0); //MC=0, data=1
    // m_selec.SetYear(2016);
  }
  if(sample==2){
    chain->Add("/home/soumya/CMSSW_files/NTuples/root_files/QCDpT150to200_ntuple.root");
    hstfilename = "hst_QCD.root";
    sumfilename = "sum_QCD.txt";
    treefilename = "myTree_rootfiles/myTree_QCD.root";
    m_selec.SetData(0); //MC=0, data=1
    // m_selec.SetYear(2016);
  }
  
  if(sample==3){
    chain->Add("/home/soumya/CMSSW_files/NTuples/root_files/EGamma_Run2018A_UL2018_v1_data_ntuple_6.root");
    hstfilename = "hst_EGamma_2018A_UL_data.root";
    sumfilename = "sum_EGamma_2018A_UL_data.txt";
    treefilename = "myTree_rootfiles/myTree_EGamma_2018A.root";
    m_selec.SetData(1); //MC=0, data=1
    // m_selec.SetYear(2016);
  }
  
  if(sample==4){
    chain->Add("/home/soumya/CMSSW_files/NTuples/root_files/HNL_trilepton_M2_V0p01_e_ntuple.root");
    hstfilename = "hst_HNL_trilepton_M2_V0p01_e.root";
    sumfilename = "sum_HNL_trilepton_M2_V0p01_e.txt";
    treefilename = "myTree_rootfiles/myTree_HNL_trilepton_M2_V0p01_e.root";
    m_selec.SetData(0); //MC=0, data=1
    // m_selec.SetYear(2016);
  }

  //Kept only if you want to have to execute a data file
  if(sample==5){
    chain->Add("inputs/SingleElectron_2016FpostVFP_1.root");
    hstfilename = "hst_SingleElectron_RunII_2016postVFP.root";
    sumfilename = "sum_SingleElectron_RunII_2016postVFP.txt";
    //m_selec.SetData(1); //MC=0, data=1
    // m_selec.SetYear(2016);
  }
  
   std::cout<<"Output files are "<<hstfilename<<" and "<<sumfilename<<std::endl;
   // Set some more options.. set the output file names.
   m_selec.SetHstFileName(hstfilename);
   m_selec.SetSumFileName(sumfilename);
   m_selec.SetTreeFileName(treefilename);
   // m_selec.SetVerbose(10);//set verbosity level for output.
   // Call the process function which runs the code.
   chain->Process(&m_selec);
   
}
