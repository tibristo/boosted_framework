#include "plotter.hpp"
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <utility> 
using namespace std;

// Read in the list of algorithms that we are going to run over.  This is obtained from a file called plotter.config
std::vector<std::string> getAlgs()
{
  ifstream f("plotter.config");
  std::string line;
  std::vector<std::string> algs;
  while (getline(f,line))
    {
      algs.push_back(rtrim(line));
    }
  f.close();
  return algs;
} // getAlgs


// Main method reads in the algorithms we're going to run over using getAlgs(), we can then set which ptbins we want to run over, but ideally this should also be read in from a config file.  An important thing
// to note here is that the program assumes that the names of the input files are in the form ALGORITHM_ptbin_sig(or bkg).root.  It gets the physics tree from these files and creates a new plotter object and sets it running.
// It also reads in the pt reweighting histograms from each file and uses this in the plotter constructor
int main(int argc, char * argv[])
{
  gROOT->ProcessLine("#include <vector>"); // With certain versions of ROOT this line is necessary

  std::string working_dir_suf = "/";
  if (argc > 1)
    working_dir_suf = std::string(argv[1])+"/";
  std::vector<std::string> algs = getAlgs();
  std::map<std::string, std::string> algIdx = {{"TopoSplitFilteredMu67SmallR0YCut9","1"},{"TopoSplitFilteredMu100SmallR30YCut4","2"},{"TopoTrimmedPtFrac5SmallR30","3"},{"TopoTrimmedPtFrac5SmallR20","4"},{"TopoPrunedCaRcutFactor50Zcut10","5"},{"TopoPrunedCaRcutFactor50Zcut20","6"},{"AntiKt2LCTopo","7"},{"AntiKt3LCTopo","8"},{"AntiKt4LCTopo","9"}};
  std::vector<std::string> ptbins = {"inclusive"};//,"gtr1000","500to1000","1000to1500","1500to2000","gtr2000"};

  std::cout << "in main" << std::endl;
  for (std::vector<std::string>::iterator it = algs.begin(); it != algs.end(); it++)
    {
      std::cout << (*it) << std::endl;
      for (std::vector<std::string>::iterator it2 = ptbins.begin(); it2!= ptbins.end(); it2++)
	{
	  // for filename: algname_algIndex_ptbins_sig/bkg.root
	  std::string working_dir = (*it)+working_dir_suf;
	  std::string sig_fname = working_dir + (*it)+"_"+algIdx[(*it)]+"_"+(*it2)+"_sig.root";
	  TFile * f_sig = new TFile(sig_fname.c_str(),"READ");
	  TH1F * Wp_PtReweight = (TH1F*)f_sig->Get(std::string("Wp_ptreweight_"+(*it)).c_str());
	  //Wp_PtReweight->Scale(1./Wp_PtReweight->Integral());
	  std::cout << Wp_PtReweight << std::endl;
	  TTree * sig = (TTree *)f_sig->Get("physics");
	  std::string bkg_fname = working_dir+(*it)+"_"+algIdx[(*it)]+"_"+(*it2)+"_bkg.root";
	  TFile * f_bkg = new TFile(bkg_fname.c_str(),"READ");
	  TH1F * qcd_PtReweight = (TH1F*)f_bkg->Get(std::string("qcd_ptreweight_"+(*it)).c_str());
	  //qcd_PtReweight->Scale(1./qcd_PtReweight->Integral());
	  TTree * bkg = (TTree *)f_bkg->Get("physics");
	  plotter * p;
	  std::string name_suffix = algIdx[(*it)]+"_"+(*it2);
	  bool applyPt = true, applyEta = false;
	  p = new plotter((*it), bkg, sig, name_suffix, working_dir, qcd_PtReweight, Wp_PtReweight, applyPt, applyEta);
	  p->runPlotter();
	  delete sig;
	  delete bkg;
	  delete f_sig;
	  delete f_bkg;
	  delete p;
	}
    }
} // main


plotter::plotter(std::string & algorithm_in, TTree * bkg_in, TTree * sig_in, std::string & params, std::string & dir, TH1F * bkg_hist, TH1F * sig_hist, bool applyPt, bool applyEta)
  {
    /*
      The constructor for plotter.  It sets the signal and bkg TTrees, gets the variables we are going to run over, reads in the config file for the grooming algorithm we are using, clones the truth pt hists and 
      sets up the pt reweighting histogram (bkg / signal).
      
    // Keyword arguments:
    algorithm_in --- the name of the algorithm we are running on
    bkg_in --- background TTree
    sig_in --- signal TTree
    params --- The parameters here refer to the pt bin we are running over - inclusive, gtr2000, etc.
    dir --- the working directory
    bkg_hist --- The background TH1F containing the truth pt from the full sample (before mass window cuts applied) for reweighting the signal
    sig_hist --- The signal TH1F containing the truth pt from the full sample before mass window cuts are applied
    applyPtCut --- Whether or not to apply the pt window cut, 200 GeV < pt < 350 GeV
    applyEta --- whether or not to apply the eta cut of fabs(eta) > 1.2.  

    Note:  the pt cut and eta cut will eventually be changed to accept the pt and eta ranges, I'll keep you posted!
    */

    readAlgorithmConfig(algorithm_in); // reads in the variables and the parameters for the TH1Fs for each one - bin numbers, limits, etc.  This is all read into a map called config
    getWeights(algorithm_in, params); // get all of the weights for the different samples - k-factors, xs, efficiency
    sig = sig_in;
    bkg = bkg_in;
    algorithm = algorithm_in;
    algorithm_params = params;
    working_dir = dir;
    Wp_PtReweight = (TH1F*)sig_hist->Clone("wp_reweight");
    qcd_PtReweight = (TH1F*)bkg_hist->Clone("qcd_reweight");
    applyPtCut = applyPt;
    applyEtaCut = applyEta;
    pt_reweight_hist = new TH1F ("Pt Reweighting", "Pt Reweighting", 20, 0, 3500 ); // kinda random, but legacy...... using the same parameters as the input pt histograms.
    // set up the pt reweighting histogram.
    for (int j = 1; j <= pt_reweight_hist->GetNbinsX() ; j++)
      pt_reweight_hist->SetBinContent(j,qcd_PtReweight->GetBinContent(j)/Wp_PtReweight->GetBinContent(j));
  } // plotter()



void plotter::getVariables()
  {
    // loop through config (from readAlgorithmConfig() ) and get the variable names and set up a map, called vars.
    // this map just stores the name of the variable, so E, m , pt and then the name this has in the ntuple, so jet_CamKt12Topo_E, for example

    for (std::map<std::string, std::vector<std::string> >::iterator it = config.begin(); it != config.end() ; it++)
      {
	vars[(*it).first] = (*it).second[0];
      }
  
  } // getVariables()


void plotter::runPlotter()
{
  // This is really the driver method for the plotter.  After it has been constructed, all of the relevant methods are run to create the histograms and fill the histograms,
  // produce the nice looking plots and the ROC curves.

  // set up all of the histograms.  We have a vector of TH1Fs, so first we call initHists to push_back(0) for each one, then we create each Histogram.
  initHists();
  createHistograms(config);

  // setBranchAddresses for all of the variables
  setBranches();

  // fill the signal and background histograms
  fillHistograms(sig, true); // true is signal
  fillHistograms(bkg, false); // false is bkg

  plotHistograms(); // plot bkg, signal and combined
  produceAllROCs(std::string(algorithm+algorithm_params)); // algorithm+ptbin
  
  clearHists();
} // runPlotter()




void plotter::produceAllROCs(std::string algorithm)
{
  /* 
     produce the ROC for each variable we're interested in
     
     Keyword Args:
     algorithm --- Name of the algorithm+ptbin
  */

  int type = 1; // this sets the type of roc curve we want to create
  bool draw = true;

  int vec_idx = jet_type::GROOMED;  // only interested in the groomed ones at the moment, maybe this will change?

  TCanvas * c = new TCanvas();
  //bkg_m[2]->Scale(1/bkg_m[2]->Integral());
  //signal_m[2]->Scale(1/signal_m[2]->Integral());
  //bkg_pt[2]->Draw();
  //signal_pt[2]->Draw("same");
  signal_PlanarFlow[2]->SetLineColor(kRed);
  signal_PlanarFlow[2]->Draw();
  c->SaveAs("TEST_PFSIG.png");
  bkg_PlanarFlow[2]->Draw();
  c->SaveAs("TEST_PFBKG.png");
  delete c;
  //TH1F * sigTau2Tau1 = ;
  //TH1F * bkgTau2Tau1 = ;
  compareHistograms(signal_Tau1[vec_idx],bkg_Tau1[vec_idx],"Test_Tau1");
  //compareHistograms(sigTau2Tau1[vec_idx],bkgTau2Tau1[vec_idx],"Test_Tau2Tau1");
  compareHistograms(signal_SPLIT12[vec_idx],bkg_SPLIT12[vec_idx],"Test_SPLIT12");
  compareHistograms(signal_PlanarFlow[vec_idx],bkg_PlanarFlow[vec_idx],"Test_PlanarFlow");
  compareHistograms(signal_m[vec_idx],bkg_m[vec_idx],"Test_m");
  
  
  std::string var_prefix = "groomed_";

  //if (config.find(var_prefix+"E")!=config.end())
  produceROC(type, signal_E[vec_idx], bkg_E[vec_idx], algorithm+"_E", draw);
  produceROC(type, signal_pt[vec_idx], bkg_pt[vec_idx], algorithm+"_pt", draw);
  produceROC(type, signal_m[vec_idx], bkg_m[vec_idx], algorithm+"_m", draw);
  produceROC(type, signal_eta[vec_idx], bkg_eta[vec_idx], algorithm+"_eta", draw);
  produceROC(type, signal_phi[vec_idx], bkg_phi[vec_idx], algorithm+"_phi", draw);
  ///// BECAREFUL WITH THIS ONE!!! doesn't exist for truth..... [0]
  produceROC(type, signal_emfrac[vec_idx], bkg_emfrac[vec_idx], algorithm+"_emfrac", draw);

  produceROC(type, signal_Tau1[vec_idx], bkg_Tau1[vec_idx], algorithm+"_Tau1", draw);
  produceROC(type, signal_Tau2[vec_idx], bkg_Tau2[vec_idx], algorithm+"_Tau2", draw);
  produceROC(type, signal_Tau3[vec_idx], bkg_Tau3[vec_idx], algorithm+"_Tau3", draw);
  produceROC(type, signal_WIDTH[vec_idx], bkg_WIDTH[vec_idx], algorithm+"_WIDTH", draw);
  produceROC(type, signal_SPLIT12[vec_idx], bkg_SPLIT12[vec_idx], algorithm+"_SPLIT12", draw);
  produceROC(type, signal_SPLIT23[vec_idx], bkg_SPLIT23[vec_idx], algorithm+"_SPLIT34", draw);
  produceROC(type, signal_SPLIT34[vec_idx], bkg_SPLIT34[vec_idx], algorithm+"_SPLIT34", draw);
  produceROC(type, signal_Dip12[vec_idx], bkg_Dip12[vec_idx], algorithm+"_Dip12", draw);
  produceROC(type, signal_Dip13[vec_idx], bkg_Dip13[vec_idx], algorithm+"_Dip13", draw);
  produceROC(type, signal_Dip23[vec_idx], bkg_Dip23[vec_idx], algorithm+"_Dip23", draw);
  produceROC(type, signal_DipExcl12[vec_idx], bkg_DipExcl12[vec_idx], algorithm+"_DipExcl12", draw);
  produceROC(type, signal_PlanarFlow[vec_idx], bkg_PlanarFlow[vec_idx], algorithm+"_PlanarFlow", draw);
  produceROC(type, signal_Angularity[vec_idx], bkg_Angularity[vec_idx], algorithm+"_Angularity", draw);
  produceROC(type, signal_QW[vec_idx], bkg_QW[vec_idx], algorithm+"_QW", draw);
  produceROC(type, signal_PullMag[vec_idx], bkg_PullMag[vec_idx], algorithm+"_PullMag", draw);
  produceROC(type, signal_PullPhi[vec_idx], bkg_PullPhi[vec_idx], algorithm+"_PullPhi", draw);
  produceROC(type, signal_Pull_C00[vec_idx], bkg_Pull_C00[vec_idx], algorithm+"_Pull_C00", draw);
  produceROC(type, signal_Pull_C01[vec_idx], bkg_Pull_C01[vec_idx], algorithm+"_Pull_C01", draw);
  produceROC(type, signal_Pull_C10[vec_idx], bkg_Pull_C10[vec_idx], algorithm+"_Pull_C10", draw);
  produceROC(type, signal_Pull_C11[vec_idx], bkg_Pull_C11[vec_idx], algorithm+"_Pull_C11", draw);
  
} //produceAllROC()


TGraph plotter::produceROC(int type, TH1F *&S, TH1F *&B, TString name, bool draw)//, TGraph &curve_up, TGraph &curve_do){
  //Usage:
  //TH1D *S = new TH1D();
  //TH1D *B = new TH1D();
  //...Fill S and B with signal (S) and background (B)
  //Need S and B to have the same number of bins!
  //TGraph *curve = new TGraph(M); where M is however fine you want the ROC curve binning to be.
  //ROC(S,B,curve);

// This will output a TGraph to name.png with the ROC if draw is set to true.  Ideally we want this to return the TGraph as well, which I will do soon!
{
  const int n = S->GetNbinsX();
  vector<double> s;
  vector<double> b;
  vector<double> r;
  vector<double> serr;
  vector<double> berr;

  //cout<<"Getting init ROC"<<endl;
  for (int i=0;i<n;i++){
    s.push_back(S->GetBinContent(i+1));
    b.push_back(B->GetBinContent(i+1));
    serr.push_back(S->GetBinError(i+1));
    berr.push_back(B->GetBinError(i+1));
    if (B->GetBinContent(i+1)>0){
      r.push_back(S->GetBinContent(i+1)/B->GetBinContent(i+1));
    }
    else{
      r.push_back(-1);
    }
  }

  //sort by ascending order
  float temp_s=1;
  float temp_b=1;
  float temp_r=1;
  float temp_serr=1;
  float temp_berr=1;
  int sizes = s.size();
  for(int isort=sizes; isort>1; isort=isort-1){
    for(int i=0; i<isort-1; i++){
      if( r.at(i)<r.at(i+1) ){
	temp_s  = s.at(i);
	temp_b  = b.at(i);
	temp_serr  = serr.at(i);
	temp_berr  = berr.at(i);
	temp_r  = r.at(i);

	s.at(i) = s.at(i+1);
	b.at(i) = b.at(i+1);
	serr.at(i) = serr.at(i+1);
	berr.at(i) = berr.at(i+1);
	r.at(i) = r.at(i+1);

	s.at(i+1) = temp_s;
	b.at(i+1) = temp_b;
	serr.at(i+1) = temp_serr;
	berr.at(i+1) = temp_berr;
	r.at(i+1) = temp_r;
      }
    }
  }

  double totalB = B->Integral();
  double totalS = S->Integral();
  //put into graph
  //cout<<"NPoints: "<<n<<endl;
  TGraph gr(n);
  TGraph gr_up(n);
  TGraph gr_do(n);
  ofstream log;
  log.open(name+".txt");
  for (int i=0; i<n; i++){
    double myS = 0.;
    double myB = 0.;
    double mySerr = 0.;
    double myBerr = 0.;
    for (int j=0; j<i; j++){
      myS += s.at(j)/totalS;
      myB += b.at(j)/totalB;

      mySerr += pow(pow(mySerr, 2.0) + pow((serr.at(j)/totalS), 2.0), 0.5);
      myBerr += pow(pow(myBerr, 2.0) + pow((berr.at(j)/totalB), 2.0), 0.5);

    }
    if(type==1){
      log << myS << " , " << myB << endl;
      gr.SetPoint(i, myS, (1-myB));
      gr_up.SetPoint(i, myS+mySerr, (1-(myB-myBerr)));
      gr_do.SetPoint(i, myS-mySerr, (1-(myB+myBerr)));
    }
    if(type==2){
      if(myB==1)
	gr.SetPoint(i, myS, 100000);
      else{
	gr.SetPoint(i, myS, 1/(1-myB));
      }
    }
  }
  log.close();
  //curve=gr;
  if (draw==true)
    {
      TCanvas *c1 = new TCanvas();
      gr.Draw("");
      c1->SaveAs(working_dir+name+".png");
      delete c1;
    }
  //curve_up=gr_up;
  //curve_do=gr_do;
  return gr;
  
} //produceROC



void plotter::createHistograms(std::map<std::string, std::vector<std::string> > & config)
{
  /*
    Instantiate all of the histograms we are using.  At the moment this is a little stupid in that it creates _all_ of the variables, and we don't actually save any memory by not instantiating those histograms which we will not fill.
    Each histogram is stored in a vector.  This vector contains an entry for the truth, topo and groomed version of a certain variable.
    
    Keyword arguments:
    config --- This is the config file which is created by readAlgoConfig().  This is a global variable so we don't actually need to pass this in here.  However, the original code was structured differently and I haven't yet changed it.
    --- The config map contains the name of the variable in the ntuple and the parameters for the histogram - number of bins and axis limits.
  */
  
  for (int jet_idx = 0; jet_idx < 3 ; jet_idx++) // for truth, topo and groomed
    {
      std::string var_prefix = "";
      switch (jet_idx)
	{
	case jet_type::TRUTH:
	  var_prefix = "truth_";
	  break;
	case jet_type::TOPO:
	  var_prefix = "topo_";
	  break;
	case jet_type::GROOMED:
	  var_prefix = "groomed_";
	  break;
	}
      if (config.find(var_prefix+"E")!=config.end())
	{
	  signal_E[jet_idx] = new TH1F(std::string(config[var_prefix+"E"][0]+"_WP").c_str(), std::string(config[var_prefix+"E"][0]+"_WP").c_str(), std::stoi(config[var_prefix+"E"][1]), std::stod(config[var_prefix+"E"][2]) ,std::stod(config[var_prefix+"E"][3]));
	  bkg_E[jet_idx] = new TH1F(std::string(config[var_prefix+"E"][0]+"_QCD").c_str(), std::string(config[var_prefix+"E"][0]+"_QCD").c_str(), std::stoi(config[var_prefix+"E"][1]), std::stod(config[var_prefix+"E"][2]) ,std::stod(config[var_prefix+"E"][3]));
	}
	
      if (config.find(var_prefix+"pt")!=config.end())
	{
	  signal_pt[jet_idx] = new TH1F(std::string(config[var_prefix+"pt"][0]+"_WP").c_str(), std::string(config[var_prefix+"pt"][0]+"_WP").c_str(), std::stoi(config[var_prefix+"pt"][1]), std::stod(config[var_prefix+"pt"][2]) ,std::stod(config[var_prefix+"pt"][3]));
	  bkg_pt[jet_idx] = new TH1F(std::string(config[var_prefix+"pt"][0]+"_QCD").c_str(), std::string(config[var_prefix+"pt"][0]+"_QCD").c_str(), std::stoi(config[var_prefix+"pt"][1]), std::stod(config[var_prefix+"pt"][2]) , std::stod(config[var_prefix+"pt"][3]));
	}
      
      if (config.find(var_prefix+"m")!=config.end())
	{
	  signal_m[jet_idx] = new TH1F(std::string(config[var_prefix+"m"][0]+"_WP").c_str(), std::string(config[var_prefix+"m"][0]+"_WP").c_str(), std::stoi(config[var_prefix+"m"][1]), std::stod(config[var_prefix+"m"][2]) ,std::stod(config[var_prefix+"m"][3]));
	  bkg_m[jet_idx] = new TH1F(std::string(config[var_prefix+"m"][0]+"_QCD").c_str(), std::string(config[var_prefix+"m"][0]+"_QCD").c_str(), std::stoi(config[var_prefix+"m"][1]), std::stod(config[var_prefix+"m"][2]) ,std::stod(config[var_prefix+"m"][3]));
	}
      
      if (config.find(var_prefix+"eta")!=config.end())
	{
	  signal_eta[jet_idx] = new TH1F(std::string(config[var_prefix+"eta"][0]+"_WP").c_str(), std::string(config[var_prefix+"eta"][0]+"_WP").c_str(), std::stoi(config[var_prefix+"eta"][1]), std::stod(config[var_prefix+"eta"][2]) ,std::stod(config[var_prefix+"eta"][3]));
	  bkg_eta[jet_idx] = new TH1F(std::string(config[var_prefix+"eta"][0]+"_QCD").c_str(), std::string(config[var_prefix+"eta"][0]+"_QCD").c_str(), std::stoi(config[var_prefix+"eta"][1]), std::stod(config[var_prefix+"eta"][2]) ,std::stod(config[var_prefix+"eta"][3]));
	}
      
      if (config.find(var_prefix+"phi")!=config.end())
	{
	  signal_phi[jet_idx] = new TH1F(std::string(config[var_prefix+"phi"][0]+"_WP").c_str(), std::string(config[var_prefix+"phi"][0]+"_WP").c_str(), std::stoi(config[var_prefix+"phi"][1]), std::stod(config[var_prefix+"phi"][2]) ,std::stod(config[var_prefix+"phi"][3]));
	  bkg_phi[jet_idx] = new TH1F(std::string(config[var_prefix+"phi"][0]+"_QCD").c_str(), std::string(config[var_prefix+"phi"][0]+"_QCD").c_str(), std::stoi(config[var_prefix+"phi"][1]), std::stod(config[var_prefix+"phi"][2]) ,std::stod(config[var_prefix+"phi"][3]));
	}
      
      if (config.find(var_prefix+"emfrac")!=config.end() && jet_idx != jet_type::TRUTH) // no emfrac for truth 
	{
	  signal_emfrac[jet_idx] = new TH1F(std::string(config[var_prefix+"emfrac"][0]+"_WP").c_str(), std::string(config[var_prefix+"emfrac"][0]+"_WP").c_str(), std::stoi(config[var_prefix+"emfrac"][1]), std::stod(config[var_prefix+"emfrac"][2]) ,std::stod(config[var_prefix+"emfrac"][3]));
	  bkg_emfrac[jet_idx] = new TH1F(std::string(config[var_prefix+"emfrac"][0]+"_QCD").c_str(), std::string(config[var_prefix+"emfrac"][0]+"_QCD").c_str(), std::stoi(config[var_prefix+"emfrac"][1]), std::stod(config[var_prefix+"emfrac"][2]) ,std::stod(config[var_prefix+"emfrac"][3]));	  
	}
      
      if (config.find(var_prefix+"Tau1")!=config.end())
	{
	  signal_Tau1[jet_idx] = new TH1F(std::string(config[var_prefix+"Tau1"][0]+"_WP").c_str(), std::string(config[var_prefix+"Tau1"][0]+"_WP").c_str(), std::stoi(config[var_prefix+"Tau1"][1]), std::stod(config[var_prefix+"Tau1"][2]) ,std::stod(config[var_prefix+"Tau1"][3]));
	  bkg_Tau1[jet_idx] = new TH1F(std::string(config[var_prefix+"Tau1"][0]+"_QCD").c_str(), std::string(config[var_prefix+"Tau1"][0]+"_QCD").c_str(), std::stoi(config[var_prefix+"Tau1"][1]), std::stod(config[var_prefix+"Tau1"][2]) ,std::stod(config[var_prefix+"Tau1"][3]));
	}
      
      if (config.find(var_prefix+"Tau2")!=config.end())
	{
	  signal_Tau2[jet_idx] = new TH1F(std::string(config[var_prefix+"Tau2"][0]+"_WP").c_str(), std::string(config[var_prefix+"Tau2"][0]+"_WP").c_str(), std::stoi(config[var_prefix+"Tau2"][1]), std::stod(config[var_prefix+"Tau2"][2]) ,std::stod(config[var_prefix+"Tau2"][3]));
	  bkg_Tau2[jet_idx] = new TH1F(std::string(config[var_prefix+"Tau2"][0]+"_QCD").c_str(), std::string(config[var_prefix+"Tau2"][0]+"_QCD").c_str(), std::stoi(config[var_prefix+"Tau2"][1]), std::stod(config[var_prefix+"Tau2"][2]) ,std::stod(config[var_prefix+"Tau2"][3]));
	}
      
      if (config.find(var_prefix+"Tau3")!=config.end())
	{
	  signal_Tau3[jet_idx] = new TH1F(std::string(config[var_prefix+"Tau3"][0]+"_WP").c_str(), std::string(config[var_prefix+"Tau3"][0]+"_WP").c_str(), std::stoi(config[var_prefix+"Tau3"][1]), std::stod(config[var_prefix+"Tau3"][2]) ,std::stod(config[var_prefix+"Tau3"][3]));
	  bkg_Tau3[jet_idx] = new TH1F(std::string(config[var_prefix+"Tau3"][0]+"_QCD").c_str(), std::string(config[var_prefix+"Tau3"][0]+"_QCD").c_str(), std::stoi(config[var_prefix+"Tau3"][1]), std::stod(config[var_prefix+"Tau3"][2]) ,std::stod(config[var_prefix+"Tau3"][3]));
	}
      
      if (config.find(var_prefix+"WIDTH")!=config.end())
	{
	  signal_WIDTH[jet_idx] = new TH1F(std::string(config[var_prefix+"WIDTH"][0]+"_WP").c_str(), std::string(config[var_prefix+"WIDTH"][0]+"_WP").c_str(), std::stoi(config[var_prefix+"WIDTH"][1]), std::stod(config[var_prefix+"WIDTH"][2]) ,std::stod(config[var_prefix+"WIDTH"][3]));
	  bkg_WIDTH[jet_idx] = new TH1F(std::string(config[var_prefix+"WIDTH"][0]+"_QCD").c_str(), std::string(config[var_prefix+"WIDTH"][0]+"_QCD").c_str(), std::stoi(config[var_prefix+"WIDTH"][1]), std::stod(config[var_prefix+"WIDTH"][2]) ,std::stod(config[var_prefix+"WIDTH"][3]));
	}
      
      if (config.find(var_prefix+"SPLIT12")!=config.end())
	{
	  signal_SPLIT12[jet_idx] = new TH1F(std::string(config[var_prefix+"SPLIT12"][0]+"_WP").c_str(), std::string(config[var_prefix+"SPLIT12"][0]+"_WP").c_str(), std::stoi(config[var_prefix+"SPLIT12"][1]), std::stod(config[var_prefix+"SPLIT12"][2]) ,std::stod(config[var_prefix+"SPLIT12"][3]));
	  bkg_SPLIT12[jet_idx] = new TH1F(std::string(config[var_prefix+"SPLIT12"][0]+"_QCD").c_str(), std::string(config[var_prefix+"SPLIT12"][0]+"_QCD").c_str(), std::stoi(config[var_prefix+"SPLIT12"][1]), std::stod(config[var_prefix+"SPLIT12"][2]) ,std::stod(config[var_prefix+"SPLIT12"][3]));
	}
      if (config.find(var_prefix+"SPLIT23")!=config.end())
	{
	  signal_SPLIT23[jet_idx] = new TH1F(std::string(config[var_prefix+"SPLIT23"][0]+"_WP").c_str(), std::string(config[var_prefix+"SPLIT23"][0]+"_WP").c_str(), std::stoi(config[var_prefix+"SPLIT23"][1]), std::stod(config[var_prefix+"SPLIT23"][2]) ,std::stod(config[var_prefix+"SPLIT23"][3]));
	  bkg_SPLIT23[jet_idx] = new TH1F(std::string(config[var_prefix+"SPLIT23"][0]+"_QCD").c_str(), std::string(config[var_prefix+"SPLIT23"][0]+"_QCD").c_str(), std::stoi(config[var_prefix+"SPLIT23"][1]), std::stod(config[var_prefix+"SPLIT23"][2]) ,std::stod(config[var_prefix+"SPLIT23"][3]));
	}	  
      if (config.find(var_prefix+"SPLIT34")!=config.end())
	{
	  signal_SPLIT34[jet_idx] = new TH1F(std::string(config[var_prefix+"SPLIT34"][0]+"_WP").c_str(), std::string(config[var_prefix+"SPLIT34"][0]+"_WP").c_str(), std::stoi(config[var_prefix+"SPLIT34"][1]), std::stod(config[var_prefix+"SPLIT34"][2]) ,std::stod(config[var_prefix+"SPLIT34"][3]));
	  bkg_SPLIT34[jet_idx] = new TH1F(std::string(config[var_prefix+"SPLIT34"][0]+"_QCD").c_str(), std::string(config[var_prefix+"SPLIT34"][0]+"_QCD").c_str(), std::stoi(config[var_prefix+"SPLIT34"][1]), std::stod(config[var_prefix+"SPLIT34"][2]) ,std::stod(config[var_prefix+"SPLIT34"][3]));
	}
      
      if (config.find(var_prefix+"Dip12")!=config.end())
	{
	  signal_Dip12[jet_idx] = new TH1F(std::string(config[var_prefix+"Dip12"][0]+"_WP").c_str(), std::string(config[var_prefix+"Dip12"][0]+"_WP").c_str(), std::stoi(config[var_prefix+"Dip12"][1]), std::stod(config[var_prefix+"Dip12"][2]) ,std::stod(config[var_prefix+"Dip12"][3]));
	  bkg_Dip12[jet_idx] = new TH1F(std::string(config[var_prefix+"Dip12"][0]+"_QCD").c_str(), std::string(config[var_prefix+"Dip12"][0]+"_QCD").c_str(), std::stoi(config[var_prefix+"Dip12"][1]), std::stod(config[var_prefix+"Dip12"][2]) ,std::stod(config[var_prefix+"Dip12"][3]));
	}
      
      if (config.find(var_prefix+"Dip13")!=config.end())
	{
	  signal_Dip13[jet_idx] = new TH1F(std::string(config[var_prefix+"Dip13"][0]+"_WP").c_str(), std::string(config[var_prefix+"Dip13"][0]+"_WP").c_str(), std::stoi(config[var_prefix+"Dip13"][1]), std::stod(config[var_prefix+"Dip13"][2]) ,std::stod(config[var_prefix+"Dip13"][3]));
	  bkg_Dip13[jet_idx] = new TH1F(std::string(config[var_prefix+"Dip13"][0]+"_QCD").c_str(), std::string(config[var_prefix+"Dip13"][0]+"_QCD").c_str(), std::stoi(config[var_prefix+"Dip13"][1]), std::stod(config[var_prefix+"Dip13"][2]) ,std::stod(config[var_prefix+"Dip13"][3]));
	}
      
      if (config.find(var_prefix+"Dip23")!=config.end())
	{
	  signal_Dip23[jet_idx] = new TH1F(std::string(config[var_prefix+"Dip23"][0]+"_WP").c_str(), std::string(config[var_prefix+"Dip23"][0]+"_WP").c_str(), std::stoi(config[var_prefix+"Dip23"][1]), std::stod(config[var_prefix+"Dip23"][2]) ,std::stod(config[var_prefix+"Dip23"][3]));
	  bkg_Dip23[jet_idx] = new TH1F(std::string(config[var_prefix+"Dip23"][0]+"_QCD").c_str(), std::string(config[var_prefix+"Dip23"][0]+"_QCD").c_str(), std::stoi(config[var_prefix+"Dip23"][1]), std::stod(config[var_prefix+"Dip23"][2]) ,std::stod(config[var_prefix+"Dip23"][3]));
	}
      
      if (config.find(var_prefix+"DipExcl12")!=config.end())
	{
	  signal_DipExcl12[jet_idx] = new TH1F(std::string(config[var_prefix+"DipExcl12"][0]+"_WP").c_str(), std::string(config[var_prefix+"DipExcl12"][0]+"_WP").c_str(), std::stoi(config[var_prefix+"DipExcl12"][1]), std::stod(config[var_prefix+"DipExcl12"][2]) ,std::stod(config[var_prefix+"DipExcl12"][3]));
	  bkg_DipExcl12[jet_idx] = new TH1F(std::string(config[var_prefix+"DipExcl12"][0]+"_QCD").c_str(), std::string(config[var_prefix+"DipExcl12"][0]+"_QCD").c_str(), std::stoi(config[var_prefix+"DipExcl12"][1]), std::stod(config[var_prefix+"DipExcl12"][2]) ,std::stod(config[var_prefix+"DipExcl12"][3]));
	}
      
      if (config.find(var_prefix+"PlanarFlow")!=config.end())
	{
	  signal_PlanarFlow[jet_idx] = new TH1F(std::string(config[var_prefix+"PlanarFlow"][0]+"_WP").c_str(), std::string(config[var_prefix+"PlanarFlow"][0]+"_WP").c_str(), std::stoi(config[var_prefix+"PlanarFlow"][1]), std::stod(config[var_prefix+"PlanarFlow"][2]) ,std::stod(config[var_prefix+"PlanarFlow"][3]));
	  bkg_PlanarFlow[jet_idx] = new TH1F(std::string(config[var_prefix+"PlanarFlow"][0]+"_QCD").c_str(), std::string(config[var_prefix+"PlanarFlow"][0]+"_QCD").c_str(), std::stoi(config[var_prefix+"PlanarFlow"][1]), std::stod(config[var_prefix+"PlanarFlow"][2]) ,std::stod(config[var_prefix+"PlanarFlow"][3]));
	}
      
      if (config.find(var_prefix+"Angularity")!=config.end())
	{
	  signal_Angularity[jet_idx] = new TH1F(std::string(config[var_prefix+"Angularity"][0]+"_WP").c_str(), std::string(config[var_prefix+"Angularity"][0]+"_WP").c_str(), std::stoi(config[var_prefix+"Angularity"][1]), std::stod(config[var_prefix+"Angularity"][2]) ,std::stod(config[var_prefix+"Angularity"][3]));
	  bkg_Angularity[jet_idx] = new TH1F(std::string(config[var_prefix+"Angularity"][0]+"_QCD").c_str(), std::string(config[var_prefix+"Angularity"][0]+"_QCD").c_str(), std::stoi(config[var_prefix+"Angularity"][1]), std::stod(config[var_prefix+"Angularity"][2]) ,std::stod(config[var_prefix+"Angularity"][3]));
	}
      
      if (config.find(var_prefix+"QW")!=config.end())
	{
	  signal_QW[jet_idx] = new TH1F(std::string(config[var_prefix+"QW"][0]+"_WP").c_str(), std::string(config[var_prefix+"QW"][0]+"_WP").c_str(), std::stoi(config[var_prefix+"QW"][1]), std::stod(config[var_prefix+"QW"][2]) ,std::stod(config[var_prefix+"QW"][3]));
	  bkg_QW[jet_idx] = new TH1F(std::string(config[var_prefix+"QW"][0]+"_QCD").c_str(), std::string(config[var_prefix+"QW"][0]+"_QCD").c_str(), std::stoi(config[var_prefix+"QW"][1]), std::stod(config[var_prefix+"QW"][2]) ,std::stod(config[var_prefix+"QW"][3]));
	}
      
      if (config.find(var_prefix+"PullMag")!=config.end())
	{
	  signal_PullMag[jet_idx] = new TH1F(std::string(config[var_prefix+"PullMag"][0]+"_WP").c_str(), std::string(config[var_prefix+"PullMag"][0]+"_WP").c_str(), std::stoi(config[var_prefix+"PullMag"][1]), std::stod(config[var_prefix+"PullMag"][2]) ,std::stod(config[var_prefix+"PullMag"][3]));
	  bkg_PullMag[jet_idx] = new TH1F(std::string(config[var_prefix+"PullMag"][0]+"_QCD").c_str(), std::string(config[var_prefix+"PullMag"][0]+"_QCD").c_str(), std::stoi(config[var_prefix+"PullMag"][1]), std::stod(config[var_prefix+"PullMag"][2]) ,std::stod(config[var_prefix+"PullMag"][3]));
	}
      
      if (config.find(var_prefix+"PullPhi")!=config.end())
	{
	  signal_PullPhi[jet_idx] = new TH1F(std::string(config[var_prefix+"PullPhi"][0]+"_WP").c_str(), std::string(config[var_prefix+"PullPhi"][0]+"_WP").c_str(), std::stoi(config[var_prefix+"PullPhi"][1]), std::stod(config[var_prefix+"PullPhi"][2]) ,std::stod(config[var_prefix+"PullPhi"][3]));
	  bkg_PullPhi[jet_idx] = new TH1F(std::string(config[var_prefix+"PullPhi"][0]+"_QCD").c_str(), std::string(config[var_prefix+"PullPhi"][0]+"_QCD").c_str(), std::stoi(config[var_prefix+"PullPhi"][1]), std::stod(config[var_prefix+"PullPhi"][2]) ,std::stod(config[var_prefix+"PullPhi"][3]));
	}
	  
      if (config.find(var_prefix+"Pull_C00")!=config.end())
	{
	  signal_Pull_C00[jet_idx] = new TH1F(std::string(config[var_prefix+"Pull_C00"][0]+"_WP").c_str(), std::string(config[var_prefix+"Pull_C00"][0]+"_WP").c_str(), std::stoi(config[var_prefix+"Pull_C00"][1]), std::stod(config[var_prefix+"Pull_C00"][2]) ,std::stod(config[var_prefix+"Pull_C00"][3]));
	  bkg_Pull_C00[jet_idx] = new TH1F(std::string(config[var_prefix+"Pull_C00"][0]+"_QCD").c_str(), std::string(config[var_prefix+"Pull_C00"][0]+"_QCD").c_str(), std::stoi(config[var_prefix+"Pull_C00"][1]), std::stod(config[var_prefix+"Pull_C00"][2]) ,std::stod(config[var_prefix+"Pull_C00"][3]));
	}
      if (config.find(var_prefix+"Pull_C01")!=config.end())
	{
	  signal_Pull_C01[jet_idx] = new TH1F(std::string(config[var_prefix+"Pull_C01"][0]+"_WP").c_str(), std::string(config[var_prefix+"Pull_C01"][0]+"_WP").c_str(), std::stoi(config[var_prefix+"Pull_C01"][1]), std::stod(config[var_prefix+"Pull_C01"][2]) ,std::stod(config[var_prefix+"Pull_C01"][3]));
	  bkg_Pull_C01[jet_idx] = new TH1F(std::string(config[var_prefix+"Pull_C01"][0]+"_QCD").c_str(), std::string(config[var_prefix+"Pull_C01"][0]+"_QCD").c_str(), std::stoi(config[var_prefix+"Pull_C01"][1]), std::stod(config[var_prefix+"Pull_C01"][2]) ,std::stod(config[var_prefix+"Pull_C01"][3]));
	}
      if (config.find(var_prefix+"Pull_C10")!=config.end())
	{
	  signal_Pull_C10[jet_idx] = new TH1F(std::string(config[var_prefix+"Pull_C10"][0]+"_WP").c_str(), std::string(config[var_prefix+"Pull_C10"][0]+"_WP").c_str(), std::stoi(config[var_prefix+"Pull_C10"][1]), std::stod(config[var_prefix+"Pull_C10"][2]) ,std::stod(config[var_prefix+"Pull_C10"][3]));
	  bkg_Pull_C10[jet_idx] = new TH1F(std::string(config[var_prefix+"Pull_C10"][0]+"_QCD").c_str(), std::string(config[var_prefix+"Pull_C10"][0]+"_QCD").c_str(), std::stoi(config[var_prefix+"Pull_C10"][1]), std::stod(config[var_prefix+"Pull_C10"][2]) ,std::stod(config[var_prefix+"Pull_C10"][3]));
	}
      if (config.find(var_prefix+"Pull_C11")!=config.end())
	{
	  signal_Pull_C11[jet_idx] = new TH1F(std::string(config[var_prefix+"Pull_C11"][0]+"_WP").c_str(), std::string(config[var_prefix+"Pull_C11"][0]+"_WP").c_str(), std::stoi(config[var_prefix+"Pull_C11"][1]), std::stod(config[var_prefix+"Pull_C11"][2]) ,std::stod(config[var_prefix+"Pull_C11"][3]));
	  bkg_Pull_C11[jet_idx] = new TH1F(std::string(config[var_prefix+"Pull_C11"][0]+"_QCD").c_str(), std::string(config[var_prefix+"Pull_C11"][0]+"_QCD").c_str(), std::stoi(config[var_prefix+"Pull_C11"][1]), std::stod(config[var_prefix+"Pull_C11"][2]) ,std::stod(config[var_prefix+"Pull_C11"][3]));
	}
    }// end for loop
} // createHistograms()



void plotter::readAlgorithmConfig(std::string &algorithm)
{
    // Should read in from xml but no time, so txt file it is....
    /*
      Reads in the config file for an algorithm (csv file named ALGORITHM.config).  Gets the name of the variable, so E, m, pt, and then the name of that variable in the ntuple, and the histogram parameters.

      Keyword args:
      algorithm --- name of the algorithm we are running on

     */

    ifstream in(std::string(algorithm+".config"));
    std::string line;
    std::string delimeter = ",";
    while (getline(in,line))
      {
	// line will contain the variable type (E, m, etc), varname, the number of bins, the min and max - all comma separated
	line = rtrim(line); // remove end of line whitespace
	size_t pos = 0;
	std::string token;
	std::string curr_token;
	bool first = true; // initialise map entry correctly
	while ((pos = line.find(delimeter)) != std::string::npos) {
	  token = line.substr(0, pos);
	  if (!first) // need to init the map for token.  However, we could just do the first thing of the while loop to do this... sigh
	    {
	      config[curr_token].push_back(token);  // map of vectors of strings
	    }
	  else
	    {
	      config[token];
	      curr_token = token;
	      first = false;
	    }
	  
	  line.erase(0, pos + delimeter.length());
	}
      }
    in.close();
    getVariables(); // variables we are running over - tau1, tau2, mass, etc.  These come from data read in by the readAlgorithmConfig() method, so it _must_ be run after that method.  This should be a protected method....
  } // readAlgorithmConfig()


void plotter::getWeights(std::string & algorithm, std::string & params)
{
  /*
    Read in the weights for all of the samples.  Weights include cross section, k-factor and filter efficiency weights.  The are stored according to Run Number.
    All read in from weightings.csv.

    It creates and fills the weighting map, which is a map of Int_t (the RunNumber) containing a vector of weights, which are stored as strings here, later converted to float using std::stod().

    The number of events, weighted by mc_event_weight is read into events[] maphere.

    Keyword args:
    algorithm --- algorithm name
    params --- algorithm parameters, like pt bin name
   */


  ifstream f("weightings.csv");
  std::string line;
  std::string delimeter = ",";
  while (getline(f,line))
    {
      line = rtrim(line);
      size_t pos = 0;
      std::string token;
      Int_t curr_token;
      bool first = true;
      while ((pos = line.find(delimeter)) != std::string::npos) {
	token = line.substr(0, pos);
	if (!first) // need to init the map for token.  However, we could just do the first thing of the while loop to do this... sigh
	  {
	    // assume lumi = 20100 (from Sam's code)
	    weighting[curr_token].push_back(token); // contains vector of name,xs,k,eff
	  }
	else
	  {
	    Int_t idx = std::stoi(token);
	    weighting[idx]; // use the RunNumber as the map index
	    curr_token = idx;
	    first = false;
	  }
	
	line.erase(0, pos + delimeter.length());
      }
    }
  f.close();


  std::cout << algorithm << std::endl;
  std::vector<std::string> files = {algorithm+"_"+params+"_sig",algorithm+"_"+params+"_bkg"};
  for (std::vector<std::string>::iterator it = files.begin(); it!=files.end(); it++)
    {
      std::cout << (*it) << std::endl;
      ifstream evts((*it)+".nevents");
  
      while (getline(evts,line))
	{
	  line = rtrim(line);
	  size_t pos = 0;
	  pos = line.find(delimeter);
	  //std::cout << line.substr(pos, line.length()) << std::endl;
	  events[std::stoi(line.substr(0,pos))] = std::stod(line.substr(pos+1, line.length()));
	  //std::cout << line << std::endl;
	}
      evts.close();
    }
  
} // getWeights()





// need branches and th1f stuff still

void plotter::setBranches()
  {
    // arg, need to do this differently because we have vectors of TH1F * ......
    // in fact, we have to do this for truth/topo/groomed..........
  
    sig->SetBranchAddress("normalisation",&sig_normalisation);
    sig->SetBranchAddress("leadTruthIndex",&signal_truth_index);
    sig->SetBranchAddress("leadTopoIndex",&signal_topo_index);
    sig->SetBranchAddress("leadGroomedIndex",&signal_groomed_index);
    sig->SetBranchAddress("mc_event_weight",&signal_mc_event_weight);
    sig->SetBranchAddress("RunNumber",&signal_RunNumber);
    sig->SetBranchAddress("mc_channel_number",&signal_RunNumber);
    //sig->SetBranchAddress("NEvents_weighted",&signal_NEvents_weighted);
    sig->SetBranchAddress("NEvents",&signal_NEvents);

    bkg->SetBranchAddress("normalisation",&bkg_normalisation);
    bkg->SetBranchAddress("leadTruthIndex",&bkg_truth_index);
    bkg->SetBranchAddress("leadTopoIndex",&bkg_topo_index);
    bkg->SetBranchAddress("leadGroomedIndex",&bkg_groomed_index);
    bkg->SetBranchAddress("mc_event_weight",&bkg_mc_event_weight);
    bkg->SetBranchAddress("RunNumber",&bkg_RunNumber);
    bkg->SetBranchAddress("mc_channel_number",&bkg_RunNumber);
    //bkg->SetBranchAddress("NEvents_weighted",&bkg_NEvents_weighted);
    bkg->SetBranchAddress("NEvents",&bkg_NEvents);

    const char * str[] = {"truth","topo", "groomed"};
    // create a vector containing the jetTypes so we can make an iterator...
    std::vector<std::string> jetTypes(str, std::end(str));
    // and a map for easy referencing
    std::map<std::string, int> jetMap = {{"truth",jet_type::TRUTH},{"topo",jet_type::TOPO},{"groomed",jet_type::GROOMED}};

    for (std::vector<std::string>::iterator it = jetTypes.begin(); it != jetTypes.end(); it++)
      {
	sig->SetBranchAddress(vars[(*it)+"_E"].c_str(), &signal_E_vec[jetMap.at((*it))]);
	sig->SetBranchAddress(vars[(*it)+"_pt"].c_str(), &signal_pt_vec[jetMap.at((*it))]);
	sig->SetBranchAddress(vars[(*it)+"_m"].c_str(), &signal_m_vec[jetMap.at((*it))]);
	sig->SetBranchAddress(vars[(*it)+"_eta"].c_str(), &signal_eta_vec[jetMap.at((*it))]);
	sig->SetBranchAddress(vars[(*it)+"_phi"].c_str(), &signal_phi_vec[jetMap.at((*it))]);
	if (jetMap.at((*it)) != 0) // there is no truth emfrac...
	  {
	    sig->SetBranchAddress(vars[(*it)+"_emfrac"].c_str(), &signal_emfrac_vec[jetMap.at((*it))]);
	  }
	sig->SetBranchAddress(vars[(*it)+"_Tau1"].c_str(), &signal_Tau1_vec[jetMap.at((*it))]);
	sig->SetBranchAddress(vars[(*it)+"_Tau2"].c_str(), &signal_Tau2_vec[jetMap.at((*it))]);
	sig->SetBranchAddress(vars[(*it)+"_Tau3"].c_str(), &signal_Tau3_vec[jetMap.at((*it))]);
	sig->SetBranchAddress(vars[(*it)+"_WIDTH"].c_str(), &signal_WIDTH_vec[jetMap.at((*it))]);
	sig->SetBranchAddress(vars[(*it)+"_SPLIT12"].c_str(), &signal_SPLIT12_vec[jetMap.at((*it))]);
	sig->SetBranchAddress(vars[(*it)+"_SPLIT23"].c_str(), &signal_SPLIT23_vec[jetMap.at((*it))]);
	sig->SetBranchAddress(vars[(*it)+"_SPLIT34"].c_str(), &signal_SPLIT34_vec[jetMap.at((*it))]);
	sig->SetBranchAddress(vars[(*it)+"_Dip12"].c_str(), &signal_Dip12_vec[jetMap.at((*it))]);
	sig->SetBranchAddress(vars[(*it)+"_Dip13"].c_str(), &signal_Dip13_vec[jetMap.at((*it))]);
	sig->SetBranchAddress(vars[(*it)+"_Dip23"].c_str(), &signal_Dip23_vec[jetMap.at((*it))]);
	sig->SetBranchAddress(vars[(*it)+"_DipExcl12"].c_str(), &signal_DipExcl12_vec[jetMap.at((*it))]);
	sig->SetBranchAddress(vars[(*it)+"_PlanarFlow"].c_str(), &signal_PlanarFlow_vec[jetMap.at((*it))]);
	sig->SetBranchAddress(vars[(*it)+"_Angularity"].c_str(), &signal_Angularity_vec[jetMap.at((*it))]);
	sig->SetBranchAddress(vars[(*it)+"_QW"].c_str(), &signal_QW_vec[jetMap.at((*it))]);
	sig->SetBranchAddress(vars[(*it)+"_PullMag"].c_str(), &signal_PullMag_vec[jetMap.at((*it))]);
	sig->SetBranchAddress(vars[(*it)+"_PullPhi"].c_str(), &signal_PullPhi_vec[jetMap.at((*it))]);
	sig->SetBranchAddress(vars[(*it)+"_Pull_C00"].c_str(), &signal_Pull_C00_vec[jetMap.at((*it))]);
	sig->SetBranchAddress(vars[(*it)+"_Pull_C01"].c_str(), &signal_Pull_C01_vec[jetMap.at((*it))]);
	sig->SetBranchAddress(vars[(*it)+"_Pull_C10"].c_str(), &signal_Pull_C10_vec[jetMap.at((*it))]);
	sig->SetBranchAddress(vars[(*it)+"_Pull_C11"].c_str(), &signal_Pull_C11_vec[jetMap.at((*it))]);

	bkg->SetBranchAddress(vars[(*it)+"_E"].c_str(), &bkg_E_vec[jetMap.at((*it))]);
	bkg->SetBranchAddress(vars[(*it)+"_pt"].c_str(), &bkg_pt_vec[jetMap.at((*it))]);
	bkg->SetBranchAddress(vars[(*it)+"_m"].c_str(), &bkg_m_vec[jetMap.at((*it))]);
	bkg->SetBranchAddress(vars[(*it)+"_eta"].c_str(), &bkg_eta_vec[jetMap.at((*it))]);
	bkg->SetBranchAddress(vars[(*it)+"_phi"].c_str(), &bkg_phi_vec[jetMap.at((*it))]);
	if (jetMap.at((*it))!=0) // ther eis no truth emfrac
	  bkg->SetBranchAddress(vars[(*it)+"_emfrac"].c_str(), &bkg_emfrac_vec[jetMap.at((*it))]);
	bkg->SetBranchAddress(vars[(*it)+"_Tau1"].c_str(), &bkg_Tau1_vec[jetMap.at((*it))]);
	bkg->SetBranchAddress(vars[(*it)+"_Tau2"].c_str(), &bkg_Tau2_vec[jetMap.at((*it))]);
	bkg->SetBranchAddress(vars[(*it)+"_Tau3"].c_str(), &bkg_Tau3_vec[jetMap.at((*it))]);
	bkg->SetBranchAddress(vars[(*it)+"_WIDTH"].c_str(), &bkg_WIDTH_vec[jetMap.at((*it))]);
	bkg->SetBranchAddress(vars[(*it)+"_SPLIT12"].c_str(), &bkg_SPLIT12_vec[jetMap.at((*it))]);
	bkg->SetBranchAddress(vars[(*it)+"_SPLIT23"].c_str(), &bkg_SPLIT23_vec[jetMap.at((*it))]);
	bkg->SetBranchAddress(vars[(*it)+"_SPLIT34"].c_str(), &bkg_SPLIT34_vec[jetMap.at((*it))]);
	bkg->SetBranchAddress(vars[(*it)+"_Dip12"].c_str(), &bkg_Dip12_vec[jetMap.at((*it))]);
	bkg->SetBranchAddress(vars[(*it)+"_Dip13"].c_str(), &bkg_Dip13_vec[jetMap.at((*it))]);
	bkg->SetBranchAddress(vars[(*it)+"_Dip23"].c_str(), &bkg_Dip23_vec[jetMap.at((*it))]);
	bkg->SetBranchAddress(vars[(*it)+"_DipExcl12"].c_str(), &bkg_DipExcl12_vec[jetMap.at((*it))]);
	bkg->SetBranchAddress(vars[(*it)+"_PlanarFlow"].c_str(), &bkg_PlanarFlow_vec[jetMap.at((*it))]);
	bkg->SetBranchAddress(vars[(*it)+"_Angularity"].c_str(), &bkg_Angularity_vec[jetMap.at((*it))]);
	bkg->SetBranchAddress(vars[(*it)+"_QW"].c_str(), &bkg_QW_vec[jetMap.at((*it))]);
	bkg->SetBranchAddress(vars[(*it)+"_PullMag"].c_str(), &bkg_PullMag_vec[jetMap.at((*it))]);
	bkg->SetBranchAddress(vars[(*it)+"_PullPhi"].c_str(), &bkg_PullPhi_vec[jetMap.at((*it))]);
	bkg->SetBranchAddress(vars[(*it)+"_Pull_C00"].c_str(), &bkg_Pull_C00_vec[jetMap.at((*it))]);
	bkg->SetBranchAddress(vars[(*it)+"_Pull_C01"].c_str(), &bkg_Pull_C01_vec[jetMap.at((*it))]);
	bkg->SetBranchAddress(vars[(*it)+"_Pull_C10"].c_str(), &bkg_Pull_C10_vec[jetMap.at((*it))]);
	bkg->SetBranchAddress(vars[(*it)+"_Pull_C11"].c_str(), &bkg_Pull_C11_vec[jetMap.at((*it))]);
      } // end of iterator over varname  
  } //setBranches()


void plotter::compareHistograms(TH1F * hist1, TH1F * hist2, std::string file_name)
{
  /*
    Method to compare two histograms and draw them on the same canvas.  This will save the canvas to a .png.  The name of the
    file is given by the argument file_name.

    Keyword args:
    TH1F * hist1 --- first histogram to be plotted
    TH1F * hist2 --- second histogram to be plotted
    string file_name --- file name of the output file
   */

  hist1->SetLineColor(kBlue);
  hist1->Sumw2();
  //hist1->Scale(1./hist1->Integral());
  //hist1->GetXaxis()->SetLabel();
  hist1->GetYaxis()->SetTitle("# Entries");
  hist1->SetMarkerSize(3);
  hist1->SetLineWidth(2);
  hist2->SetLineColor(kRed);
  hist2->Sumw2();
  //hist2->Scale(1./hist2->Integral());
  //hist2->GetXaxis()->SetLabel();
  hist2->GetYaxis()->SetTitle("# Entries");
  hist2->SetMarkerSize(3);
  hist2->SetLineWidth(2);


  int draw_first = 1;
  if (hist2->GetMaximum() > hist1->GetMaximum())
    draw_first = 2;

  TLegend * leg = new TLegend();
  leg->AddEntry(hist1);
  leg->AddEntry(hist2);

  TCanvas * c = new TCanvas();

  if (draw_first == 1)
    {
      hist1->Draw("E");
      hist2->Draw("Esame");
    }
  else
    {
      hist2->Draw("E");
      hist1->Draw("Esame");
    }
  leg->Draw("same");
  c->SaveAs(std::string(working_dir+file_name+".png").c_str());
  delete c;
  delete leg;

}// compareHistograms


void plotter::fillHistograms(TTree * tree, bool signal)
  {
    /*
      This runs over a tree and gets all of the entries and fills the histograms.  The weights are applied at this stage.

      All TH1Fs for the variables are filled here.

      Keyword args:
      tree --- TTree we're running over
      signal --- whether this is a signal or bkg sample
     */

    int nentries = (int)tree->GetEntries();
    double weight = 1;

    std::cout <<"nentries: " << nentries << std::endl;

    for (int i = 0; i < nentries; i++)
      {
	tree->GetEntry(i);

	weight = signal ? signal_mc_event_weight : bkg_mc_event_weight;
	Int_t runNum = signal ? signal_RunNumber : bkg_RunNumber;

	// the number of events weighted by the mc_evt for the FULL sample, before mass window cuts
	//double nevents_weighted = signal ? signal_NEvents_weighted : bkg_NEvents_weighted;
	double nevents_original = signal ? signal_NEvents : bkg_NEvents;
	//std::cout << "NEVENTS: " << nevents_original << std::endl;
	//std::cout << "NEVENTS_W: " << nevents_weighted << std::endl;
	//if (runNum == 195847) // this is a special W+jets run of a few hundred events.  There is no weighting info on the twiki and I'm not sure how to handle it
	    //continue;
	
	// 20100 is the lumi, weighting[runNum][1] is the cross section, weighting[runNum][2] is the k-factor, weighting[runNum][3] is the filter efficiency
	if (weighting.find(runNum)!=weighting.end())
	  {
	    weight*=(20100*std::stod(weighting[runNum][3])*std::stod(weighting[runNum][2]));
	    weight*= std::stod(weighting[runNum][1]);
	    // normalise the weighting
	    //std::cout << events[runNum] << std::endl;
	    weight/=events[runNum];
	    
	  }

	

	double pt = (*signal_pt_vec[jet_type::TRUTH])[signal_truth_index]/1000.0;
	
	if (applyPtCut &&  (pt < 200 || pt > 350))
	  {
	    continue;
	  }	

	/*if (applyEtaCut && fabs((*signal_eta_vec[jet_type::TRUTH])[sig_truth_idx]) < 1.2)
	  {
	    continue;
	    }*/

	// Apply pt_reweighting....
	if (signal)
	  {
	    //double pt = signal ? (*signal_pt_vec[0])[signal_truth_index]:(*bkg_pt_vec[0])[bkg_truth_index];
	    weight*=getPtWeight(pt);
	    // does this need to be unnormalised somehow? :/ *rewight->Integral()?
	    //Int_t pt_idx = pt_reweight_hist->GetXaxis()->FindBin(pt);
	    //std::cout << "BEFORE " << weight << std::endl;
	    //weight *= (pt_reweight_hist->GetBinContent(pt_idx));///nevents_original);
	    //std::cout << "AFTER " << weight << std::endl;
	  }
	// fill all the histograms for truth,topo.groomed
	for (int j = 0; j < 3; j++) // for truth/topo/groomed
	  {
	    // The lead index is already calculated and stored.
	    int sig_lead_idx = 0;
	    int bkg_lead_idx = 0;
	    std::string var_prefix = "";
	    switch (j)
	      {
	      case jet_type::TRUTH:
		sig_lead_idx = signal_truth_index;
		bkg_lead_idx = bkg_truth_index;
		var_prefix = "truth_";
		break;
	      case jet_type::TOPO:
		sig_lead_idx = signal_topo_index;
		bkg_lead_idx = bkg_topo_index;
		var_prefix = "topo_";
		break;
	      case jet_type::GROOMED:
		sig_lead_idx = signal_groomed_index;
		bkg_lead_idx = bkg_groomed_index;
		var_prefix = "groomed_";
		break;
	      }
	    //TODO arg, apparently having no topo jet index when not matching truth to reco is the only time it means anything, the rest of the time we don't care!!! So we don't even really need the topo jets... FML
	    sig_lead_idx = sig_lead_idx < 0 ? 0 : sig_lead_idx;
	    bkg_lead_idx = bkg_lead_idx < 0 ? 0 : bkg_lead_idx;




	    if (signal)
	    {
	      	if (config.find(var_prefix+"E")!=config.end())
		  signal_E[j]->Fill((*signal_E_vec[j])[sig_lead_idx]/1000.0, weight);
	      	if (config.find(var_prefix+"pt")!=config.end())
		  signal_pt[j]->Fill((*signal_pt_vec[j])[sig_lead_idx]/1000.0, weight);
	      	if (config.find(var_prefix+"m")!=config.end())
		  signal_m[j]->Fill((*signal_m_vec[j])[sig_lead_idx]/1000.0, weight);
	      	if (config.find(var_prefix+"eta")!=config.end())
		  signal_eta[j]->Fill((*signal_eta_vec[j])[sig_lead_idx], weight);
	      	if (config.find(var_prefix+"phi")!=config.end())
		  signal_phi[j]->Fill((*signal_phi_vec[j])[sig_lead_idx], weight);
	      	if (config.find(var_prefix+"E")!=config.end() && j!= jet_type::TRUTH)
		    signal_emfrac[j]->Fill((*signal_emfrac_vec[j])[sig_lead_idx]);
	      	if (config.find(var_prefix+"Tau1")!=config.end())
		  signal_Tau1[j]->Fill((*signal_Tau1_vec[j])[sig_lead_idx], weight);
		if (config.find(var_prefix+"Tau2")!=config.end())
		  signal_Tau2[j]->Fill((*signal_Tau2_vec[j])[sig_lead_idx], weight);
		if (config.find(var_prefix+"Tau3")!=config.end())
		  signal_Tau3[j]->Fill((*signal_Tau3_vec[j])[sig_lead_idx], weight);
		if (config.find(var_prefix+"WIDTH")!=config.end())
		  signal_WIDTH[j]->Fill((*signal_WIDTH_vec[j])[sig_lead_idx], weight);
		if (config.find(var_prefix+"SPLIT12")!=config.end())
		  signal_SPLIT12[j]->Fill((*signal_SPLIT12_vec[j])[sig_lead_idx], weight);
		if (config.find(var_prefix+"SPLIT23")!=config.end())
		  signal_SPLIT23[j]->Fill((*signal_SPLIT23_vec[j])[sig_lead_idx], weight);
		if (config.find(var_prefix+"SPLIT34")!=config.end())
		  signal_SPLIT34[j]->Fill((*signal_SPLIT34_vec[j])[sig_lead_idx], weight);
		if (config.find(var_prefix+"Dip12")!=config.end())
		  signal_Dip12[j]->Fill((*signal_Dip12_vec[j])[sig_lead_idx], weight);
		if (config.find(var_prefix+"Dip13")!=config.end())
		  signal_Dip13[j]->Fill((*signal_Dip13_vec[j])[sig_lead_idx], weight);
		if (config.find(var_prefix+"Dip23")!=config.end())
		  signal_Dip23[j]->Fill((*signal_Dip23_vec[j])[sig_lead_idx], weight);
		if (config.find(var_prefix+"DipExcl12")!=config.end())
		  signal_DipExcl12[j]->Fill((*signal_DipExcl12_vec[j])[sig_lead_idx], weight);
		if (config.find(var_prefix+"PlanarFlow")!=config.end())
		  signal_PlanarFlow[j]->Fill((*signal_PlanarFlow_vec[j])[sig_lead_idx], weight);
		if (config.find(var_prefix+"Angularity")!=config.end())
		  signal_Angularity[j]->Fill((*signal_Angularity_vec[j])[sig_lead_idx], weight);
		if (config.find(var_prefix+"QW")!=config.end())
		  signal_QW[j]->Fill((*signal_QW_vec[j])[sig_lead_idx], weight);
		if (config.find(var_prefix+"PullMag")!=config.end())
		  signal_PullMag[j]->Fill((*signal_PullMag_vec[j])[sig_lead_idx], weight);
		if (config.find(var_prefix+"PullPhi")!=config.end())
		  signal_PullPhi[j]->Fill((*signal_PullPhi_vec[j])[sig_lead_idx], weight);
		if (config.find(var_prefix+"Pull_C00")!=config.end())
		  signal_Pull_C00[j]->Fill((*signal_Pull_C00_vec[j])[sig_lead_idx], weight);
		if (config.find(var_prefix+"Pull_C01")!=config.end())
		  signal_Pull_C01[j]->Fill((*signal_Pull_C01_vec[j])[sig_lead_idx], weight);
		if (config.find(var_prefix+"Pull_C10")!=config.end())
		  signal_Pull_C10[j]->Fill((*signal_Pull_C10_vec[j])[sig_lead_idx], weight);
		if (config.find(var_prefix+"Pull_C11")!=config.end())
		  signal_Pull_C11[j]->Fill((*signal_Pull_C11_vec[j])[sig_lead_idx], weight);
		}
	    else
	    {
	      if (config.find(var_prefix+"E")!=config.end())
		bkg_E[j]->Fill((*bkg_E_vec[j])[bkg_lead_idx]/1000.0, weight);
	      if (config.find(var_prefix+"pt")!=config.end())
		bkg_pt[j]->Fill((*bkg_pt_vec[j])[bkg_lead_idx]/1000.0, weight);
	      if (config.find(var_prefix+"m")!=config.end())
		bkg_m[j]->Fill((*bkg_m_vec[j])[bkg_lead_idx]/1000.0, weight);
	      if (config.find(var_prefix+"eta")!=config.end())
		bkg_eta[j]->Fill((*bkg_eta_vec[j])[bkg_lead_idx], weight);
	      if (config.find(var_prefix+"phi")!=config.end())
		bkg_phi[j]->Fill((*bkg_phi_vec[j])[bkg_lead_idx], weight);
	      if (config.find(var_prefix+"emfrac")!=config.end() && j!=jet_type::TRUTH)
		bkg_emfrac[j]->Fill((*bkg_emfrac_vec[j])[bkg_lead_idx]);
	      if (config.find(var_prefix+"Tau1")!=config.end())
		bkg_Tau1[j]->Fill((*bkg_Tau1_vec[j])[bkg_lead_idx], weight);
	      if (config.find(var_prefix+"Tau2")!=config.end())
		bkg_Tau2[j]->Fill((*bkg_Tau2_vec[j])[bkg_lead_idx], weight);
	      if (config.find(var_prefix+"Tau3")!=config.end())
		bkg_Tau3[j]->Fill((*bkg_Tau3_vec[j])[bkg_lead_idx], weight);
	      if (config.find(var_prefix+"WIDTH")!=config.end())
		bkg_WIDTH[j]->Fill((*bkg_WIDTH_vec[j])[bkg_lead_idx], weight);
	      if (config.find(var_prefix+"SPLIT12")!=config.end())
		bkg_SPLIT12[j]->Fill((*bkg_SPLIT12_vec[j])[bkg_lead_idx], weight);
	      if (config.find(var_prefix+"SPLIT23")!=config.end())
		bkg_SPLIT23[j]->Fill((*bkg_SPLIT23_vec[j])[bkg_lead_idx], weight);
	      if (config.find(var_prefix+"SPLIT34")!=config.end())
		bkg_SPLIT34[j]->Fill((*bkg_SPLIT34_vec[j])[bkg_lead_idx], weight);
	      if (config.find(var_prefix+"Dip12")!=config.end())
		bkg_Dip12[j]->Fill((*bkg_Dip12_vec[j])[bkg_lead_idx], weight);
	      if (config.find(var_prefix+"Dip13")!=config.end())
		bkg_Dip13[j]->Fill((*bkg_Dip13_vec[j])[bkg_lead_idx], weight);
	      if (config.find(var_prefix+"Dip23")!=config.end())
		bkg_Dip23[j]->Fill((*bkg_Dip23_vec[j])[bkg_lead_idx], weight);
	      if (config.find(var_prefix+"DipExcl12")!=config.end())
		bkg_DipExcl12[j]->Fill((*bkg_DipExcl12_vec[j])[bkg_lead_idx], weight);
	      if (config.find(var_prefix+"PlanarFlow")!=config.end())
		bkg_PlanarFlow[j]->Fill((*bkg_PlanarFlow_vec[j])[bkg_lead_idx], weight);
	      if (config.find(var_prefix+"Angularity")!=config.end())
		bkg_Angularity[j]->Fill((*bkg_Angularity_vec[j])[bkg_lead_idx], weight);
	      if (config.find(var_prefix+"QW")!=config.end())
		bkg_QW[j]->Fill((*bkg_QW_vec[j])[bkg_lead_idx], weight);
	      if (config.find(var_prefix+"PullMag")!=config.end())
		bkg_PullMag[j]->Fill((*bkg_PullMag_vec[j])[bkg_lead_idx], weight);
	      if (config.find(var_prefix+"PullPhi")!=config.end())
		bkg_PullPhi[j]->Fill((*bkg_PullPhi_vec[j])[bkg_lead_idx], weight);
	      if (config.find(var_prefix+"Pull_C00")!=config.end())
		bkg_Pull_C00[j]->Fill((*bkg_Pull_C00_vec[j])[bkg_lead_idx], weight);
	      if (config.find(var_prefix+"Pull_C01")!=config.end())
		bkg_Pull_C01[j]->Fill((*bkg_Pull_C01_vec[j])[bkg_lead_idx], weight);
	      if (config.find(var_prefix+"Pull_C10")!=config.end())
		bkg_Pull_C10[j]->Fill((*bkg_Pull_C10_vec[j])[bkg_lead_idx], weight);
	      if (config.find(var_prefix+"Pull_C11")!=config.end())
		bkg_Pull_C11[j]->Fill((*bkg_Pull_C11_vec[j])[bkg_lead_idx], weight);
	    }
	  } // end of loop j
      }

  } // fillHistograms()



void plotter::plotHistograms()
  {
    //pretty plots will be made here :)
  

  }


void plotter::initHists()
{
  // Initialise all the histogram vectors...
  
  for (int i = 0; i < 3; i++) //truth/topo.groomed
    {
  signal_E.push_back(0);
  signal_pt.push_back(0);
  signal_m.push_back(0);
  signal_eta.push_back(0);
  signal_phi.push_back(0);
  signal_emfrac.push_back(0);
  signal_Tau1.push_back(0);
  signal_Tau2.push_back(0);
  signal_Tau3.push_back(0);
  signal_WIDTH.push_back(0);
  signal_SPLIT12.push_back(0);
  signal_SPLIT23.push_back(0);
  signal_SPLIT34.push_back(0);
  signal_Dip12.push_back(0);
  signal_Dip13.push_back(0);
  signal_Dip23.push_back(0);
  signal_DipExcl12.push_back(0);
  signal_PlanarFlow.push_back(0);
  signal_Angularity.push_back(0);
  signal_QW.push_back(0);
  signal_PullMag.push_back(0);
  signal_PullPhi.push_back(0);
  signal_Pull_C00.push_back(0);
  signal_Pull_C01.push_back(0);
  signal_Pull_C10.push_back(0);
  signal_Pull_C11.push_back(0);


  bkg_E.push_back(0);
  bkg_pt.push_back(0);
  bkg_m.push_back(0);
  bkg_eta.push_back(0);
  bkg_phi.push_back(0);
  bkg_emfrac.push_back(0);
  bkg_Tau1.push_back(0);
  bkg_Tau2.push_back(0);
  bkg_Tau3.push_back(0);
  bkg_WIDTH.push_back(0);
  bkg_SPLIT12.push_back(0);
  bkg_SPLIT23.push_back(0);
  bkg_SPLIT34.push_back(0);
  bkg_Dip12.push_back(0);
  bkg_Dip13.push_back(0);
  bkg_Dip23.push_back(0);
  bkg_DipExcl12.push_back(0);
  bkg_PlanarFlow.push_back(0);
  bkg_Angularity.push_back(0);
  bkg_QW.push_back(0);
  bkg_PullMag.push_back(0);
  bkg_PullPhi.push_back(0);
  bkg_Pull_C00.push_back(0);
  bkg_Pull_C01.push_back(0);
  bkg_Pull_C10.push_back(0);
  bkg_Pull_C11.push_back(0);
    }  
} //initHists

void plotter::clearHists()
{
  
  // clear all of the histogram vectors

  signal_E.clear();
  signal_pt.clear();
  signal_m.clear();
  signal_eta.clear();
  signal_phi.clear();
  signal_emfrac.clear();
  signal_Tau1.clear();
  signal_Tau2.clear();
  signal_Tau3.clear();
  signal_WIDTH.clear();
  signal_SPLIT12.clear();
  signal_SPLIT23.clear();
  signal_SPLIT34.clear();
  signal_Dip12.clear();
  signal_Dip13.clear();
  signal_Dip23.clear();
  signal_DipExcl12.clear();
  signal_PlanarFlow.clear();
  signal_Angularity.clear();
  signal_QW.clear();
  signal_PullMag.clear();
  signal_PullPhi.clear();
  signal_Pull_C00.clear();
  signal_Pull_C01.clear();
  signal_Pull_C10.clear();
  signal_Pull_C11.clear();


  bkg_E.clear();
  bkg_pt.clear();
  bkg_m.clear();
  bkg_eta.clear();
  bkg_phi.clear();
  bkg_emfrac.clear();
  bkg_Tau1.clear();
  bkg_Tau2.clear();
  bkg_Tau3.clear();
  bkg_WIDTH.clear();
  bkg_SPLIT12.clear();
  bkg_SPLIT23.clear();
  bkg_SPLIT34.clear();
  bkg_Dip12.clear();
  bkg_Dip13.clear();
  bkg_Dip23.clear();
  bkg_DipExcl12.clear();
  bkg_PlanarFlow.clear();
  bkg_Angularity.clear();
  bkg_QW.clear();
  bkg_PullMag.clear();
  bkg_PullPhi.clear();
  bkg_Pull_C00.clear();
  bkg_Pull_C01.clear();
  bkg_Pull_C10.clear();
  bkg_Pull_C11.clear();
  
} // clearHists()

void plotter::initVectors()
  {

    // have the vectors for the above histograms so we can do the reading in stuff from the TTree
    for (int i = 0; i < 3; i++) // for jetType::TRUTH/TOPO/GROOMED
      {
	signal_E_vec[i] = 0;
	signal_pt_vec[i] = 0;
	signal_m_vec[i] = 0;
	signal_eta_vec[i] = 0;
	signal_phi_vec[i] = 0;
	signal_emfrac_vec[i] = 0;
	signal_Tau1_vec[i] = 0;
	signal_Tau2_vec[i] = 0;
	signal_Tau3_vec[i] = 0;
	signal_WIDTH_vec[i] = 0;
	signal_SPLIT12_vec[i] = 0;
	signal_SPLIT23_vec[i] = 0;
	signal_SPLIT34_vec[i] = 0;
	signal_Dip12_vec[i] = 0;
	signal_Dip13_vec[i] = 0;
	signal_Dip23_vec[i] = 0;
	signal_DipExcl12_vec[i] = 0;
	signal_PlanarFlow_vec[i] = 0;
	signal_Angularity_vec[i] = 0;
	signal_QW_vec[i] = 0;
	signal_PullMag_vec[i] = 0;
	signal_PullPhi_vec[i] = 0;
	signal_Pull_C00_vec[i] = 0;
	signal_Pull_C01_vec[i] = 0;
	signal_Pull_C10_vec[i] = 0;
	signal_Pull_C11_vec[i] = 0;



	bkg_E_vec[i] = 0;
	bkg_pt_vec[i] = 0;
	bkg_m_vec[i] = 0;
	bkg_eta_vec[i] = 0;
	bkg_phi_vec[i] = 0;
	bkg_Tau1_vec[i] = 0;
	bkg_Tau2_vec[i] = 0;
	bkg_Tau3_vec[i] = 0;
	bkg_WIDTH_vec[i] = 0;
	bkg_SPLIT12_vec[i] = 0;
	bkg_SPLIT23_vec[i] = 0;
	bkg_SPLIT34_vec[i] = 0;
	bkg_Dip12_vec[i] = 0;
	bkg_Dip13_vec[i] = 0;
	bkg_Dip23_vec[i] = 0;
	bkg_DipExcl12_vec[i] = 0;
	bkg_PlanarFlow_vec[i] = 0;
	bkg_Angularity_vec[i] = 0;
	bkg_QW_vec[i] = 0;
	bkg_PullMag_vec[i] = 0;
	bkg_PullPhi_vec[i] = 0;
	bkg_Pull_C00_vec[i] = 0;
	bkg_Pull_C01_vec[i] = 0;
	bkg_Pull_C10_vec[i] = 0;
	bkg_Pull_C11_vec[i]= 0;//std::map<i, 0>;
  
      }  
  } // initVectors()

float plotter::getPtWeight(double & pt_for_binning)
{
  double scalept = 1.0;
  if     (pt_for_binning <= 225.0) scalept = 0.27;

  else if(pt_for_binning <= 250.0) scalept = 0.32;

  else if(pt_for_binning <= 275.0) scalept = 0.73;

  else if(pt_for_binning <= 300.0) scalept = 1.72;

  else if(pt_for_binning <= 335.0) scalept = 2.09;

  else if(pt_for_binning <= 350.0) scalept = 2.34;

  else if(pt_for_binning <= 400.0) scalept = 0.42;

  else if(pt_for_binning <= 450.0) scalept = 0.73;

  else if(pt_for_binning <= 500.0) scalept = 1.33;

  else if(pt_for_binning <= 600.0) scalept = 0.35;

  else if(pt_for_binning <= 700.0) scalept = 1.00;

  else if(pt_for_binning <= 800.0) scalept = 2.55;

  else if(pt_for_binning <= 900.0) scalept = 5.16;

  else if(pt_for_binning <= 1000.0) scalept = 15.15;

  return scalept;
} //getPtWeight
