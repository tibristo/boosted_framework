///
/// STILL TO DO: 
/// 
/// 
/// * add other jet algorithms once available
/// 

#include <TH1F.h>
#include <TH1I.h>
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TGraph.h"
#include "Rtypes.h"
#include <string>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <map>
#include <cmath>
#include "TH1F.h"
#include "TROOT.h"
#include "TInterpreter.h"
#include "TPostScript.h"
#include "TCanvas.h"
#include "TTreeIndex.h"
#include "TMath.h"
#include "AtlasStyle.C"
#include "TLegend.h"
#include "TAxis.h"
#include "TLatex.h"
#include "TLorentzVector.h"
//FASTJET INCLUDES
/*#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/PseudoJetStructureBase.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/"
#include "fastjet/"
#include "fastjet/"*/

#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequence1GhostPassiveArea.hh"
#include "fastjet/ClusterSequenceActiveArea.hh"
#include "fastjet/ClusterSequenceActiveAreaExplicitGhosts.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/ClusterSequenceAreaBase.hh"
#include "fastjet/ClusterSequencePassiveArea.hh"
#include "fastjet/ClusterSequenceStructure.hh"
#include "fastjet/ClusterSequenceVoronoiArea.hh"
#include "fastjet/CompositeJetStructure.hh"
#include "fastjet/EECambridgePlugin.hh"
#include "fastjet/Error.hh"
#include "fastjet/FunctionOfPseudoJet.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/GridJetPlugin.hh"
#include "fastjet/JadePlugin.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/LimitedWarning.hh"
#include "fastjet/NNH.hh"
#include "fastjet/NestedDefsPlugin.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/PseudoJetStructureBase.hh"
#include "fastjet/RangeDefinition.hh"
#include "fastjet/SISConeBasePlugin.hh"
#include "fastjet/SISConePlugin.hh"
#include "fastjet/SISConeSphericalPlugin.hh"
#include "fastjet/Selector.hh"
#include "fastjet/SharedPtr.hh"
#include "fastjet/WrappedStructure.hh"
#include "fastjet/version.hh"



//Declarations
//bool CheckBranchExists( string BranchName );
float DeltaR(float eta1,float phi1,float eta2,float phi2);
double mpv(TH1F* histo);
void Qw(double &minWidth, double &topEdge, TH1F* histo, double frac);

//void cleanBranches();
//void defineStrings(TString *AlgoList, TString *binLabel, TString *pTbins, TString *finePtBins);
void createHistos();
void getMassHistograms(TTree *inputTree, TTree *inputTree1, TString groomAlgo, int groomAlgoIndex);
void initializeVariables();
void getBranches(TTree *inputTree, TTree *inputTree1, TString groomAlgo, int groomAlgoIndex);
void deleteVectors();
void getNormSherpaW(TString inputTree, unsigned long & evnum, double & weight);
void makeROC(int type, TH1F *&S,TH1F *&B,TGraph &curve, TString name="", bool draw=false);
void makePlots();
void makePtPlots();
void setMassBranch(TTree * tree, std::string &algorithm);
void plotVariables(TTree * tree, vector<std:string> & branches);
std::vector<std::string> getListOfBranches(std::string &algorithm);
void make68Plots(int algidx, TTree * bkg, TTree * sig);

enum class groomAlgo{groomZero, TopoSplitFilteredMu67SmallR0YCut9, TopoSplitFilteredMu100SmallR30YCut4, TopoTrimmedPtFrac5SmallR30, TopoTrimmedPtFrac5SmallR20, TopoPrunedCaRcutFactor50Zcut10, TopoPrunedCaRcutFactor50Zcut20, AntiKt2LCTopo, AntiKt3LCTopo, AntiKt4LCTopo};
//groomAlgo options:
//TopoSplitFilteredMu67SmallR0YCut9   - 1
//TopoSplitFilteredMu100SmallR30YCut4  - 2
//TopoTrimmedPtFrac5SmallR30 - 3
// TopoTrimmedPtFrac5SmallR20 - 4
// TopoPrunedCaRcutFactor50Zcut10 - 5
// TopoPrunedCaRcutFactor50Zcut20 - 6

//AntiKt2LCTopo - 7
//AntiKt3LCTopo - 8
//AntiKt4LCTopo - 9

vector<int> algoMap; // stores the algorithm used for inputTree[x] in the case that we are not using all the input types - split/filter, trim, prune and recluster
vector<int> fileMap; // stores the index of the input file for each grooming algorithm

//DECLARATIONS FOR RECLUSTERING FUNCTIONS

vector<fastjet::PseudoJet> ObjsToPJ(vector<TLorentzVector> jets);
vector<TLorentzVector> Recluster(vector<TLorentzVector> small_jets, double PTcut=15, double fcut=0.05, double jetRad=1.0);

/////



TFile *inputFile[20];
TTree *inputTree[20];

//get the branches we want to use

//QCD split filtering with Y cut 9
vector<float> * qcd_CA12_truth_pt = 0;
vector<float> * qcd_CA12_truth_eta = 0;
vector<float> * qcd_CA12_truth_phi = 0;
vector<float> * qcd_CA12_truth_mass = 0;
vector<float> * qcd_CA12_truth_E = 0;

vector<float> * qcd_CA12_topo_pt = 0;
vector<float> * qcd_CA12_topo_eta = 0;
vector<float> * qcd_CA12_topo_phi = 0;
vector<float> * qcd_CA12_topo_mass = 0;
vector<float> * qcd_CA12_topo_E = 0;
vector<float> * qcd_CA12_topo_emfrac = 0;

vector<float> * qcd_CA12_groomed_pt = 0;
vector<float> * qcd_CA12_groomed_eta = 0;
vector<float> * qcd_CA12_groomed_phi = 0;
vector<float> * qcd_CA12_groomed_mass = 0;
vector<float> * qcd_CA12_groomed_E = 0;
vector<float> * qcd_CA12_groomed_emfrac = 0;

UInt_t qcd_mc_channel_number = 0; 
Float_t qcd_mc_event_weight = 0; 



//Wprime split filtering both cuts

vector<float> * Wp_CA12_truth_pt = 0;
vector<float> * Wp_CA12_truth_eta = 0;
vector<float> * Wp_CA12_truth_phi = 0;
vector<float> * Wp_CA12_truth_mass = 0;
vector<float> * Wp_CA12_truth_E = 0;
UInt_t Wp_mc_channel_number = 0; 
Float_t Wp_mc_event_weight = 0; 

vector<float> * Wp_CA12_topo_pt = 0;
vector<float> * Wp_CA12_topo_eta = 0;
vector<float> * Wp_CA12_topo_phi = 0;
vector<float> * Wp_CA12_topo_mass = 0;
vector<float> * Wp_CA12_topo_E = 0;
vector<float> * Wp_CA12_topo_emfrac = 0;

vector<float> * Wp_CA12_groomed_pt = 0;
vector<float> * Wp_CA12_groomed_eta = 0;
vector<float> * Wp_CA12_groomed_phi = 0;
vector<float> * Wp_CA12_groomed_mass = 0;
vector<float> * Wp_CA12_groomed_E = 0;
vector<float> * Wp_CA12_groomed_emfrac = 0;

vector<float> * currTree_mass = 0;

// create a hashtable/ map for all the variables we want to plot.... Index on variable name - we can read this in from xml config

//counters for efficiency issues
int nEvt_0=0; //total events
int nEvt_1=0; // events with leading proper topo jet
int nEvt_2=0; // events with matched truth to proper topo jet
int nEvt_3=0; // events with split cut 9 matched to truth jet matches to proper topo jet
int nEvt_4=0;

//counters for efficiency issues
int nEvt1_0=0; //total events
int nEvt1_1=0; // events with leading proper topo jet
int nEvt1_2=0; // events with matched truth to proper topo jet
int nEvt1_3=0; // events with split cut 9 matched to truth jet matches to proper topo jet
int nEvt1_4=0;


const int nAlgosMax=12;
int nAlgos=0;
TString AlgoList[nAlgosMax];
TString binLabel[nAlgosMax-2];
const int nPtBins=6;
TString pTbins[nPtBins];
const int nFineBins=12;
TString finePtBins[nFineBins];


TH1F *qcd_finePtBin_mass[nAlgosMax][nFineBins];
TH1F *Wprime_finePtBin_mass[nAlgosMax][nFineBins];

TH1F *qcd_Lead_CA12_mass[nAlgosMax][nPtBins]; 
TH1F *Wprime_Lead_CA12_mass[nAlgosMax][nPtBins];
TH1F *qcd_Lead_CA12_pt[nAlgosMax];
TH1F *Wprime_Lead_CA12_pt[nAlgosMax];
TH1F *Wprime_Lead_CA12_scaled_pt[nAlgosMax];
TH1F *pTweights[nAlgosMax];
TH1F *qcd_PtReweight[nAlgosMax];
TH1F *Wp_PtReweight[nAlgosMax];
TString AlgoListN[nAlgosMax];
TString pTbinsN[nPtBins];
TString pTFinebinsN[nFineBins];

double myMPV[nAlgosMax][nPtBins];
double WidthMassWindow[nAlgosMax][nPtBins];
double TopEdgeMassWindow[nAlgosMax][nPtBins];
double BottomEdgeMassWindow[nAlgosMax][nPtBins];
double QCDfrac[nAlgosMax][nPtBins];

double myMPV_finePt[nAlgosMax][nFineBins];
double WidthMassWindow_finePt[nAlgosMax][nFineBins];
double TopEdgeMassWindow_finePt[nAlgosMax][nFineBins];
double BottomEdgeMassWindow_finePt[nAlgosMax][nFineBins];
double QCDfrac_finePt[nAlgosMax][nFineBins];


TH1F * hMassLow[nPtBins];
TH1F * hMassHigh[nPtBins];
TH1F * hMPV[nPtBins];
TH1F * hWmass[nPtBins];
TH1F * hQCDeff[nPtBins];

TH1F * hQCDeff_finePt[nAlgosMax];


TH1F * windowsVsPt[nAlgosMax];

// ROC
TGraph *finePtBin_mass_curve[nAlgosMax][nFineBins];
TGraph *Lead_CA12_mass_curve[nAlgosMax][nPtBins];
TGraph *Lead_CA12_pt_curve[nAlgosMax];
TGraph *Lead_CA12_scaled_pt_curve[nAlgosMax];
TGraph *pTweights_curve[nAlgosMax];
TGraph *PtReweight_curve[nAlgosMax];


TCanvas * c1[nAlgosMax][nPtBins]; 
TCanvas * c3[nAlgosMax][nFineBins];
TCanvas * c2[nPtBins];
TPad *pad1[nPtBins];
TPad *pad2[nPtBins];


void defineStrings(TString *AlgoList, TString *binLabel, TString *pTbins, TString *finePtBins){
  
  AlgoList[0]="TruthJet_RecoMatch";
  //SplitFilteredMu67SmallR0YCut9
  AlgoList[1]="SF67r0Y9";
  //SplitFilteredMu100SmallR30YCut4
  AlgoList[2]="SF100r30Y4";
  //TrimmedPtFrac5SmallR30
  AlgoList[3]="TrimPt5r30";
  //TrimmedPtFrac5SmallR20
  AlgoList[4]="TrimPt5r20";
  //PrunedCaRcutFactor50Zcut10
  AlgoList[5]="PrunRf50Z10";
  //PrunedCaRcutFactor50Zcut20
  AlgoList[6]="PrunRf50Z20";
  //reclustering
  AlgoList[7]="ReclusAK2";
  AlgoList[8]="ReclusAK3";
  AlgoList[9]="ReclusAK4";
  //
  AlgoList[10]="LeadTruthJet"; //for pt reweighting
  AlgoList[11]="LeadTruthJet_finebin"; //to show off
    
  binLabel[0]="CA12 Truth RecoMatch";
  binLabel[1]="CA12 Split-Filt (0.67, 0.09, variable)";
  binLabel[2]="CA12 Split-Filt (1.0, 0.04, 0.3)";
  binLabel[3]="CA12 Trimmed (0.05, 0.3)";
  binLabel[4]="CA12 Trimmed (0.05, 0.2)";
  binLabel[5]="CA12 Pruned (0.5, 0.1)";
  binLabel[6]="CA12 Pruned (0.5, 0.2)";
  binLabel[7]="AK10 Reclustered (from AK2)";
  binLabel[8]="AK10 Reclustered (from AK3)";
  binLabel[9]="AK10 Reclustered (from AK4)";

  pTbins[0]="inclusive";
  pTbins[1]="gtr1000";
  pTbins[2]="500to1000";
  pTbins[3]="1000to1500";
  pTbins[4]="1500to2000";
  pTbins[5]="gtr2000";
  
  finePtBins[0]="0to250";
  finePtBins[1]="250to500";
  finePtBins[2]="500to750";
  finePtBins[3]="750to1000";
  finePtBins[4]="1000to1250";
  finePtBins[5]="1250to1500";
  finePtBins[6]="1500to1750";
  finePtBins[7]="1750to2000";
  finePtBins[8]="2000to2250";
  finePtBins[9]="2250to2500";
  finePtBins[10]="2500to2750";
  finePtBins[11]="2750to3000";
  
  AlgoListN[0]="Truth CA12 matched to a reco ungroomed jet";
  //SplitFilteredMu67SmallR0YCut9
  AlgoListN[1]="CA12 Split Filtered jets, #mu=0.67, y_{cut}=0.09, R_{subjet} variable";
  //SplitFilteredMu100SmallR30YCut4
  AlgoListN[2]="CA12 Split Filtered jets, #mu=1.0, y_{cut}=0.04, R_{subjet}=0.3";
  //TrimmedPtFrac5SmallR30
  AlgoListN[3]="CA12 Trimmed jets, f_{cut}=0.05, R_{subjet}=0.3";
  //TrimmedPtFrac5SmallR20
  AlgoListN[4]="CA12 Trimmed jets, f_{cut}=0.05, R_{subjet}=0.2";
  //PrunedCaRcutFactor50Zcut10
  AlgoListN[5]="CA12 Pruned jets, R_{cut}=0.5, Z_{cut}=0.1";
  //PrunedCaRcutFactor50Zcut20
  AlgoListN[6]="CA12 Pruned jets, R_{cut}=0.5, Z_{cut}=0.2";
  AlgoListN[7]="Reclustering AntiKt10 from AntiKt2 jets";
  AlgoListN[8]="Reclustering AntiKt10 from AntiKt3 jets";
  AlgoListN[9]="Reclustering AntiKt10 from AntiKt4 jets";

  AlgoListN[10]="Leading truth CA12 jet (unmatched)"; //for pt reweighting
  AlgoListN[11]="Leading truth CA12 jet (unmatched)"; //to show off

  pTbinsN[0]="Inclusive p_{T}^{CA12}";
  pTbinsN[1]="p_{T}^{CA12} > 1000 GeV";
  pTbinsN[2]="500 < p_{T}^{CA12} < 1000 GeV";
  pTbinsN[3]="1000 < p_{T}^{CA12} < 1500 GeV";
  pTbinsN[4]="1500 < p_{T}^{CA12} < 2000 GeV";
  pTbinsN[5]="p_{T}^{CA12} > 2000 GeV";
  
  pTFinebinsN[0]="0 < p_{T}^{CA12} < 250 GeV";
  pTFinebinsN[1]="250 < p_{T}^{CA12} < 500 GeV";
  pTFinebinsN[2]="500 < p_{T}^{CA12} < 750 GeV";
  pTFinebinsN[3]="750 < p_{T}^{CA12} < 1000 GeV";
  pTFinebinsN[4]="1000 < p_{T}^{CA12} < 1250 GeV";
  pTFinebinsN[5]="1250 < p_{T}^{CA12} < 1500 GeV";
  pTFinebinsN[6]="1500 < p_{T}^{CA12} < 1750 GeV";
  pTFinebinsN[7]="1750 < p_{T}^{CA12} < 2000 GeV";
  pTFinebinsN[8]="2000 < p_{T}^{CA12} < 2250 GeV";
  pTFinebinsN[9]="2250 < p_{T}^{CA12} < 2500 GeV";
  pTFinebinsN[10]="2500 < p_{T}^{CA12} < 2750 GeV";
  pTFinebinsN[11]="2750 < p_{T}^{CA12} < 3000 GeV";



}


static inline std::string &rtrim(std::string &s) {
  s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
  return s;
}
