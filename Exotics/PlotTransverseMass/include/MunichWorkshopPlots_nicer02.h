
#include <TH1F.h>
#include <TH2F.h>
#include "TChain.h"
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
#include <unordered_map>
#include <cmath>

#include "TROOT.h"
#include "TError.h"
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



//Declarations of alllllll the methods. check cpp file for details of methods.

float DeltaR(float eta1,float phi1,float eta2,float phi2);
double mpv(TH1F* histo);
void Qw(double &minWidth, double &topEdge, TH1F* histo, double frac);

void makeROC(int type, TH1F *&S,TH1F *&B,TGraph &curve, TString name="", bool draw=false);
std::vector<std::string> getListOfBranches(std::string &algorithm);
void makeMassWindowFile(bool applyMassWindow, std::string & algorithm);
void setJetsBranches(TChain * tree, std::string & algorithm, std::string & groomIdx, std::unordered_map<std::string,bool> & current_branchmap);
void addInfoBranches(TTree * tree);
void addSubJets(TTree * tree, std::string & algorithm, std::string & groomIdx);
vector<std::pair<std::string,bool> > getListOfJetBranches(std::string & algorithm, std::unordered_map<std::string, bool> & brancharray);
void initVectors();
std::pair<int,int> getTwoLeadingSubjets(std::vector<int> & indices, std::vector<float> *& subjets);
std::string returnJetType(std::string & samplePrefix, std::string & groomalgo, bool addLC, int idx, bool underscore = true);
std::string returnSubJetType(std::string & samplePrefix, std::string & groomalgo, bool addLC);
void clearVectors();
void clearOutputVariables();
void setOutputVariables(int idx1, int idx2, int idx3, int idxca, int idxcatopo, int subidx, std::string & algorithm, std::string & groomalgo, std::string & groomIdx);
void setOutputBranches(TTree* tree, std::string & algorithm, std::string & groomIdx);
void resetOutputVariables();
void setAddress(TChain * tree, std::string  name, std::vector<Float_t> * var_vec);
void readWeights();
void createPtReweightFile(TH1F * bkg, TH1F * sig, std::string & fname);
void overlapRemoval();
int eventSelection();
bool leptonSelection(int lepType);
void setLeptons(TChain * tree, TObjArray * list);
void addLeptonBranches(std::string & jetString, TChain * tree);
void setLeptonVectors();
std::vector<float> dummyCharge(int size);
bool setVector(TChain *& tree, std::unordered_map<std::string, int> & list, vector<TLorentzVector> *& vec, std::string branch);

bool setVector(TChain *& tree, std::unordered_map<std::string, int> & list, vector<Float_t> *& vec, std::string branch);
bool setVector(TChain *& tree, std::unordered_map<std::string, int> & list, vector<Int_t> *& vec, std::string branch);
bool useBranch(std::string const & branch, bool partialmatch = false);
void setLLJMass(int jetidx);
double calculateFoxWolfram20(vector<TLorentzVector>& clusters);
int calculateSoftDropTag(vector<TLorentzVector> & cluster);
vector<double> calculateQJetsVol(vector<TLorentzVector> & EventParticles, float radius);
double calculateQJetsVol_v2(vector<TLorentzVector> & clusters);
bool createClusters(int jettype, int jetidx, vector<TLorentzVector> & cluster);
double calculateEEC(int jettype, float beta=0.3, float exp=2);
void calculateECF(vector<TLorentzVector> & clusters, int jettype, float beta=0.3);
void setEEC(int jettype, int jetidx);
void setRadius(std::string & prefix);
void printTLV(vector<TLorentzVector> & tlv);
std::unordered_map<std::string,bool> createBranchMap(TObjArray *& arr);
void SignalHandlerMapAccess(int signal);
void calculateResponseValues();
void setCa12Vectors(bool truth, bool topo);

// typedef for the exception when accessing a missing element from a map
typedef void (*SignalHandlerPointer)(int);

// return dummy vectors
void tlvvec(vector<TLorentzVector> *& tmp);
void tlvvecvec(vector<vector<TLorentzVector> > *& tmp);
void floatvec(vector<float> * & tmp);
void intvec(vector<int> * & tmp);
void vecintvec(vector< vector<int> > *& tmp);


bool replace(std::string& str, const std::string& from, const std::string& to) {
  size_t start_pos = str.find(from);
  if(start_pos == std::string::npos)
    return false;
  str.replace(start_pos, from.length(), to);
  return true;
}

// calculate the mean of a vector of masses
inline double mean(vector<double>& masses){
  double ret(0.);
  for(vector<double>::iterator it = masses.begin(); it != masses.end(); it++)
    ret += (*it);
  return ret/masses.size();
}

// calculate the RMS of a vector of masses. 
inline double var(vector<double>& masses){
  double ret(0.), avg(mean(masses));
  for(vector<double>::iterator it = masses.begin(); it != masses.end(); it++)
    ret += ((*it)-avg)*((*it)-avg);
  return ret/masses.size();
}

// overloaded version where the avg is passed as well and doesn't need to be recalculated
inline double var(vector<double>& masses, double avg){
  double ret(0.);
  // average is mean(masses)
  for(vector<double>::iterator it = masses.begin(); it != masses.end(); it++)
    ret += ((*it)-avg)*((*it)-avg);
  return ret/masses.size();
}


// enums used for sample types, jet types, histogram types and lepton types
enum sampleType{BACKGROUND,SIGNAL,DATA};
enum jetType{TRUTH,TOPO,GROOMED,MAX};
enum leptonType{FAIL,ELECTRON,MUON};

//DECLARATIONS FOR RECLUSTERING FUNCTIONS
vector<fastjet::PseudoJet> ObjsToPJ(vector<TLorentzVector> jets);
vector<TLorentzVector> Recluster(vector<TLorentzVector> small_jets, double PTcut=15, double fcut=0.05, double jetRad=1.0);

// struct that is used for reading in the xml containing the algorithm information
struct Algorithms
{
  std::unordered_map<std::string, std::string> AlgoNames; // the jet groomieng alg
  std::unordered_map<std::string, std::string> AlgoPrefix; // the jet reco algorithm
  std::unordered_map<std::string, std::string> AlgoType; // if split/filter, trimmed, etc
  std::unordered_map<std::string, std::string> AlgoList; // the abbreviated name
  std::unordered_map<std::string, std::string> AlgoListN; // the plotting labels
  std::unordered_map<std::string, std::string> subjetMap; // subjet grooming alg for jet grooming alg 
  std::unordered_map<std::string, std::string> subjetIndex;
  std::unordered_map<std::string, std::string> binLabel;
  void load(const std::string & filename);
};

// create one of the algorithm structs
struct Algorithms algorithms;
/////

// how many events pass selection
long passed_counter = 0;

// settings for histograms
const int nAlgosMax=12;
int nAlgos=0;

const int nPtBins=6;
TString pTbins[nPtBins];
const int nFineBins=12;
TString finePtBins[nFineBins];
std::vector<pair<float,float> > ptrange;
const float ptweightBins[25] = {200,250,300,350,400,450,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1800,2000,2200,2400,2600,2800,3000};


// defines
float GEV = 1000.;
float ELMASS = 0.511;
float MUMASS = 105.7;

// global variables used per analysis
std::string fileid_global;
std::string treeName;
std::string weightsfile;
std::string branchesFile;
std::unordered_map<std::string, bool> branchmap;

// settings for selection
bool calcQJets = false;
bool calcFoxWolfram20 = false;
bool preCalcFoxWolfram20 = false;
bool calcSoftDrop = false;
bool calcEEC = false;
bool preCalcEEC = false;
bool calcClusters = false;
bool calcTauWTA21 = true;
bool recluster = false;
bool calcYFilt = false;
bool truthBosonMatching = false;
bool truthBoson4vec = false;
bool beta2available = false;
bool addResponse = false;
bool clusterTLV = false;
bool ca12TLV = false;
bool ca12topoTLV = false;
bool keepTopo = true;

bool xAODJets = false;
//bool xAODemfrac = false;
bool hvtllqq = false;
bool xAOD = false;

float radius = 1.0;
int nqjets = 25;



TChain *inputTChain[2];

// histograms
TH1F * pt_reweight_arr[2];
TH2F * cluster_vs_truthpt = 0;

// variables in every tree
Float_t normalisation = 1.0;
Int_t NEvents = 0;
std::map<long, float> NEvents_weighted;
// Change this to Float_t for D3PD, vector<float>* for xaod
Float_t mc_event_weight_d3pd = 1.0;
std::vector<float> * mc_event_weight_xaod = 0;
Float_t mc_event_weight_out = 1.0;
UInt_t mc_channel_number  = 0;
UInt_t mc_channel_number_out  = 0;
UInt_t runNumberOut = 0;
UInt_t EventNumberOut = 0;
UInt_t runNumberIn = 0;
UInt_t EventNumber = 0;
Float_t avgIntpXingOut = 0;

// Change this to Float_t for D3PD, UInt_t for xaod
UInt_t avgIntpXingIn_xaod = 0; // new xaods might need this to be float
Float_t avgIntpXingIn_d3pd = 0;

Float_t actualIntPerXingIn = 0;
Float_t actualIntPerXingOut = 0;
// xaod expects UInt_t, but D3PDs Int_t. This is annoying and will cause problems
UInt_t nvtxIn = 0;
UInt_t nvtxOut = 0;

// settings for algorithm

// values calculated for signal and background mass windows
TString AlgoListN;
TString pTbinsN[nPtBins];
TString pTFinebinsN[nFineBins];

// define strings for a bunch of the pt ranges and pt bins
void defineStrings(TString *pTbins, TString *finePtBins){

  pTbins[0]="inclusive";
  pTbins[1]="gtr1000";
  pTbins[2]="500to1000";
  pTbins[3]="1000to1500";
  pTbins[4]="1500to2000";
  pTbins[5]="gtr2000";
  ptrange.push_back(std::make_pair(0,10000));
  ptrange.push_back(std::make_pair(1000,10000));
  ptrange.push_back(std::make_pair(500,1000));
  ptrange.push_back(std::make_pair(1000,1500));
  ptrange.push_back(std::make_pair(1500,2000));
  ptrange.push_back(std::make_pair(2000,10000));
  
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

// trim whitespace from rhs of string
static inline std::string &rtrim(std::string &s) {
  s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
  return s;
}


//store all the weights for different runNumbers:
std::unordered_map<long, float> k_factors;
std::unordered_map<long, float> filt_eff;
std::unordered_map<long, float> xs;
std::unordered_map<long, int> runNumber_map;

// variables used for reading in from tree.
// the double vector is to store the variable for truth, topo and groomed
std::vector<std::vector<float> *> var_E_vec;
std::vector<std::vector<float> *> var_pt_vec;
std::vector<std::vector<float> *> var_m_vec;
std::vector<std::vector<float> *> var_eta_vec;
std::vector<std::vector<float> *> var_phi_vec;
std::vector<std::vector<float> *> var_emfrac_vec;
std::vector<std::vector<float> *> var_Tau1_vec;
std::vector<std::vector<float> *> var_Tau2_vec;
std::vector<std::vector<float> *> var_SPLIT12_vec;
std::vector<std::vector<float> *> var_Dip12_vec;
std::vector<std::vector<float> *> var_PlanarFlow_vec;
std::vector<std::vector<float> *> var_Angularity_vec;
std::vector<std::vector<float> *> var_TauWTA1_vec; 
std::vector<std::vector<float> *> var_TauWTA2_vec; 
std::vector<std::vector<float> *> var_TauWTA2TauWTA1_vec; 
std::vector<std::vector<float> *> var_ZCUT12_vec;
std::vector<std::vector<float> *> var_Aplanarity_vec;
std::vector<std::vector<float> *> var_Sphericity_vec;
std::vector<std::vector<float> *> var_ThrustMaj_vec;
std::vector<std::vector<float> *> var_ThrustMin_vec;
std::vector<std::vector<float> *> var_FoxWolfram0_vec;
std::vector<std::vector<float> *> var_FoxWolfram2_vec;
std::vector<std::vector<float> *> var_FoxWolfram20_vec;
std::vector<std::vector<int> *> var_SoftDropTag_vec;
std::vector<std::vector<float> *> var_ECF1_vec;
std::vector<std::vector<float> *> var_ECF2_vec;
std::vector<std::vector<float> *> var_ECF3_vec;
std::vector<std::vector<float> *> var_ECF1_2_vec;
std::vector<std::vector<float> *> var_ECF2_2_vec;
std::vector<std::vector<float> *> var_ECF3_2_vec;
std::vector<std::vector<float> *> var_Mu12_vec;


std::vector<TLorentzVector> * var_ca12_tlv_vec;
std::vector<float> * var_ca12_pt_vec;
std::vector<float> * var_ca12_phi_vec;
std::vector<float> * var_ca12_eta_vec;
std::vector<float> * var_ca12_m_vec;
std::vector<TLorentzVector> * var_ca12topo_tlv_vec;
std::vector<float> * var_ca12topo_pt_vec;
std::vector<float> * var_ca12topo_phi_vec;
std::vector<float> * var_ca12topo_eta_vec;
std::vector<float> * var_ca12topo_m_vec;
std::vector<float> * var_YFilt_vec;

// reading in truth boson info when running on xAOD
std::vector<float> * var_truthboson_pt_vec;
std::vector<float> * var_truthboson_eta_vec;
std::vector<float> * var_truthboson_phi_vec;
std::vector<int> * var_truthboson_ID_vec;
// the new mc15 samples have truthboson info in a vector<TLV>
std::vector<TLorentzVector> * var_truthboson_tlv_vec;

// reading in jet clusters
Int_t var_cl_n;
std::vector<float> * var_cl_pt_vec;
std::vector<float> * var_cl_eta_vec;
std::vector<float> * var_cl_phi_vec;

// extra variables to read in

std::vector<int> * vxp_nTracks;
Float_t var_massFraction_vec;
Float_t var_ktycut2_vec;



// electrons in
vector<TLorentzVector> * var_electrons_vec;
std::vector<TLorentzVector> electrons;
std::vector<float> * var_electronPt_vec;
std::vector<float> * var_electronEta_vec;
std::vector<float> * var_electronPhi_vec;
std::vector<float> * var_el_ptcone20_vec;
std::vector<float> * var_el_etcone20_vec;

// muons in
std::vector<TLorentzVector> * var_muons_vec;
std::vector<TLorentzVector> muons;
std::vector<float> * var_muonPt_vec;
std::vector<float> * var_muonEta_vec;
std::vector<float> * var_muonPhi_vec;
std::vector<float> * var_mu_ptcone20_vec;
std::vector<float> * var_mu_etcone20_vec;
std::vector<float> * var_mu_charge_vec;

// jet constituents
std::vector<std::vector<std::vector < int>  > * > var_constit_index;
std::vector< std::vector < int>  * > var_constit_n;

// subjet variables
std::vector<std::vector <int> > * subjet_index;
std::vector<float> * var_subjets_E_vec;
std::vector<float> * var_subjets_pt_vec;
std::vector<float> * var_subjets_m_vec;
std::vector<float> * var_subjets_eta_vec;
std::vector<float> * var_subjets_phi_vec;

// subjet clusters

// variables that we write out to the outfile
// store one value each for truth, topo, groomed - use a vector
std::vector<float> var_E;
std::vector<float> var_pt;
std::vector<float> var_m;
std::vector<float> var_eta;
std::vector<float> var_phi;
std::vector<float> var_emfrac;
std::vector<float> var_Tau1;
std::vector<float> var_Tau2;
std::vector<float> var_SPLIT12;
std::vector<float> var_Dip12;
std::vector<float> var_PlanarFlow;
std::vector<float> var_Angularity;
std::vector<float> var_Tau21;
std::vector<float> var_Mu12;
Float_t var_ca12_pt;
Float_t var_ca12_eta;
Float_t var_ca12_phi;
Float_t var_ca12_m;
Float_t var_ca12topo_pt;
Float_t var_ca12topo_eta;
Float_t var_ca12topo_phi;
Float_t var_ca12topo_m;
Float_t var_YFilt; // only for groomed

// leptons out
std::vector<TLorentzVector> var_leptons;
std::vector<Float_t> var_ptcone20;
std::vector<Float_t> var_etcone20;
Float_t var_mllj;
Float_t var_mll;
Float_t var_ptll;

// store leading jet pt, before any truth matching or selections
Float_t var_leadingJetPt;

// if it is an electron or muon event
Int_t var_isElectronEvent = false;

std::vector<float> var_charge;

// these variables are only stored for the subjets of the groomed jets, so we don't need a vector
Float_t var_subjets_E;
Float_t var_subjets_pt;
Float_t var_subjets_m;
Float_t var_subjets_eta;
Float_t var_subjets_phi;
Float_t var_massdrop;
Float_t var_yt;

// extra output variables
std::vector<Float_t> var_TauWTA1; 
std::vector<Float_t> var_TauWTA2; 
std::vector<Float_t> var_TauWTA2TauWTA1; 
std::vector<Float_t> var_ZCUT12;

std::vector<Float_t> var_Aplanarity;
std::vector<Float_t> var_Sphericity;
std::vector<Float_t> var_ThrustMaj;
std::vector<Float_t> var_ThrustMin;

std::vector<Float_t> var_QjetVol;
std::vector<Float_t> var_FoxWolfram20;
std::vector<Int_t> var_softdrop;
std::vector<Float_t> var_EEC_C2_1;
std::vector<Float_t> var_EEC_C2_2;
std::vector<Float_t> var_EEC_D2_1;
std::vector<Float_t> var_EEC_D2_2;


// read in clusters from xAODs
std::vector<std::vector<TLorentzVector> > * var_clusters_truth_vec;
std::vector<std::vector<TLorentzVector> > * var_subjets_truth_vec;
std::vector<std::vector<TLorentzVector> > * var_clusters_groomed_vec;
std::vector<std::vector<TLorentzVector> > * var_subjets_groomed_vec;
std::vector<std::vector<TLorentzVector> > * var_clusters_ca12_vec;
std::vector<std::vector<TLorentzVector> > * var_subjets_ca12_vec;

// store clusters
std::vector<TLorentzVector> var_clusters_truth;
std::vector<TLorentzVector> var_subjets_truth;
std::vector<TLorentzVector> var_clusters_groomed;
std::vector<TLorentzVector> var_subjets_groomed;
std::vector<TLorentzVector> var_clusters_ca12;
std::vector<TLorentzVector> var_subjets_ca12;

// store the weights for the samples
Float_t var_k_factor;
Float_t var_filter_eff;
Float_t var_xs;
// some have the weights pre-calculated
Float_t evt_kfactor = 0;
//Long_t evt_nEvts = 0; // this is on the dc14 samples
Float_t evt_nEvts = 0; // this is on the mc15 samples
Float_t evt_filtereff = 0;
Float_t evt_sumWeights = 0;
Float_t evt_xsec = 0;

Float_t evt_kfactor_out;
Float_t evt_nEvts_out;
Float_t evt_filtereff_out;
Float_t evt_sumWeights_out;
Float_t evt_xsec_out;
// for reading in
Float_t var_filtereff_in;
Float_t var_kfactor_in;
Float_t scale1fb = 0;
Float_t scale1fbOut;
// if subjets are being used
bool subjetscalc;
bool subjetspre;

// response values
Float_t response_E;
Float_t response_pt;
Float_t response_m;
Float_t response_eta;
Float_t response_phi;
Float_t response_Tau1;
Float_t response_Tau2;
Float_t response_SPLIT12;
Float_t response_Dip12;
Float_t response_PlanarFlow;
Float_t response_Angularity;
Float_t response_Tau21;
Float_t response_Mu12;

// extra output responseiables
Float_t response_TauWTA1; 
Float_t response_TauWTA2; 
Float_t response_TauWTA2TauWTA1;
Float_t response_ZCUT12;

Float_t response_Aplanarity;
Float_t response_Sphericity;
Float_t response_ThrustMaj;
Float_t response_ThrustMin;

Float_t response_QjetVol; //
Float_t response_FoxWolfram20;
Int_t response_softdrop;
Float_t response_EEC_C2_1; 
Float_t response_EEC_C2_2; 
Float_t response_EEC_D2_1; 
Float_t response_EEC_D2_2; 

Int_t vxp_nTracks_out;
