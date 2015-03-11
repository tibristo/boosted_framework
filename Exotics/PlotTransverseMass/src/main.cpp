///
///T. Bristow and F. A. Dias (flavia.dias@cern.ch)
///Analysis code for Jet Substructure 
///
///
///v0.1 - August 14th, 2014
///v0.2 - August 25th, 2014
///v0.3 - August 26th, 2014
///v0.4 - August 30th, 2014
///
///v whatevs 
///

#include "MunichWorkshopPlots_nicer02.h"
#include <string>
#include <utility>
#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/regex.hpp>
#include <algorithm>
#include "QjetsPlugin.h"
#include <signal.h>


//#include "LinkDef.h"
using namespace boost::algorithm;
using namespace std;
namespace po = boost::program_options;

// used for controlling cout statements for debugging
#define DEBUG 0


template<class T>
ostream& operator<<(ostream& os, const vector<T>& v)
{
  copy(v.begin(), v.end(), ostream_iterator<T>(os, " ")); 
  return os;
  }

int main( int argc, char * argv[] ) {

  //gROOT->ProcessLine("#include <vector>");
  // so that we can use vector<TLorentzVector> in ROOT
  gROOT->LoadMacro( "include/TLorentzVectorDict.h+" );
  AtlasStyle();
  
  bool makePtPlotsFlag = false;
  bool scaleHistsFlag = false;  
  bool makePlotsFlag = false;
  bool getMPVFlag = false;
  bool checkBkgFrac = false;
  bool applyMassWindowFlag = false;
  std::string fileid = "";
  po::variables_map vm;
  // set up the argument management
  try {
    int opt;
    vector<string> config_file;
    
    // Declare a group of options that will be 
    // allowed only on command line
    po::options_description generic("Generic options");
    generic.add_options()
      ("help", "produce help message")
      ("config,c", po::value<vector<string> >()->multitoken(),
       "name of a file of a configuration.")
      ;
    
    // Declare a group of options that will be 
    // allowed both on command line and in
    // config file
    std::cout << " setting configuration" << std::endl;
    po::options_description config("Configuration");
    config.add_options()
      ("signal-file", po::value<vector<string> >()->multitoken(), 
       "signal input files")
      ("background-file", po::value<vector<string> >()->multitoken(), 
       "background input files")
      ("algorithm", 
       po::value< string>(), 
       "algorithm type")
      ("subjets-calc", po::value<bool>(&subjetscalc)->default_value(false),"calculate mass drop and momentum balance from subjet kinematics")
      ("subjets-pre", po::value<bool>(&subjetspre)->default_value(false),"use massFraction and ktycut2 from the input files")

      ("mass-window", po::value<bool>(&applyMassWindowFlag)->default_value(false),"apply mass window cuts")
      ("make-plots", po::value<bool>(&makePlotsFlag)->default_value(false),"create plots")
      ("make-ptplots", po::value<bool>(&makePtPlotsFlag)->default_value(false),"create pT plots")
      ("scale-hists", po::value<bool>(&scaleHistsFlag)->default_value(false),"scale mass/ pt histograms")
      ("mpv", po::value<bool>(&getMPVFlag)->default_value(false),"Find the MPV")
      ("fileid", po::value<string>(&fileid)->default_value(""),"Add an identifier to the output folder/ file names")
      ("bkg-frac", po::value<bool>(&checkBkgFrac)->default_value(false),"Check the background fraction in the signal")
      ("tree-name", po::value<string>(&treeName)->default_value("physics"),"Name of tree to be read in from input file")
      ("weights-file", po::value<string>(&weightsfile)->default_value("weightings_input.csv"),"File with weights info - xsec, filter eff, k-factors.")
      ("branches-file", po::value<string>(&branchesFile)->default_value(""),"Name of file containing branches, otherwise Alg_branches.txt is used. Becareful with this, because if the branches are not read in then when any of them are used later on a segfault will occur, so make sure there is one that it will use.")
      ("xAOD-jets", po::value<bool>(&xAODJets)->default_value(false),"Indicate if we are running over xAOD output and there is the word Jets appended to the algorithm name.")
      //("xAOD-emfrac", po::value<bool>(&xAODemfrac)->default_value(false),"Indicate if we are running over xAOD output and emfrac is not available.")
      ("hvtllqq-selection", po::value<bool>(&hvtllqq)->default_value(false),"Indicate if we are running the HVTllqq analysis selection.")
      ("calcQjets", po::value<bool>(&calcQJets)->default_value(false),"Indicate if we are calculating the Qjet volatility.  Note this is potentially quite CPU intensive.")
      ("calcFoxWolfram", po::value<bool>(&calcFoxWolfram20)->default_value(false),"Indicate if we are calculating FoxWolfram20.  preCalcFoxWolfram and calcFoxWolfram cannot both be true. Note this is potentially quite CPU intensive. ")
      ("preCalcFoxWolfram", po::value<bool>(&preCalcFoxWolfram20)->default_value(false),"Indicate if we are using pre-calculated FoxWolfram20. preCalcFoxWolfram and calcFoxWolfram cannot both be true.")
      ("calcSoftDrop", po::value<bool>(&calcSoftDrop)->default_value(false),"Indicate if we are calculating the soft drop tag.  Note this is potentially quite CPU intensive.")
      ("calcEEC", po::value<bool>(&calcEEC)->default_value(false),"Indicate if we are calculating EEC.  preCalcEEC and this cannot both be true.  Note this is potentially quite CPU intensive.")
      ("calcYFilt", po::value<bool>(&calcYFilt)->default_value(false),"Indicate if we are calculating YFilt.  This is done for split filtered samples already.  If YFilt is not in the samples already set this to true.")
      ("preCalcEEC", po::value<bool>(&preCalcEEC)->default_value(false),"Indicate if we are using pre-calculated EEC. calcEEC and this cannot both be true.")
      ("calcClusters", po::value<bool>(&calcClusters)->default_value(false),"Reconstruct TLVs of the topo clusters.  This is done automatically if doing qjets, foxwolfram, softdrop or EEC.")
      ("calcTauWTA21", po::value<bool>(&calcTauWTA21)->default_value(true),"Calculate TauWTA2/TauWTA1.")
      ("xAOD", po::value<bool>(&xAOD)->default_value(false),"Set if running on xAOD, otherwise D3PD is assumed.  This affects the way the mc_event_weight and avgIntPerXing variables are read, and if set incorrectly segfaults will occur.")
      ("truthBosonMatching", po::value<bool>(&truthBosonMatching)->default_value(false),"Set if running on xAOD.  This matches jets to their parent particles - W or Z.  This is needed since no emfrac variable exists in the xAODs.")
      ("ecf-beta2", po::value<bool>(&beta2available)->default_value(false),"Set if running on xAOD and the ECF (beta2) values are available. Not all xAODs have this variable.")
      
      ;
        
    po::options_description cmdline_options;
    cmdline_options.add(generic).add(config);//.add(hidden);

    po::options_description config_file_options;
    config_file_options.add(config);//.add(hidden);

    po::options_description visible("Allowed options");
    visible.add(generic).add(config);

    po::positional_options_description p;
    p.add("background-file", 1);
    p.add("signal-file", 1);
    p.add("algorithm",1);

    store(po::command_line_parser(argc, argv).
	  options(cmdline_options).positional(p).run(), vm);
    notify(vm);

    if (vm.count("help")) {
      cout << visible << "\n";
      return 0;
    }

    // need to check what happens if we don't actually add anything here!

    if (vm.count("config"))
      {
	config_file = vm["config"].as<vector<string> >();
	for (vector<std::string>::iterator it = config_file.begin(); it!= config_file.end(); it++)
	  {
	    ifstream ifs((*it).c_str());
	    if (!ifs)
	      {
		cout << "can not open config file: " << config_file << "\n";
		return 0;
	      }
	    else
	      {
		store(parse_config_file(ifs, config_file_options), vm);
		notify(vm);
	      }	    
		
	  }
      }    


    if (vm.count("signal-file"))
      {
	cout << "Signal files are: " 
	     << vm["signal-file"].as< vector<string> >() << "\n";
      }
    else
      {
	cout << "No signal files!" << std::endl;
	return 0;
      }

    if (vm.count("background-file"))
      {
	cout << "Background files are: " 
	     << vm["background-file"].as< vector<string> >() << "\n";
      }
    else
      {
	cout << "No background files!" << std::endl;
	return 0;
      }
	
    if (vm.count("algorithm"))
      {
	cout << "algorithm is: " 
	     << vm["algorithm"].as< string >() << "\n";
      }
    else
      {
	cout << "No algorithm specified!" << std::endl;
	return 0;
      }
    if (vm.count("preCalcFoxWolfram") && vm.count("calcFoxWolfram"))
    {
      // check that pre and calc are not BOTH zero
      bool pre = vm["preCalcFoxWolfram"].as<bool>();
      bool calc = vm["calcFoxWolfram"].as<bool>();
      if (pre && calc)
	{
	  cout << "Both pre-calculation and calculated FoxWolfram20 are set to true - change one and re-run." << endl;
	  return 0;
	}
    }

    if (vm.count("preCalcEEC") && vm.count("calcEEC"))
    {
      // check that pre and calc are not BOTH zero
      bool pre = vm["preCalcEEC"].as<bool>();
      bool calc = vm["calcEEC"].as<bool>();
      if (pre && calc)
	{
	  cout << "Both pre-calculation and calculated EEC are set to true - change one and re-run." << endl;
	  return 0;
	}
    }



  }
  catch(exception& e)
    {
      cout << e.what() << "\n";
      return 1;
    }    

  // read in algorithm config
  algorithms.load("algorithms.xml");

  // whether we want to create the mass histograms
  bool massHistos = applyMassWindowFlag;

  // file identifier that gets appended to output files.
  fileid_global = fileid;

  vector<string> inputBkgFiles = vm["background-file"].as<vector<string> > ();
  vector<string> inputSigFiles = vm["signal-file"].as<vector<string> > ();

  // shuffle the file order to increase efficiency when running multiple instances of this at once
  std::random_shuffle(inputBkgFiles.begin(), inputBkgFiles.end());
  std::random_shuffle(inputSigFiles.begin(), inputSigFiles.end());


  std::string alg_in = vm["algorithm"].as<string>();
  std::string alg_in_orig = vm["algorithm"].as<string>();

  // set up the TChains for the signal and background files.
  inputTChain[sampleType::BACKGROUND] = new TChain(treeName.c_str());
  for (vector<string>::iterator it = inputBkgFiles.begin(); it != inputBkgFiles.end(); it++)
    {
      inputTChain[sampleType::BACKGROUND]->Add((*it).c_str());
      std::cout << "bkg file added: " << (*it) << std::endl;
    }

  inputTChain[sampleType::SIGNAL] = new TChain(treeName.c_str());
  for (vector<string>::iterator it = inputSigFiles.begin(); it != inputSigFiles.end(); it++)
    {
      inputTChain[sampleType::SIGNAL]->Add((*it).c_str());
      std::cout << "sig file added: " << (*it) << std::endl;
    }

  // create a lower case version of the algorithm name
  std::transform(alg_in.begin(), alg_in.end(), alg_in.begin(), ::tolower);
  std::cout << alg_in_orig << std::endl;

  int nArg = 1; // TODO: this will need to be changed when moving from a system where only one algorithm is run at once

  // right now the idea of filteridx, Xidx is not fully complete... The code as it stands runs on one type of algorithm at once, but 
  // I would like it to run on a number of types.  This is why this is here.  It used to do this in fact, before switching stuff around

  // check if we have a configuration for this algorithm
  if (algorithms.AlgoNames.find(alg_in_orig) == algorithms.AlgoNames.end())
    {
      std::cout << "Algorithm not defined in config file!" << std::endl;
      return 1;
    }
 
  // read in the event weights
  readWeights();
  // set up some of the strings for the ptBins
  defineStrings(pTbins, finePtBins);

  // make an output file
  makeMassWindowFile(applyMassWindowFlag, alg_in_orig);

  return 0;
}


/*
 * Return the delta R between two objects in eta/phi space.
 *
 * @param eta1 Eta of first object
 * @param phi1 Phi of first object
 * @param eta2 Eta of second object
 * @param phi2 Phi of second object
 *
 * @return delta R between two objects.
 */
float DeltaR(float eta1,float phi1,float eta2,float phi2){ 
  // we could actually inline this function
  float deltaPhi = TMath::Abs(phi1-phi2);
  float deltaEta = eta1-eta2;
  if(deltaPhi > TMath::Pi())
    deltaPhi = TMath::TwoPi() - deltaPhi;
  return TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
}

/*
 * Calculate the most probable value in a histogram.
 *
 * @param histo Pointer to the TH1F.
 *
 * @return most probable value.
 */
double mpv(TH1F* histo){
  double maxc=0.0001;
  int bin_maxc=0;
  int Nbins = histo->GetNbinsX(); 
  for( int i=0; i!=Nbins; i++){
    double tc = histo->GetBinContent(i);
    if(tc > maxc){
      maxc=tc;
      bin_maxc=i;
    }
  }
  double mostprob = histo->GetBinCenter(bin_maxc);
  return mostprob;
}

/*
 *
 */
void Qw(double &minWidth, double &topEdge, TH1F* histo, double frac){
  
  minWidth = 100000.0;
  topEdge = 0.0;
  
  int Nbins = histo->GetNbinsX();
  double integral = histo->Integral();

  for(int i=0; i!=Nbins; i++){
    double tempFrac = 0;
    int imax = i;
    
    // loop through until the tempFrac is above the frac (68%) criteria,
    // but making sure not to go out of range.
    while(tempFrac<frac && imax != Nbins){
      //fraction in bin imax=0,1,2,...                                                                                       
      tempFrac+=histo->GetBinContent(imax)/integral;
      imax+=1;
    }                                          
    
    double width = histo->GetBinCenter(imax) - histo->GetBinCenter(i); 
    
    double top_edge = histo->GetBinCenter(imax);
    
    // by applying this we say that the window we have just calculate MUST have
    // at least 68%.
    if(tempFrac >= frac && width<minWidth){ 
      minWidth = width;
      topEdge = top_edge;
    }
  }  
}

//FOR RECLUSTERING


///=========================================                                                                                                                                                                
/// Function for adding a new collection of jets which are re-clustered out of smaller jets.                                                                                                                  
///                                                                                                                                
/// Requires the function ObjsToPJ (given below)
/// You can also easily save the small radius jets - easiest if your
/// framework has a nice way to assign pointers to jets (can use the PJ_Info in FJ)
///
///=========================================     

vector<fastjet::PseudoJet> ObjsToPJ(vector<TLorentzVector> jets){
  vector<fastjet::PseudoJet> out;
  for(unsigned int i = 0; i < jets.size(); i++){
    fastjet::PseudoJet newjet(jets[i].Px(), jets[i].Py(), jets[i].Pz(), jets[i].E());
    out.push_back(newjet);                                                                           
    //out.push_back(clus);
  }
  return out;
}

vector<TLorentzVector> Recluster(vector<TLorentzVector> small_jets, double PTcut, double fcut, double jetRad){
  vector<fastjet::PseudoJet> particles = ObjsToPJ(small_jets);
  fastjet::JetDefinition fJetDef(fastjet::antikt_algorithm, jetRad, fastjet::Best);
  fastjet::ClusterSequence hardClustSeq(particles, fJetDef);
  vector<fastjet::PseudoJet> StandardJets = fastjet::sorted_by_pt(hardClustSeq.inclusive_jets(PTcut));
  vector<TLorentzVector> RCjets;
  for (unsigned int i=0; i<StandardJets.size(); i++){
    TLorentzVector sub = TLorentzVector();   
    sub.SetPtEtaPhiE(StandardJets[i].pt(), StandardJets[i].eta(), StandardJets[i].phi(), StandardJets[i].e());
    vector<fastjet::PseudoJet> constituents = StandardJets[i].constituents();
    for(unsigned int iCons = 0; iCons < constituents.size(); iCons++){
      //Do something with the small radius jets.
    }
    RCjets.push_back(sub);
  } 

  vector<TLorentzVector> RTjets;   
  //Now for my trimming on the re-clustered jets 
  for (unsigned int i=0; i<StandardJets.size(); i++){ 
    TLorentzVector trimmedjet = TLorentzVector(); 

    int trimmedjet_subjets_n = 0;
    
    vector<fastjet::PseudoJet> constituents = StandardJets[i].constituents();
    for(unsigned int iCons = 0; iCons < constituents.size(); iCons++){
      TLorentzVector subjet = TLorentzVector();
      subjet.SetPtEtaPhiE(constituents[iCons].pt(), constituents[iCons].eta(), constituents[iCons].phi(), constituents[iCons].e());    
      if (subjet.Pt() > fcut*RCjets[i].Pt()){     
        trimmedjet+=subjet;
	trimmedjet_subjets_n++;
      }
    }
    //if (trimmedjet_subjets_n>1)
    RTjets.push_back(trimmedjet);
  }
  
  return RTjets;
}



/*
 * Create a ROC curve for signal and background histograms.  This is stored in a TGraph which is passed by reference.
 *
 * @param type 1 is signal vs bkg rej and 2 is the inverse of this.
 * @param S Input signal histogram.
 * @param B Input background histogram.
 * @param curve Reference to TGraph.  This will store the ROC.
 * @param name Name of ROC.
 * @param draw Flag indicating if the ROC should be drawn and saved to file.
 */
void makeROC(int type, TH1F *&S, TH1F *&B, TGraph &curve, TString name, bool draw){

  //Need S and B to have the same number of bins!
  //TH1D *curve = new TH1D("","",M,0,1); where M is however fine you want the ROC curve binning to be.
  //ROC(S,B,curve);

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
      c1->SaveAs(name+".png");
      delete c1;
    }
  //curve_up=gr_up;
  //curve_do=gr_do;

  return;

} //makeROC



/*
 * Create an ntuple storing variables related to the grooming algorithm being run.  Some selections are run here too - event selection, 
 * jet and lepton selection, pT cuts.
 *
 * @param applyMassWindow Flag if any mass window cuts should be applied.
 * @param algorithm The name of the algorithm being run.
 */
void makeMassWindowFile(bool applyMassWindow,std::string & algorithm)
{
  std::cout << "starting make mass window output files " << std::endl;
  // set inital mass window values
  double mass_max = 10000000.0;
  double mass_min = 0.0;
  // is this background or signal?
  bool signal = false;

  // need branches for specific algorithms because otherwise we end up with a million branches that we don't need
  vector<std::pair<std::string,bool> > branches;

  // set up the algorithm name
  std::string i = algorithm;
  std::cout << "Doing mass window plots for " << algorithms.AlgoNames[algorithm] << std::endl;
  recluster = algorithms.AlgoType[algorithm].find("recluster") == std::string::npos ? false: true;

  
  // set the algorithm name details so that the lookup is not done in the loop
  std::string groomAlgoIndex = algorithm;
  std::string prefix = algorithms.AlgoPrefix[algorithm];
  std::string algorithmName = algorithms.AlgoNames[algorithm];
  std::string algorithmType = algorithms.AlgoType[groomAlgoIndex];

  // set the radius for the jet algorithm
  setRadius(prefix);

  // get an initial list of branches from one of the input files
  TObjArray * brancharray = inputTChain[0]->GetListOfBranches();
  // create an unordered map from this to be used for lookups.
  // tobjarray.findobject() returns a tvirtualobject which is slow
  std::unordered_map<string, bool> brancharray_initial = createBranchMap(brancharray);
  
  // read in the list of branches we want to use
  if (branchesFile != "")
    {
      std::cout << "branchesFile: " << branchesFile << std::endl;
      branches = getListOfJetBranches(branchesFile, brancharray_initial);
    }
  else
    {
      std::cout << "branchesFile not defined" << std::endl;
      branches = getListOfJetBranches(algorithmName, brancharray_initial);
    }
  
  // loop through the different pt bins. j == 0 is the inclusive one
    for (int j=0; j<1; j++){//nPtBins; j++){
      // set signal to false, so running background
      signal = false;
      
      // set up a stringstream to use for file names that can easily have integers added to it
      std::stringstream ss;
      ss << algorithmName << "_" << pTbins[j];
      

      // loop through background and signal
      for (int k = 0; k < 2; k++)
	{
	  // if k = 1 we are running signal
	  signal = k == 1? false : true;
	  // Initialise all the vectors to.. something 
	  initVectors();

	  // set up the mass window - no longer using the top/bottomedge mass window hists
	  mass_max = 300*GEV;//TopEdgeMassWindow[j];
	  mass_min = 0;//BottomEdgeMassWindow[j];
	  
	  // set the tchain pointer to the signal or background sample
	  int tchainIdx = signal ? sampleType::SIGNAL : sampleType::BACKGROUND;
	  std::stringstream ss_fname; // store the name of the output file and include the i and j indices!
	  std::string bkg = signal ? "sig": "bkg";
	  // create the output file name
	  ss_fname << ss.str() << "_" << bkg << ".root";

	  // make sure output directory exists
	  boost::filesystem::path dir(std::string(algorithm+fileid_global));
	  boost::filesystem::create_directory(dir);

	  // create a ttree cache to speed up branch reads
	  inputTChain[tchainIdx]->SetCacheSize(20000000);

	  // the branches map was done using a branch array from the background sample. Instead of creating a union of the
	  // branch array from signal and background the check is just done here.  The reason is two fold - I don't want to 
	  // fight with ROOT and what is inside the TObjArray, and since we only run this once for the signal and once for the
	  // background the overhead is minimal.
	  // If using TObjArray::FindObject() it will return a virtual TObject*, which will then run a comparison if testing for
	  // existence.  A solution is to add all of the names of the objects in the tree to a map, which will have faster 
	  // lookup, because no comparisons are run.
	  brancharray = 0;
	  brancharray = inputTChain[tchainIdx]->GetListOfBranches();
	  // unordered map for the current tree - this changes on iteration over k and j
	  std::unordered_map<std::string,bool> current_branchmap = createBranchMap(brancharray);

	  // turn off the branches we're not interested in

	  gROOT->ProcessLine("gErrorIgnoreLevel = 5000;");
	  for (std::vector<std::pair<std::string, bool > >::iterator it = branches.begin(); it != branches.end(); it++)
	    {
	      UInt_t * found = new UInt_t(1);// = 0;
	      //*found = 1;
	      // if the branch is not in this sample just skip it
	      if (current_branchmap.find((*it).first) == current_branchmap.end())
		continue;
	      // if the branch is in the samlple but not being used, turn it off
	      if ((*it).second == 0)
		{
		  inputTChain[tchainIdx]->SetBranchStatus((*it).first.c_str(),0,found);
		  // if the branch is not found in the tchain, then set it up in the branches file.
		  {
		    if ((&found) != 0)
		      (*it).second = -1;
		  }
		}
	      // otherwise add it to the ttree cache
	      else if ((*it).second == 1)
		{
		  //std::cout << (*it).first << std::endl;
		  // add branches we are using to branch cache
		  inputTChain[tchainIdx]->AddBranchToCache((*it).first.c_str());
		}
	      
	    }
	  // load all files and get number of entries
	  long entries = (long)inputTChain[tchainIdx]->GetEntries();
	  //set all of the branches for the output tree for the jets	  
	  setJetsBranches(inputTChain[tchainIdx], algorithmName, algorithm, current_branchmap); 

	  // output file to store cluster info
	  ofstream clusterfile(string("clusterfile"+bkg+".csv"));
	  // output file
	  TFile * outfile = new TFile(std::string(algorithm+fileid_global+"/"+ss_fname.str()).c_str(),"RECREATE");   
	  // output tree
	  TTree * outTree = new TTree(treeName.c_str(),treeName.c_str());

	  // reset all of the output variables so that they have default values
	  resetOutputVariables();
	  // set up all of the branches for the output tree
	  setOutputBranches(outTree, algorithmName, algorithm);

	  // if we are going to add the subjet branches we set them up here
	  if (subjetscalc || subjetspre)
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
	    addSubJets(outTree, AlgoNames[i], signal, i);
	  addJetsBranches(outTree, AlgoNames[i], signal, i);
<<<<<<< HEAD
=======
	  std::cout << "added subjets" << std::endl;
>>>>>>> Working single jet output
=======
=======
	    addSubJets(outTree, algorithms.AlgoNames[i], signal, i);
	  addJetsBranches(outTree, algorithms.AlgoNames[i], signal);
>>>>>>> Adding different algorithms to run on by using algorithms.xml.  It compiles, but still needs to be tested.
=======
	    addSubJets(outTree, algorithms.AlgoNames[i], i);
<<<<<<< HEAD
<<<<<<< HEAD
	  addJetsBranches(outTree, algorithms.AlgoNames[i]);
>>>>>>> Working version running over 14, lily, note and nc29.
=======
=======
=======
	    addSubJets(outTree, algorithmName, algorithm);
>>>>>>> Changed algorithm config maps to unordered maps, removed calls from loop.
	  // add branches for runnumber, mc_event_number etc.
>>>>>>> Updated weights file.  TLV works. Lots of protection put in against reading branches that are not valid.  Still need to add the option to turn off event and lepton selection tlection to do jet grooming algorithm analysis.
	  addInfoBranches(outTree);
>>>>>>> Fully commented.

<<<<<<< HEAD
<<<<<<< HEAD
>>>>>>> Some general cleanup. Set vars to -999 as default so that we can see when something doesn't get set. Particularly important for Topo jets which are not often found.
=======
	  // get number of entries to run over in the tree
>>>>>>> Updated weights file.  TLV works. Lots of protection put in against reading branches that are not valid.  Still need to add the option to turn off event and lepton selection tlection to do jet grooming algorithm analysis.
	  long entries = (long)inputTChain[tchainIdx]->GetEntries();

=======
>>>>>>> Changed some maps to unordered
	  // Setting up this with a high limit of 3.5 TeV so we don't miss anything.  Lots of bins - 200, so we can
	  // do a lot of tuning of the scale factor regions later on!
	  pt_reweight = new TH1F(std::string("pt_reweight"+bkg).c_str(),std::string("pt_reweight_"+bkg).c_str(), 200, 0, 3500);
	  cluster_vs_truthpt = new TH2F("ClusterVsTruthJetPt","ClusterVsTruthJetPt",40, 100, 1000, 40 ,0 , 1000);
	  cluster_vs_truthpt->GetXaxis()->SetTitle("Truth pT (GeV)");
	  cluster_vs_truthpt->GetYaxis()->SetTitle("Leading cluster pT (GeV)");
	  
	  // set up the number of events in total
	  NEvents = entries;
	  // we keep track of the weighted (by mc_event_weight) number of events too
	  NEvents_weighted.clear();

	  std::cout << "total entries: " << entries << std::endl;
	  
	  // counter for how many events pass selection
	  passed_counter = 0;

	  // mass variable
	  double mass = 0;
	  for (long n = 0; n < entries; n++)
	    {
	      // get next event
	      inputTChain[tchainIdx]->GetEntry(n);

	      // reset all of the output variables for next event
	      resetOutputVariables(); 
	     

	      if (n%1000==0)
		std::cout << "Entry: "<< n << " / " << entries <<  std::endl;

	      // increment event counter for MC channel
	      // Note that mc_event_weight is stored differently in xaod and d3pd
	      if (NEvents_weighted.find(mc_channel_number) != NEvents_weighted.end())
		{
		  if (xAOD)
		    NEvents_weighted[mc_channel_number] += mc_event_weight_xaod->at(0); // remove ->at(0) when running on D3PD
		  else
		    NEvents_weighted[mc_channel_number] += mc_event_weight_d3pd;// remove ->at(0) when running on D3PD
		}
	      else
		{
		  if (xAOD)
		    NEvents_weighted[mc_channel_number] = mc_event_weight_xaod->at(0); // remove ->at(0) when running on D3PD
		  else
		    NEvents_weighted[mc_channel_number] = mc_event_weight_d3pd; // remove ->at(0) when running on D3PD
		}

	      // initial values for the leading jet indices
	      int chosenLeadTruthJetIndex=-99;
	      int chosenLeadTopoJetIndex=-99;

	      // check that the leading truth jet is > 100 GeV
	       if ((*var_pt_vec[jetType::TRUTH])[0] / 1000.0 < 100)
		{
       		  continue;
		}
	      
	       // if we're doing a truth matching algorithm then don't apply pt range cuts on truth - it will be done on groomed
	      if (algorithmType.find("truthmatch") == std::string::npos) // check the pt is in the correct bin
		{
		  if ((*var_pt_vec[jetType::TRUTH])[0]/1000.0 < ptrange[j].first || (*var_pt_vec[jetType::TRUTH])[0]/1000.0 > ptrange[j].second)
		    continue;
		}

	      // apply extra lepton selections if running the hvtllqq selection
	      if (hvtllqq)
		{
		  // do overlap removal before looking for jets
		  overlapRemoval();
		  //apply event selection
		  int lepType = eventSelection();
		  if (lepType == leptonType::FAIL)
		    {
		      continue;
		    }
		  
		  // check if the lepton selection is passed for electron/ muon
		  if (!leptonSelection(lepType))
		    {
		      continue;
		    }
		}
	      
	      // do truth matching with topo jets
	      if (algorithmType.find("truthmatch") != std::string::npos)
		{
		  //check which is the CA12 ungroomed reco jet with EMfrac<0.99 and |eta|<1.2
		  //loop over all topo jets and get the leading with cuts
		  bool hasTopoJet=false;
		  
		  for (int jet_i=0; jet_i<(*var_pt_vec[jetType::TOPO]).size(); jet_i++)
		    {
		      // some samples don't have emfrac info, so this is a safeguard
		      float emfractmp = 0.5;
		      if (!xAOD)//emfrac)
			emfractmp = (*var_emfrac_vec[jetType::TOPO])[jet_i];
		      // check quality
		      if (!hasTopoJet && emfractmp<0.99 && fabs((*var_eta_vec[jetType::TOPO])[jet_i])<1.2) 
			{
			  chosenLeadTopoJetIndex=jet_i;
			  hasTopoJet=true;
			  // quit loop
			  continue;
			} 
		    } // end loop over var_pt_vec[jetType::TOPO]
		  
		  //now match truth jet with the chosen topo jet
		  if (hasTopoJet)
		    {
		      for (int jet_i=0; jet_i<(*var_pt_vec[jetType::TRUTH]).size(); jet_i++)
			{
			  // truth match on deltaR
			  if (chosenLeadTruthJetIndex<0 && DeltaR((*var_eta_vec[jetType::TOPO])[chosenLeadTopoJetIndex],(*var_phi_vec[jetType::TOPO])[chosenLeadTopoJetIndex],(*var_eta_vec[jetType::TRUTH])[jet_i],(*var_phi_vec[jetType::TRUTH])[jet_i])<0.9)
			    { 
			      chosenLeadTruthJetIndex=jet_i;
			      // quit loop
			      continue;
			    }	  
			}	// end loop over var_pt_vec[jetType::TRUTH]
		    } // end if(hasTopoJet)
		} // end if(groomAlgoIndex==0)
	      
	      // if truth match on topo jets was done, leave as is, otherwise use 0
	      chosenLeadTruthJetIndex = algorithmType.find("truthmatch") != std::string::npos ? chosenLeadTruthJetIndex : 0;
	      if (chosenLeadTruthJetIndex == 0)
		{
		  float maxpt = 0;
		  int count = 0;
		  for (vector<float>::iterator it = (*var_pt_vec[jetType::TRUTH]).begin(); it != (*var_pt_vec[jetType::TRUTH]).end() ; it++)
		    {
		      if ((*it)>maxpt)
			{
			  maxpt = (*it);
			  chosenLeadTruthJetIndex = count;
			}
		      count+=1;
		    }
		}

	      if (chosenLeadTruthJetIndex > 0)
		std::cout << "Found new leading truth jet!!!!" << std::endl;
	      
	      // veto the event if we have no good truth jets.
	      if (chosenLeadTruthJetIndex < 0)
		continue;

	      // truth boson matching done on the leading truth jet (signal only) to check if it is a Z or W boson.
	      int truthBosonIndex = -1;
	      if (truthBosonMatching && signal)
		{
		  int count = 0;
		  for (int jet_i = 0; jet_i < (*var_truthboson_eta_vec).size(); jet_i++)
		    {
		      // delta R matching of truth jet with truth boson
		      if (truthBosonIndex<0 && DeltaR((*var_eta_vec[jetType::TRUTH])[chosenLeadTruthJetIndex],(*var_phi_vec[jetType::TRUTH])[chosenLeadTruthJetIndex],(*var_truthboson_eta_vec)[jet_i],(*var_truthboson_phi_vec)[jet_i])<0.75*radius)
			{
			  truthBosonIndex = count;
			}
		      count++;
		    }
		  // check that a parent was found, if not, veto event
		  if (truthBosonIndex < 0)
		    {
		      std::cout << "No truth boson parent found, vetoing event." << std::endl;
		      continue;
		    }
		  // check that the parent is a W
		  if (abs((*var_truthboson_ID_vec)[truthBosonIndex]) != 24)
		    {
		      std::cout << "Truth boson parent is not a W, vetoing event." << std::endl;
		      continue;
		    }
		}

	      

	      // find groomed jets
	      int chosenLeadGroomedIndex=-99;
	      for (int jet_i=0; jet_i<(*var_pt_vec[jetType::GROOMED]).size(); jet_i++)
		{
		  // emfrac is not always available
		  float emfractmp = 0.5;
		  if (!xAOD)//emfrac)
		    emfractmp = (*var_emfrac_vec[jetType::GROOMED])[jet_i];
		  // truth match on dR with some eta and emfrac cuts
		  // is the dR < 0.75 ot 0.9???
		  if (chosenLeadGroomedIndex<0 && DeltaR((*var_eta_vec[jetType::TRUTH])[chosenLeadTruthJetIndex],(*var_phi_vec[jetType::TRUTH])[chosenLeadTruthJetIndex],(*var_eta_vec[jetType::GROOMED])[jet_i],(*var_phi_vec[jetType::GROOMED])[jet_i])<(0.75*radius) && emfractmp<0.99 && fabs((*var_eta_vec[jetType::GROOMED])[jet_i])<1.2)
		    {
		      chosenLeadGroomedIndex=jet_i;
		      continue;
		    }     
		} // end loop over var_pt_vec[jetType::GROOMED]
	      
	      if (chosenLeadGroomedIndex == -99) // failed selection
		{
		  continue;	      
		}
	      
	      // set output values
	      //leadGroomedIndex = chosenLeadGroomedIndex;
	      //leadTruthIndex = chosenLeadTruthJetIndex;
	      //leadTopoIndex = chosenLeadTopoJetIndex;

	      // set the mass to the leading groomed jet
	      mass = (*var_m_vec[2])[chosenLeadGroomedIndex]/1000.0 ;
	      if (mass < 0.0) // not sure why this happens, but it seems to happen
		continue;
	      // if we are not applying a mass window we do not apply any mass cuts
	      if (applyMassWindow && (mass > mass_max && mass < mass_min))
		{
		  continue;
		}
	      // set the mass of the two leptons and leading groomed jet
	      if (hvtllqq)
		setLLJMass(chosenLeadGroomedIndex);

	      // leading subjet index
	      int lead_subjet = 0;

	      // if we are calculating the subjet variables: massdrop and momentum balance
	      if (subjetscalc)
		{

		  std::vector<int> subjet_idx = (*subjet_index).at(chosenLeadGroomedIndex); // only groomed ones.....
		  // get the two leading subjets
		  std::pair<int,int> subjet_leading = getTwoLeadingSubjets(subjet_idx,var_subjets_pt_vec);

		  // calculate the mass drop
		  if (subjet_idx.size() <= 1)
		    {
		      var_massdrop = -999;
		      var_yt = -999;
		    }
		  else
		    {
		      lead_subjet = subjet_leading.first;
		    
		      // get kinematics for the subjets
		      double pt2 =  (*var_subjets_pt_vec)[subjet_leading.second];
		      double eta_1 =  (*var_subjets_eta_vec)[subjet_leading.first];
		      double eta_2 =  (*var_subjets_eta_vec)[subjet_leading.second];
		      double phi_1 = (*var_subjets_phi_vec)[subjet_leading.first];
		      double phi_2 = (*var_subjets_phi_vec)[subjet_leading.second];
		      // subjet mass
		      double subjet_mass =  (*var_subjets_m_vec)[subjet_leading.first];
		    
		      // calculate the mass drop - leading subjet/mass
		      double mu12 = subjet_mass/(mass*1000);
		      var_massdrop = mu12;
		    
		      // momentum balance
		      double dRsub12 = DeltaR (eta_1, phi_1, eta_2, phi_2);
		      Float_t yt = (pt2*dRsub12)/(mass*1000);		    
		      yt*=yt; // yt^2
		      var_yt = yt;
		    }
		} // end if(subjets)
	      // if we have the mass fraction and momentum balance in the input file already
	      else if (subjetspre)
		{
		  var_massdrop = var_massFraction_vec;//(*var_Mu12_vec[jetType::GROOMED])[chosenLeadGroomedIndex];//var_massFraction_vec;
		  var_yt = var_ktycut2_vec;
		}

	      // tau21 for truth, topo and groomed
	      var_Tau21[0]=(*var_Tau2_vec[0])[chosenLeadTruthJetIndex]/(*var_Tau1_vec[0])[chosenLeadTruthJetIndex];
	      var_Tau21[1]=(*var_Tau2_vec[1])[chosenLeadTopoJetIndex]/(*var_Tau1_vec[1])[chosenLeadTopoJetIndex];
	      var_Tau21[2]=(*var_Tau2_vec[2])[chosenLeadGroomedIndex]/(*var_Tau1_vec[2])[chosenLeadGroomedIndex];
	      
	    
	      // set up tauwta variables and zcut12
	      if (calcTauWTA21)//useBranch(string("TauWTA2TauWTA1"),true) && useBranch(string("TauWTA2"), true) && useBranch(string("TauWTA1"), true) )
		{
		  // tauwta21 for truth, topo and groomed
		  // veto events which have buggy values
		  if ((*var_TauWTA1_vec[0])[chosenLeadTruthJetIndex] > 0 && (*var_TauWTA2_vec[0])[chosenLeadTruthJetIndex] > 0)
		    var_TauWTA2TauWTA1[0]=(*var_TauWTA2_vec[0])[chosenLeadTruthJetIndex]/(*var_TauWTA1_vec[0])[chosenLeadTruthJetIndex];
		  else
		    var_TauWTA2TauWTA1[0]=-999;
		  if ((*var_TauWTA1_vec[2])[chosenLeadGroomedIndex] > 0 && (*var_TauWTA2_vec[2])[chosenLeadGroomedIndex] > 0)
		    var_TauWTA2TauWTA1[2]=(*var_TauWTA2_vec[2])[chosenLeadGroomedIndex]/(*var_TauWTA1_vec[2])[chosenLeadGroomedIndex];
		  else
		    var_TauWTA2TauWTA1[2]=-999;
		  if (chosenLeadTopoJetIndex == -99)
		    var_TauWTA2TauWTA1[1]=-99;
		  else if ((*var_TauWTA1_vec[1])[chosenLeadTopoJetIndex] > 0 && (*var_TauWTA2_vec[1])[chosenLeadTopoJetIndex] > 0)
		    var_TauWTA2TauWTA1[1]=(*var_TauWTA2_vec[1])[chosenLeadTopoJetIndex]/(*var_TauWTA1_vec[1])[chosenLeadTopoJetIndex];
		  else
		    var_TauWTA2TauWTA1[1]=-999;
		  response_TauWTA2TauWTA1 = var_TauWTA2TauWTA1[2]/var_TauWTA2TauWTA1[2];
		}
		

	      // calculate cluster based variables
	      // set up a vector<TLV> to store the cluster information
	      std::vector<TLorentzVector> groomedclusters;
	      // the following need clusters if they are being calculated
	      //std::cout << "about to create some clusters!" << std::endl;		  
	      if (calcQJets || calcFoxWolfram20 || calcSoftDrop || calcClusters)
		{
		  bool runGroomedClusters = false;
		  // create clusters, the method returns if the cluster creation was successful
		  runGroomedClusters = createClusters(jetType::GROOMED, chosenLeadGroomedIndex, groomedclusters);
		  //std::cout << "created clusters in cluster if statement" << std::endl;
		  if (calcQJets)
		    {
		      // groomed jet if clusters were found
		      if (runGroomedClusters)
			var_QjetVol[jetType::GROOMED] = calculateQJetsVol_v2(groomedclusters);
		    }
		  if (calcFoxWolfram20)
		    {
		      // groomed jet if clusters were found
		      if (runGroomedClusters)
			var_FoxWolfram20[jetType::GROOMED] = calculateFoxWolfram20(groomedclusters);
		    }
		  if (calcSoftDrop)
		    {
		      // groomed jet if clusters were found
		      if (runGroomedClusters)
			var_softdrop[jetType::GROOMED] = calculateSoftDropTag(groomedclusters);
		    }
		} // end if (calcQjets || calcFW || calcSD || calcClusters)

	      // if the EEC values need to be calculated by hand, unfortunately, we can
	      // no longer do this for truth jets
	      if (calcEEC)
	      {
		  // first calculate the ECF variables
		  if (groomedclusters.size() < 1)
		    createClusters(jetType::GROOMED, chosenLeadGroomedIndex, groomedclusters);
		  //std::cout << "created clusters in calcEEC if statement" << std::endl;		  
		  calculateECF(groomedclusters, jetType::GROOMED, 1);
		  calculateECF(groomedclusters, jetType::GROOMED, 2);

		  //EEC:C1 - exp 2, beta 1 for truth and groomed
		  /*var_EEC_C2_1[jetType::TRUTH] = calculateEEC(jetType::TRUTH, 1, 2);
		  var_EEC_C2_1[jetType::GROOMED] = calculateEEC(jetType::GROOMED, 1, 2);
		  //EEC:C2 - exp 2, beta 2 for truth and groomed
		  var_EEC_C2_2[jetType::TRUTH] = calculateEEC(jetType::TRUTH, 2, 2);
		  var_EEC_C2_2[jetType::GROOMED] = calculateEEC(jetType::GROOMED, 2, 2);
		  //EEC:D1 - exp 3, beta 1 for truth and groomed
		  var_EEC_D2_1[jetType::TRUTH] = calculateEEC(jetType::TRUTH, 1, 3);
		  var_EEC_D2_1[jetType::GROOMED] = calculateEEC(jetType::GROOMED, 1, 3);
		  //EEC:D2 - exp3, beta 2 for truth and groomed
		  var_EEC_D2_2[jetType::TRUTH] = calculateEEC(jetType::TRUTH, 2, 3);
		  var_EEC_D2_2[jetType::GROOMED] = calculateEEC(jetType::GROOMED, 2, 3);*/
		  } // calcEEC
	      // if we have the ECF variables the calculations for EEC are much simpler
	      else if (preCalcEEC)
		{
		  // set for truth and groomed
		  setEEC(jetType::TRUTH, chosenLeadTruthJetIndex);
		}
	      
	      setEEC(jetType::GROOMED, chosenLeadGroomedIndex);
	      // if we have the truth eec values we can calculate the response
	      if (preCalcEEC)
		{
		  response_EEC_C2_1 = var_EEC_C2_1[jetType::GROOMED]/var_EEC_C2_1[jetType::TRUTH];
		  response_EEC_C2_2 = var_EEC_C2_2[jetType::GROOMED]/var_EEC_C2_2[jetType::TRUTH];
		  response_EEC_D2_1 = var_EEC_D2_1[jetType::GROOMED]/var_EEC_D2_1[jetType::TRUTH];
		  response_EEC_D2_2 = var_EEC_D2_2[jetType::GROOMED]/var_EEC_D2_2[jetType::TRUTH];
		}

	      // make sure all of the other output variables have their values set
	      setOutputVariables(chosenLeadTruthJetIndex, chosenLeadTopoJetIndex, chosenLeadGroomedIndex, lead_subjet, algorithmName , prefix);
	      
	      // count how many entries have passed selection
	      passed_counter += 1;

	      // fill variables in output tree
	      outTree->Fill();	      
	      // fill pt reweighting histogram with leading truth jet
	      pt_reweight->Fill((*var_pt_vec[jetType::TRUTH])[chosenLeadTruthJetIndex]/1000.0);

	    } // end loop over nentries

	  // write the rweight th1f things to the outfile...
	  // calculate the overall reweight with a new histogram!
	  outTree->GetCurrentFile()->Write();
	  // write histograms to file
	  cluster_vs_truthpt->Write();
	  pt_reweight->Write();
	  pt_reweight_arr[tchainIdx] = (TH1F*)pt_reweight->Clone();
	  // stupid clone method needs this so that it doesn't delete this histo when closing the file
	  // http://root.cern.ch/phpBB3/viewtopic.php?f=3&t=11486
	  pt_reweight_arr[tchainIdx]->SetDirectory(0);

	  outTree->GetCurrentFile()->Close();

	  std::stringstream ss2; // store the name of the output file and include the i and j indices!
	  std::string bkg2 = signal ? "sig": "bkg";
	  // add all of the file name elements
	  ss2 << algorithm << fileid_global << "/" << ss.str() << "_" << bkg2 << ".nevents";
	  // open the output file
	  ofstream ev_out(ss2.str());
	  // write out the weighted events to a text file
	  for (std::map<long,float>::iterator it = NEvents_weighted.begin(); it!= NEvents_weighted.end(); it++)
	    ev_out << it->first << "," << NEvents_weighted[it->first] << std::endl;
	  ev_out.close();
	  clusterfile.close();
	  delete outfile;

	} // end loop of datatype

      // create the reweighting file
      //std::string fname = algorithm+fileid_global+"/" + ss.str() + ".ptweights";
      //createPtReweightFile(pt_reweight_arr[sampleType::BACKGROUND], pt_reweight_arr[sampleType::SIGNAL], fname);
    } // end loop over pt bins

} // makeMassWindowFile()

/*
 * Return the index of the two leading subjets (based on Pt).
 *
 * @param jet_idx The vector containing the indices of the subjets for a given jet.
 * @param subjet_pt A pointer to the subjet_pt branch.
 *
 * @return pair of ints for the indices of the first and second leading subjets.
 */
std::pair<int,int> getTwoLeadingSubjets(std::vector<int> & jet_idx, std::vector<float> *& subjet_pt)
{
  // jet_idx contains the indices in the subjet vector for the subjets of jet i
  double max_pt = 0;
  double sec_pt = 0;
  int max = -1;
  int sec = -1;
  // loop through all of the indices
  for (int i = 0; i < jet_idx.size(); i++)
    {
      int idx = jet_idx[i];
      // new maximum
      if (max_pt < (*subjet_pt)[idx])
	{
	  // set old maximum to second
	  sec_pt = max_pt;
	  sec = max;
	  max_pt = (*subjet_pt)[idx];
	  max = idx;
	}
      // second leading jet
      else if (sec_pt < (*subjet_pt)[idx])
	{
	  sec_pt = (*subjet_pt)[idx];
	  sec = idx;
	}
    }
  return std::make_pair(max,sec);
}

/*
 * Read in a list of the branches we would like to use.  This comes from an input file specified in the config file. Some
 * essential branches are added in case they are not in the input file.
 *
 * @param algorithm The name of the algorithm being run.  If this is a valid branches file it will use this variable
 *        as a filename, otherwise it will use algorithm+"_branches.txt".
 * @param brancharray A map<string,int> containing all of the branches in the input tree
 * @return vector of pairs containing all of the branches and if they should be on or off.
 */
vector<std::pair<std::string,bool> > getListOfJetBranches(std::string &algorithm, std::unordered_map<std::string, bool> & brancharray)
{
  // vector storing which branches to use
  vector<pair<string,bool> > branches;

  // some variables require cluster information, list them here:
  vector<string> clustervariables {"QJetsVol"};
  // keep track of if cluster variables are being used
  bool addClusterVariables = false;
  if (calcQJets || calcFoxWolfram20 || calcSoftDrop || calcClusters)
    addClusterVariables = true;
  
  // reset branchmap
  branchmap.clear();
  // input filename
  std::string filename = algorithm.find("branches.txt") == std::string::npos ? algorithm+"_branches.txt" : algorithm;
  std::cout << "using branches file: " << filename << std::endl;
  // open input file
  ifstream in(filename);
  string line;
  // keep a record of all the branches in the config file in a map
  while (getline(in, line))
    {
      // remove whitespace
      trim(line);
      // add to branchmap
      branchmap[line] = true;
      // check to see if this variable requires cluster info
      if (!addClusterVariables)
	    {
	      int brlen = line.length();
	      for (std::vector<std::string>::iterator it = clustervariables.begin(); it != clustervariables.end(); it++)
		{
		  if ((*it).length() <= brlen && line.compare(brlen-(*it).length(), (*it).length(), (*it)) == 0)
		    addClusterVariables =true;
		}
	    }
      
    }
  //add the clusters to the branchmap if we need them
  if (addClusterVariables)
    {
      std::cout << "Adding cluster variables" << std::endl;
      branchmap["cl_lc_n"] = true;
      branchmap["cl_lc_pt"] = true;
      branchmap["cl_lc_eta"] = true;
      branchmap["cl_lc_phi"] = true;
    }
  if (truthBosonMatching)
    {
      std::cout << "Adding truth boson variables" << std::endl;
      branchmap["truthBoson_eta"] = true;
      branchmap["truthBoson_phi"] = true;
      branchmap["truthBoson_ID"] = true;
    }

  // need to add some essential ones in case they get forgotten in that config file :)
  branchmap["RunNumber"] = true;
  branchmap["mc_channel_number"] = true;
  branchmap["vxp_n"] = true;
  branchmap["nVertices"] = true;
  branchmap["averageIntPerXing"] = true;
  branchmap["mc_event_weight"] = true;

  // need to add leptons
  branchmap["electrons"] = true;
  branchmap["el_pt"] = true;
  branchmap["el_eta"] = true;
  branchmap["el_phi"] = true;
  branchmap["el_ptcone20"] = true;
  branchmap["el_etcone20"] = true;
  branchmap["muons"] = true;
  branchmap["mu_pt"] = true;
  branchmap["mu_eta"] = true;
  branchmap["mu_phi"] = true;
  branchmap["mu_ptcone20"] = true;
  branchmap["mu_etcone20"] = true;
  branchmap["mu_charge"] = true;
  in.close();


  // loop through all of the branches in the tree and check if they exist in the branch map
  for (std::unordered_map<std::string,bool>::iterator tb = brancharray.begin(); tb != brancharray.end(); tb++)
    {
      std::string name = (*tb).first;
      trim(name); // remove whitespace
      if (branchmap.find(name) != branchmap.end() ) // add this to the branches if found in branchmap
	{
	  branches.push_back(make_pair(name, true));
	  
	}
      else
	{
	  branches.push_back(make_pair(name, false));
	}


    }
  
  return branches;
} // getListOfBranches()



/*
 * Return the full name of the jet branch, excluding the variable name (pt, m, etc) at the end.  It adds the prefix and
 * the abbreviated algorithm name together.
 *
 * @param samplePrefix The sample prefix - AntiKt10 or CamKt12 for example.
 * @param groomalgo The shortened version of the grooming algorithm, excluding the prefix.
 * @param addLC If "LC" needs to be added to the name.
 * @param i The index of the type of jet - jetType::TRUTH/GROOMED/TOPO.
 */
std::string returnJetType(std::string & samplePrefix, std::string & groomalgo, bool addLC, int i)
{
  std::string jet = "";
  std::string xaod = "";
  // right now we have an issue with xAOD where it adds "Jets" before _variable, so we need a quick fix for this
  if (xAODJets)
    {
      xaod = "Jets";
    }
  switch (i)
    {
    case jetType::TRUTH: // truth
      jet="jet_CamKt12Truth"+xaod+"_";
      break;
    case jetType::TOPO: // topo
      if (addLC)
	jet = "jet_" + samplePrefix + "LCTopo"+xaod+"_";
      else
	jet = "jet_" + samplePrefix + "Topo"+xaod+"_";
      break;
    default: // groomed
      if (addLC)
	jet = "jet_" + samplePrefix +"LC" + groomalgo+"_";
      else
	jet = "jet_" + samplePrefix + groomalgo + "_"; 	  
    }
  return jet;
} // returnJetType

/*
 * Return the full name of the subjet branch, excluding the variable name (pt, m, etc) at the end.  It adds the prefix and
 * the abbreviated algorithm name together.
 *
 * @param samplePrefix The sample prefix - AntiKt10 or CamKt12 for example.
 * @param groomalgo The shortened version of the grooming algorithm, excluding the prefix.
 * @param addLC If "LC" needs to be added to the name.
 */
std::string returnSubJetType(std::string & samplePrefix, std::string & groomalgo, bool addLC)
{
  if (addLC) //
    return std::string("jet_"+samplePrefix+"LC"+algorithms.subjetMap[groomalgo]+"_");
  return std::string("jet_"+samplePrefix+algorithms.subjetMap[groomalgo]+"_");

}// return subjettype

/*
 * Check if a branch is in the branchmap so we don't try to set branch addresses we aren't using.  On the first iteration
 * this is checked for all calls on the first event and adds them to a map.  This means for future events we no longer
 * have to search the branchmap for the key because we have ensured it is there - faster!
 *
 * @param branch Name of the branch to check
 * @param partialmatch Only look for substrings of keys that match the string
 *
 * @return bool indicating if it is in the branchmap
 */
bool useBranch(std::string const& br, bool partialmatch)
{
  string branch = br;
  // look for the key in the branch map
  if (!partialmatch)
    {
      if (branchmap.find(branch) != branchmap.end() && branchmap[branch])
	{
	  return true;
	}

    }
  // look for a partial match
  else
    {
      int brlen = branch.length();
      for (std::unordered_map<std::string, bool>::iterator it = branchmap.begin(); it != branchmap.end(); it++)
	{
	  if (it->first.length() >= brlen && it->first.compare(it->first.length()-brlen, brlen, branch) == 0)
	    {
	      // add this partial match to the branchmap to make the next search faster
	      branchmap[branch] = true;
	      return true;
	    }
	}

    }
  // didn't find it if we reach this point, return false
  return false;
}

/*
 * Adds essential output branches to the output file.
 *
 * @param tree Pointer to the output TTree.
 */
void addInfoBranches(TTree * tree)
{
  //tree->Branch("leadTruthIndex",&leadTruthIndex, "leadTruthIndex/I");
  //tree->Branch("leadTopoIndex",&leadTopoIndex, "leadTopoIndex/I");
  //tree->Branch("leadGroomedIndex",&leadGroomedIndex, "leadGroomedIndex/I");

  tree->Branch("normalisation",&normalisation, "normalisation/F");
  tree->Branch("NEvents",&NEvents,"NEvents/I");
  tree->Branch("RunNumber",&runNumberOut, "RunNumber/I");
  tree->Branch("k_factor", &var_k_factor, "k_factor/F");
  tree->Branch("filter_eff", &var_filter_eff, "filter_eff/F");
  tree->Branch("xs", &var_xs, "xs/F");

}// addInfoBranches

/*
 * Set up the branch address for a vector<TLorentzVector>.
 *
 * @param tree Reference to a TChain pointer.  This is the input TTree.
 * @param list Reference map<string,int> that contains a list of all branches in the tree.
 * @param vec Reference to a vector<TLV> pointer which will be used to read in from the file.
 * @param branch The name of the branch we are setting the address for.  This is not passed as reference because of how it is set up. Really it could be a const string &.
 * @return bool indicating if the branch was successfully set
 */
bool setVector(TChain *& tree, std::unordered_map<std::string, bool> & list, vector<TLorentzVector> *& vec, string branch)//const char * branch)
{
  // first check to see if the branch is actually in the tree, then check to see if we are interested in using it
  if (list.find(branch) != list.end() && useBranch(branch))
    tree->SetBranchAddress(branch.c_str(), &vec);
  else
    {
      std::cout << "missing branch " << branch << ", might cause unexpected behaviour, removing from branchmap" << std::endl;
      // erase it to avoid issues
      branchmap.erase(branch);
      return false;
    }
  return true;
}

/*
 * Set up the branch address for a vector<Int_t>.
 *
 * @param tree Reference to a TChain pointer.  This is the input TTree.
 * @param list Reference to a map<string,int> that contains a list of all branches in the tree.
 * @param vec Reference to a vector<Int_t> pointer which will be used to read in from the file.
 * @param branch The name of the branch we are setting the address for.  This is not passed as reference because of how it is set up. Really it could be a const string &.
 * @return bool indicating if the branch was successfully set
 */
bool setVector(TChain *& tree, std::unordered_map<std::string, bool> & list, vector<Int_t> *& vec, string branch)//const char * branch)
{
  // first check to see if the branch is actually in the tree, then check to see if we are interested in using it
  if (list.find(branch) != list.end() && useBranch(branch))
    tree->SetBranchAddress(branch.c_str(), &vec);
  else
    {
      std::cout << "missing branch " << branch << ", might cause unexpected behaviour, removing from branchmap" << std::endl;
      // erase it to avoid issues
      branchmap.erase(branch);
      return false;
    }
  return true;
}

/*
 * Set up the branch address for a vector<Float_t>.  Overloaded version of the one for vector<TLV>.
 *
 * @param tree Reference to a TChain pointer.  This is the input TTree.
 * @param list Reference to a map<string,int> that contains a list of all branches in the tree.
 * @param vec Reference to a vector<Float_t> pointer which will be used to read in from the file.
 * @param branch The name of the branch we are setting the address for.  This is not passed as reference because of how it is set up. Really it could be a const string &.
 * @return bool indicating if the branch was successfully set
 */
bool setVector(TChain *& tree, std::unordered_map<std::string,bool> & list, vector<Float_t> *& vec, string branch)//const char * branch)
{
  // first check to see if the branch is actually in the tree, then check to see if we are interested in using it
  if (list.find(branch)!=list.end() && useBranch(branch))
    tree->SetBranchAddress(branch.c_str(), &vec);
  else
    {
      std::cout << "missing branch " << branch << ", might cause unexpected behaviour, removing from branchmap" << std::endl;
      // erase it to avoid issues
      branchmap.erase(branch);
      return false;
    }
  return true;
}

/*
 * Erase a jet from the groomed jet collection.  It erases this entry from all jet variable collections.
 *
 * @param jet Integer giving position of the jet to be erased in the jet collection.
 */
void eraseJet(int jet)
{
  int i = jetType::GROOMED;
  // erasing an element/ range that doesn't exist it will cause undefined behaviour
  // so first we check if it is NULL and if it is not empty
  if (var_E_vec[i] != NULL && var_E_vec[i]->size() > jet)
    var_E_vec[i]->erase(var_E_vec[i]->begin()+jet);
  if (var_pt_vec[i] != NULL && var_pt_vec[i]->size() > jet)
    var_pt_vec[i]->erase(var_pt_vec[i]->begin()+jet);
  if (var_m_vec[i] != NULL && var_m_vec[i]->size() > jet)
    var_m_vec[i]->erase(var_m_vec[i]->begin()+jet);
  if (var_eta_vec[i] != NULL && var_eta_vec[i]->size() > jet)
    var_eta_vec[i]->erase(var_eta_vec[i]->begin()+jet);
  if (var_phi_vec[i] != NULL && var_phi_vec[i]->size() > jet)
    var_phi_vec[i]->erase(var_phi_vec[i]->begin()+jet);
  if (var_emfrac_vec[i] != NULL && var_emfrac_vec[i]->size() > jet)
    var_emfrac_vec[i]->erase(var_emfrac_vec[i]->begin()+jet);
  if (var_constit_index[i] != NULL && var_constit_index[i]->size() > jet)
    var_constit_index[i]->erase(var_constit_index[i]->begin()+jet);
  if (var_constit_n[i] != NULL && var_constit_n[i]->size() > jet)
    var_constit_n[i]->erase(var_constit_n[i]->begin()+jet);
  if (var_Tau1_vec[i] != NULL && var_Tau1_vec[i]->size() > jet)
    var_Tau1_vec[i]->erase(var_Tau1_vec[i]->begin()+jet);
  if (var_Tau2_vec[i] != NULL && var_Tau2_vec[i]->size() > jet)
    var_Tau2_vec[i]->erase(var_Tau2_vec[i]->begin()+jet);
  if (var_SPLIT12_vec[i] != NULL && var_SPLIT12_vec[i]->size() > jet)
    var_SPLIT12_vec[i]->erase(var_SPLIT12_vec[i]->begin()+jet);

  if (var_Dip12_vec[i] != NULL && var_Dip12_vec[i]->size() > jet)
    var_Dip12_vec[i]->erase(var_Dip12_vec[i]->begin()+jet);

  if (var_PlanarFlow_vec[i] != NULL && var_PlanarFlow_vec[i]->size() > jet)
    var_PlanarFlow_vec[i]->erase(var_PlanarFlow_vec[i]->begin()+jet);
  if (var_Angularity_vec[i] != NULL && var_Angularity_vec[i]->size() > jet)
    var_Angularity_vec[i]->erase(var_Angularity_vec[i]->begin()+jet);

  if (var_YFilt_vec != NULL && var_YFilt_vec->size() > jet)
    var_YFilt_vec->erase(var_YFilt_vec->begin()+jet);
  

  if (var_Aplanarity_vec[i] != NULL && var_Aplanarity_vec[i]->size() > jet)
    var_Aplanarity_vec[i]->erase(var_Aplanarity_vec[i]->begin()+jet);
  if (var_Sphericity_vec[i] != NULL && var_Sphericity_vec[i]->size() > jet)
    var_Sphericity_vec[i]->erase(var_Sphericity_vec[i]->begin()+jet);
  if (var_ThrustMaj_vec[i] != NULL && var_ThrustMaj_vec[i]->size() > jet)
    var_ThrustMaj_vec[i]->erase(var_ThrustMaj_vec[i]->begin()+jet);
  if (var_ThrustMin_vec[i] != NULL && var_ThrustMin_vec[i]->size() > jet)
    var_ThrustMin_vec[i]->erase(var_ThrustMin_vec[i]->begin()+jet);

  if (var_TauWTA1_vec[i] != NULL && var_TauWTA1_vec[i]->size() > jet)
    var_TauWTA1_vec[i]->erase(var_TauWTA1_vec[i]->begin()+jet);
  if (var_TauWTA2_vec[i] != NULL && var_TauWTA2_vec[i]->size() > jet)
    var_TauWTA2_vec[i]->erase(var_TauWTA2_vec[i]->begin()+jet);
  if (var_ZCUT12_vec[i] != NULL && var_ZCUT12_vec[i]->size() > jet)
    var_ZCUT12_vec[i]->erase(var_ZCUT12_vec[i]->begin()+jet);

  if (var_Mu12_vec[i] != NULL && var_Mu12_vec[i]->size() > jet)
    var_Mu12_vec[i]->erase(var_Mu12_vec[i]->begin()+jet);

  // only erase these variables if we are not calculating them
  if (!calcFoxWolfram20 && !preCalcFoxWolfram20)
    {
      if (var_FoxWolfram0_vec[i] != NULL && var_FoxWolfram0_vec[i]->size() > jet)
	var_FoxWolfram0_vec[i]->erase(var_FoxWolfram0_vec[i]->begin()+jet);
      if (var_FoxWolfram2_vec[i] != NULL && var_FoxWolfram2_vec[i]->size() > jet)
	var_FoxWolfram2_vec[i]->erase(var_FoxWolfram2_vec[i]->begin()+jet);
    }
  else if (preCalcFoxWolfram20)
    {
      if (var_FoxWolfram20_vec[i] != NULL && var_FoxWolfram20_vec[i]->size() > jet)
	var_FoxWolfram20_vec[i]->erase(var_FoxWolfram20_vec[i]->begin()+jet);
    }


  if (!calcSoftDrop && var_SoftDropTag_vec[i] != NULL && var_SoftDropTag_vec[i]->size() > jet)
    var_SoftDropTag_vec[i]->erase(var_SoftDropTag_vec[i]->begin()+jet);
  
} // eraseJet

/*
 * Overlap removal between jets and electrons, removes any (groomed) jet that overlaps with an electron.
 *
 */
void overlapRemoval()
{
  // loop through the groomed jets
  for (int it = 0 ; it < (*var_pt_vec[jetType::GROOMED]).size(); it++)
    {
      // loop through the electron jets
      for (std::vector<TLorentzVector>::iterator el = var_electrons_vec->begin(); el != var_electrons_vec->end(); el++)
	{
	  // calc deltaR between jet and electron
	  float dR = DeltaR((*var_eta_vec[jetType::GROOMED])[it], (*var_phi_vec[jetType::GROOMED])[it], (*el).Eta(), (*el).Phi());
	  // dR < 0.8 remove jet from collection
	  if (dR < 0.8)
	    {
	      // remove jet from collection
	      eraseJet(it);
	      // account for removed jet
	      it--;
	    }
	}
    }
  
} // overlapRemoval

/*
 * Set the invariant mass of the two leptons and jet.
 */
void setLLJMass(int jetidx)
{
  // if we aren't running the lepton selection then just return
  if (!hvtllqq)
    return;
  // set the invariant mass of the two leptons and groomed jet
  int t = jetType::GROOMED;
  TLorentzVector j = TLorentzVector();
  j.SetPtEtaPhiM((*var_pt_vec[t])[jetidx], (*var_eta_vec[t])[jetidx], (*var_phi_vec[t])[jetidx] ,(*var_m_vec[t])[jetidx]);
  // add the jet and two leptons and get the mass.
  var_mllj = (j+var_leptons[0]+var_leptons[1]).M();
} // setLLJMass

/*
 * Apply basic event selection - two opposite charge muons or two electrons and the invariant mass of the 
 * two leptons minus the Z mass must be < 25 GeV
 *
 * @return leptonType::FAIL, ELECTRON or MUON.
*/
int eventSelection()
{
  if (!hvtllqq) // not running lepton selection
    {
      std::cout << "running event selection when analysis is not being done, should not call this method" << std::endl;
      return -1;
    }
  int lepType = leptonType::FAIL;
  // 2 electrons or 2 muons of opposite charge  
  if (var_electrons_vec->size() == 2)
    {
      // set mass and pt of combined lepton system
      var_mll = (var_electrons_vec->at(0)+var_electrons_vec->at(1)).M();
      var_ptll = (var_electrons_vec->at(0)+var_electrons_vec->at(1)).Pt();
      lepType = leptonType::ELECTRON;
    }
  else if (var_muons_vec->size() == 2 && (*var_mu_charge_vec)[0]*(*var_mu_charge_vec)[1] == -1)
    {
      // set mass and pt of combined lepton system
      var_mll = (var_muons_vec->at(0)+var_muons_vec->at(1)).M();
      var_ptll = (var_muons_vec->at(0)+var_muons_vec->at(1)).Pt();
      lepType = leptonType::MUON;
    }
  else
    {
      // it has failed selection
      return lepType;
    }
  // |mll-91| < 25 gev and ptll > 70 gev to be consistent with ttbar
  if (fabs(var_mll-91*GEV) >= 25*GEV || fabs(var_ptll) <= 70*GEV)
    {
      return leptonType::FAIL;
    }

  return lepType;
  
}//eventSelection

/**
 * Method to return a vector with dummy charges for electron events where we have no charge info.
 *
 * @param size Number of dummy entries to insert.  Generally this is two.
 *
 * @return vector<float> of size electrons.size()
 */
std::vector<float> dummyCharge(int size)
{
  vector<float> vec;
  for (int i  = 0; i < size; i++)
    {
      vec.push_back(-999);
    }
  return vec;
} //dummyCharge

/*
 * Apply lepton selection. 
 * If it is an electron event (2 electrons) we want to have:
 * Et >25 GeV, ptcone20/pt < 0.15, etcone20/et < 0.3 and |eta|<2.47.
 * If it is a muon event (2 muons) we want to have:
 * Pt >25 GeV, ptcone20/pt < 0.15, etcone20/et < 0.3 and |eta|<2.5.
 *
 * @param lepType The type of event, which is either leptonType::ELECTRON or leptonType::MUON
 *
 * @return Whether or not lepton selection is passed.
 */
bool leptonSelection(int lepType)
{
  if (!hvtllqq) // if we are not running the hvtllqq analysis no lep selection
    return true;
  bool pass = true;
  // clear all of the vectors that we are going to set here
  var_etcone20.clear();
  var_ptcone20.clear();
  var_charge.clear();
  var_leptons.clear();

  if (lepType == leptonType::ELECTRON)
    {
      std::vector<float> dummyChargeVec = dummyCharge(2);
      // electron selection

      for (int idx = 0; idx < 2; idx++)
	{
	  // et > 25 GeV
	  if ((*var_electrons_vec)[idx].Et() <= 25*GEV)
	    {
	      return false;
	    }
	  // eta < 2.47
	  if (fabs((*var_electrons_vec)[idx].Eta()) >= 2.47)
	    {
	      return false;
	    }
	  // ptcone20/pt < 0.15
	  if ((*var_el_ptcone20_vec)[idx]/(*var_electrons_vec)[idx].Pt() >= 0.15)
	    {
	      return false;
	    }
	  // etcone20/et < 0.3
	  if ((*var_el_etcone20_vec)[idx]/(*var_electrons_vec)[idx].Et() >= 0.3)
	    {
	      return false;
	    }
	  var_ptcone20.push_back((*var_el_ptcone20_vec)[idx]);
	  var_etcone20.push_back((*var_el_etcone20_vec)[idx]);
	  var_charge.push_back(dummyChargeVec[idx]);
	  var_leptons.push_back((*var_electrons_vec)[idx]);
	} // for loop idx
      
      var_isElectronEvent = 1;
      
    } // if leptype = electron
  else if (lepType == leptonType::MUON)
    {
      // muons
      for (int idx = 0; idx < 2; idx++)
	{
	  // pt > 25 gev
	  if ((*var_muons_vec)[idx].Pt() <= 25*GEV)
	    {
	      return false;
	    }
	  // eta < 2.5
	  if (fabs((*var_muons_vec)[idx].Eta()) >= 2.5)
	    {
	      return false;
	    }
	  // ptcone20/pt <  0.15
	  if ((*var_mu_ptcone20_vec)[idx]/(*var_muons_vec)[idx].Pt() >= 0.15)
	    {
	      return false;
	    }
	  // etcone20/et < 0.3
	  if ((*var_mu_etcone20_vec)[idx]/(*var_muons_vec)[idx].Et() >= 0.3)
	    {
	      return false;
	    }
	  var_ptcone20.push_back((*var_mu_ptcone20_vec)[idx]);
	  var_etcone20.push_back((*var_mu_etcone20_vec)[idx]);
	  var_charge.push_back((*var_mu_charge_vec)[idx]);
	  var_leptons.push_back((*var_muons_vec)[idx]);
	} // for loop idx

      var_isElectronEvent = 0;
    } // elif lepttype = muon
  else 
    return false;

  return pass;
} //leptonSelection


/*
 * Set all of the branch addresses for the jets we are reading in.  If a branch is missing in the tree it is not set and a warning is given.
 * 
 * @param tree A pointer to the TChain that is being read in.
 * @param groomalgo The shortened name of the algorithm being used, like TopoSplitFiltered
 * @param groomIdx The full name of the algorithm.
 */
void setJetsBranches(TChain * tree, std::string &groomalgo,  std::string & groomIdx, std::unordered_map<std::string,bool> & brancharray)
{
  std::string samplePrefix = "";

  bool addLC = false; // we use this variable to indicate whether "LC" should be in the algorithm string
  samplePrefix = algorithms.AlgoPrefix[groomIdx];
  
  if (algorithms.AlgoType[groomIdx].find("recluster") == std::string::npos) // we're doing reclustering
    addLC = true; // just add the LC to the name

  // these branches should always exist
  if (xAOD)
    tree->SetBranchAddress("mc_event_weight",&mc_event_weight_xaod);
  else
    tree->SetBranchAddress("mc_event_weight",&mc_event_weight_d3pd);
  tree->SetBranchAddress("mc_channel_number", &mc_channel_number);

  // as should these, but double check to make sure the branch exists in the tree
  if (brancharray.find("RunNumber") != brancharray.end() && branchmap["RunNumber"])
    tree->SetBranchAddress("RunNumber", &runNumberIn);
  if (brancharray.find("nVertices")  != brancharray.end() && branchmap["nVertices"])
    tree->SetBranchAddress("nVertices", &nvtxIn);
  if (brancharray.find("vxp_n")  != brancharray.end() && branchmap["vxp_n"])
    tree->SetBranchAddress("vxp_n", &nvtxIn);
  if (brancharray.find("averageIntPerXing") != brancharray.end() && branchmap["averageIntPerXing"])
    {
      if (xAOD)
	{
	  tree->SetBranchAddress("averageIntPerXing",&avgIntpXingIn_xaod);
	  tree->SetBranchAddress("actualIntPerXing",&actualIntPerXingIn);
	}
      else
	tree->SetBranchAddress("averageIntPerXing",&avgIntpXingIn_d3pd);
    }
  if (brancharray.find("evt_scale1fb") != brancharray.end() && branchmap["evt_scale1fb"])
    tree->SetBranchAddress("evt_scale1fb",&scale1fb);
  if (xAOD)
    {
      if (brancharray.find("evt_kfactor") != brancharray.end() && branchmap["evt_kfactor"])
	tree->SetBranchAddress("evt_kfactor",&evt_kfactor);
      if (brancharray.find("evt_filtereff") != brancharray.end() && branchmap["evt_filtereff"])
	tree->SetBranchAddress("evt_filtereff",&evt_filtereff);
      if (brancharray.find("evt_nEvts") != brancharray.end() && branchmap["evt_nEvts"])
	tree->SetBranchAddress("evt_nEvts",&evt_nEvts);
      if (brancharray.find("evt_sumWeights") != brancharray.end() && branchmap["evt_sumWeights"])
	tree->SetBranchAddress("evt_sumWeights",&evt_sumWeights);
      if (brancharray.find("evt_xsec") != brancharray.end() && branchmap["evt_xsec"])
      	tree->SetBranchAddress("evt_xsec",&evt_xsec);
    }


  
  // set up the branch addresses for the lepton variables if running hvtllqq analysis
  // if the setVector method fails, then set it to a dummy vector with default values
  if (hvtllqq)
    {
      if (!setVector(tree, brancharray, var_electrons_vec, std::string("electrons")))
	{
	  //var_electrons_vec = tlvvec();
	  tlvvec(var_electrons_vec);
	}
      if (!setVector(tree, brancharray, var_electronPt_vec, std::string("el_pt")))
	{
	  floatvec(var_electronPt_vec);
	  //var_electronPt_vec = floatvec();
	}
      if (!setVector(tree, brancharray, var_electronEta_vec, "el_eta"))
	floatvec(var_electronEta_vec);
	//var_electronEta_vec = floatvec();
      if (!setVector(tree, brancharray, var_electronPhi_vec, "el_phi"))
	floatvec(var_electronPhi_vec);
      //var_electronPhi_vec = floatvec();
      if (!setVector(tree, brancharray, var_el_ptcone20_vec, "el_ptcone20"))
	floatvec(var_el_ptcone20_vec);
	//var_el_ptcone20_vec = floatvec();
      if (!setVector(tree, brancharray, var_el_etcone20_vec, "el_etcone20"))
	floatvec(var_el_etcone20_vec);
      //var_el_etcone20_vec = floatvec();
      
      if (!setVector(tree, brancharray, var_muons_vec, "muons"))
	tlvvec(var_muons_vec);
	//var_muons_vec = tlvvec();
      
      if (!setVector(tree, brancharray, var_muonPt_vec, "mu_pt"))
	floatvec(var_muonPt_vec);
      //var_muonPt_vec = floatvec();
      if (!setVector(tree, brancharray, var_muonEta_vec, "mu_eta"))
	floatvec(var_muonEta_vec);
      //var_muonEta_vec = floatvec();
      if (!setVector(tree, brancharray, var_muonPhi_vec, "mu_phi"))
	floatvec(var_muonPhi_vec);
      //var_muonPhi_vec = floatvec();
      if (!setVector(tree, brancharray, var_mu_ptcone20_vec, "mu_ptcone20"))
	floatvec(var_mu_ptcone20_vec);
      //var_mu_ptcone20_vec = floatvec();
      if (!setVector(tree, brancharray, var_mu_etcone20_vec, "mu_etcone20"))
	floatvec(var_mu_etcone20_vec);
      //var_mu_etcone20_vec = floatvec();
      if (!setVector(tree, brancharray, var_mu_charge_vec, "mu_charge"))
	floatvec(var_mu_charge_vec);
      //var_mu_charge_vec = floatvec();


    }


  // loop through the truth, toppo and groomed jets and set up the branches for the different variables
  for (int i = 0; i < jetType::MAX; i++) // truth, topo, groomed
    {
      std::string jetString = returnJetType(samplePrefix, groomalgo, addLC,i); //set to truth/ topo/ groomed
      if (!setVector(tree, brancharray, var_E_vec.at(i), std::string(jetString+"E") ))
	floatvec(var_E_vec[i]);

      if (!setVector(tree, brancharray, var_pt_vec.at(i), std::string(jetString+"pt") ))
	floatvec(var_pt_vec[i]);

      if (!setVector(tree, brancharray, var_m_vec.at(i), std::string(jetString+"m") ))
	floatvec(var_m_vec[i]);

      if (!setVector(tree, brancharray, var_eta_vec.at(i), std::string(jetString+"eta") ))
	floatvec(var_eta_vec[i]);

      if (!setVector(tree, brancharray, var_phi_vec.at(i), std::string(jetString+"phi") ))
	floatvec(var_phi_vec[i]);

      if (!setVector(tree, brancharray, var_emfrac_vec.at(i), std::string(jetString+"emfrac") ))
	floatvec(var_emfrac_vec[i]);

      if (!setVector(tree, brancharray, var_Tau1_vec.at(i), std::string(jetString+"Tau1") ))
	floatvec(var_Tau1_vec[i]);

      if (!setVector(tree, brancharray, var_Tau2_vec.at(i), std::string(jetString+"Tau2") ))
	floatvec(var_Tau2_vec[i]);

      if (!setVector(tree, brancharray, var_constit_n.at(i), std::string(jetString+"constit_n") ))
	intvec(var_constit_n[i]);

      if (brancharray.find(std::string(jetString+"constit_index")) != brancharray.end() && useBranch(std::string(jetString+"constit_index")))
	tree->SetBranchAddress(std::string(jetString+"constit_index").c_str(), &var_constit_index.at(i));
      else
	vecintvec(var_constit_index[i]);

      if (!setVector(tree, brancharray, var_SPLIT12_vec.at(i), std::string(jetString+"SPLIT12") ))
	floatvec(var_SPLIT12_vec[i]);

      if (!setVector(tree, brancharray, var_Dip12_vec.at(i), std::string(jetString+"Dip12") ))
	floatvec(var_Dip12_vec[i]);

      if (!setVector(tree, brancharray, var_PlanarFlow_vec.at(i), std::string(jetString+"PlanarFlow") ))
	floatvec(var_PlanarFlow_vec[i]);

      if (!setVector(tree, brancharray, var_Angularity_vec.at(i), std::string(jetString+"Angularity") ))
	floatvec(var_Angularity_vec[i]);

      if (!setVector(tree, brancharray, var_Aplanarity_vec.at(i), std::string(jetString+"Aplanarity") ))
	floatvec(var_Aplanarity_vec[i]);

      if (!setVector(tree, brancharray, var_Sphericity_vec.at(i), std::string(jetString+"Sphericity") ))
	floatvec(var_Sphericity_vec[i]);

      if (!setVector(tree, brancharray, var_ThrustMaj_vec.at(i), std::string(jetString+"ThrustMaj") ))
	floatvec(var_ThrustMaj_vec[i]);

      if (!setVector(tree, brancharray, var_ThrustMin_vec.at(i), std::string(jetString+"ThrustMin") ))
	floatvec(var_ThrustMin_vec[i]);

      if (!setVector(tree, brancharray, var_TauWTA1_vec.at(i), std::string(jetString+"TauWTA1") ))
	floatvec(var_TauWTA1_vec[i]);

      if (!setVector(tree, brancharray, var_TauWTA2_vec.at(i), std::string(jetString+"TauWTA2") ))
	floatvec(var_TauWTA2_vec[i]);

      if (!setVector(tree, brancharray, var_ZCUT12_vec.at(i), std::string(jetString+"ZCUT12") ))
	floatvec(var_ZCUT12_vec[i]);

      if (!setVector(tree, brancharray, var_ECF1_vec.at(i), std::string(jetString+"ECF1") ))
	floatvec(var_ECF1_vec[i]);
      if (!setVector(tree, brancharray, var_ECF2_vec.at(i), std::string(jetString+"ECF2") ))
	floatvec(var_ECF2_vec[i]);
      if (!setVector(tree, brancharray, var_ECF3_vec.at(i), std::string(jetString+"ECF3") ))
	floatvec(var_ECF3_vec[i]);

      if (!setVector(tree, brancharray, var_ECF1_2_vec.at(i), std::string(jetString+"ECF1_beta2") ))
	floatvec(var_ECF1_vec[i]);
      if (!setVector(tree, brancharray, var_ECF2_2_vec.at(i), std::string(jetString+"ECF2_beta2") ))
	floatvec(var_ECF2_vec[i]);
      if (!setVector(tree, brancharray, var_ECF3_2_vec.at(i), std::string(jetString+"ECF3_beta2") ))
	floatvec(var_ECF3_vec[i]);

      if (!calcFoxWolfram20 && !preCalcFoxWolfram20)
	{
	  if (!setVector(tree, brancharray, var_FoxWolfram0_vec.at(i), std::string(jetString+"FoxWolfram_0") ))
	    floatvec(var_FoxWolfram0_vec[i]);
	  if (!setVector(tree, brancharray, var_FoxWolfram2_vec.at(i), std::string(jetString+"FoxWolfram_2") ))
	    floatvec(var_FoxWolfram2_vec[i]);
	}
      else if (preCalcFoxWolfram20)
	{
	  if (!setVector(tree, brancharray, var_FoxWolfram20_vec.at(i), std::string(jetString+"FoxWolfram2") ))
	    floatvec(var_FoxWolfram20_vec[i]);
	}
      
      if (!calcSoftDrop && !setVector(tree, brancharray, var_SoftDropTag_vec.at(i), std::string(jetString+"SoftDropTag") ))
	intvec(var_SoftDropTag_vec[i]);

      if (!setVector(tree,brancharray, var_Mu12_vec.at(i), std::string(jetString+"Mu12") ))
	{
	  floatvec(var_Mu12_vec[i]);
	}



    } // end for loop over topo/truth/groom

  // setup the truth boson branches - note this is only done on the xAODs
  //if (truthBosonMatching)
  //{
      if (!setVector(tree,brancharray, var_truthboson_eta_vec, std::string("truthBoson_eta")))
	floatvec(var_truthboson_eta_vec);
      if (!setVector(tree,brancharray, var_truthboson_phi_vec, std::string("truthBoson_phi")))
	floatvec(var_truthboson_phi_vec);
      if (!setVector(tree,brancharray, var_truthboson_ID_vec, std::string("truthBoson_ID")))
	intvec(var_truthboson_ID_vec);
      //}

  
      // there is only yfilt stored for groomed jets
  std::string jetstring = returnJetType(samplePrefix, groomalgo, addLC,2); //set to groomed
  if (!setVector(tree, brancharray, var_YFilt_vec, std::string(jetstring+"YFilt") ))
    floatvec(var_YFilt_vec);




  // set up cluster variables
  // if a branch is missing set up dummy values
  if (brancharray.find(std::string("cl_lc_n")) != brancharray.end() && useBranch(string("cl_lc_n")))
    tree->SetBranchAddress(std::string("cl_lc_n").c_str(), &var_cl_n );
  else
    {
      std::cout << "branch not found: cl_lc_n" << std::endl;
      var_cl_n = -99;
    }
  if (!setVector(tree, brancharray, var_cl_pt_vec, std::string("cl_lc_pt")))
    floatvec(var_cl_pt_vec);
    //var_cl_pt_vec = floatvec();
  if (!setVector(tree, brancharray, var_cl_eta_vec, std::string("cl_lc_eta")))
    floatvec(var_cl_eta_vec);
  //var_cl_eta_vec = floatvec();
  if (!setVector(tree, brancharray, var_cl_phi_vec, std::string("cl_lc_phi")))
    floatvec(var_cl_phi_vec);
  //var_cl_phi_vec = floatvec();

  std::string jetString = returnJetType( samplePrefix, groomalgo, addLC, 2); //set to truth/ topo/ groomed

  if (subjetspre) // if the sample has the pre-calculated subjet variables
    {
      if(brancharray.find(std::string(jetString+"config_massFraction"))  != brancharray.end() )
	tree->SetBranchAddress(std::string(jetString+"config_massFraction").c_str(),&var_massFraction_vec);  
      else
	{
	  std::cout << "branch not found: " << std::string(jetString+"config_massFraction/Mu12") << std::endl;
	  var_massFraction_vec = -99;
	}
      if(brancharray.find(std::string(jetString+"config_ktycut2"))  != brancharray.end())
	tree->SetBranchAddress(std::string(jetString+"config_ktycut2").c_str(),&var_ktycut2_vec);
      else
	{
	  std::cout << "branch not found: " << std::string(jetString+"config_ktycut2") << std::endl;
	  var_ktycut2_vec = -99;
	}
    }

  if (!subjetscalc) // if the sample has the subjet branches available to create massdrop and yt
    return;

  std::string subjetType = returnSubJetType(samplePrefix, groomIdx, addLC);

  // set the branch status and address for the branch name that has the subjets for the grooming algorithm
  tree->SetBranchStatus(std::string(algorithms.subjetIndex[groomIdx]).c_str(),1);
  tree->SetBranchAddress(std::string(algorithms.subjetIndex[groomIdx]).c_str(),&subjet_index);

  if(brancharray.find(std::string(algorithms.subjetIndex[groomIdx])) == brancharray.end()) {
    std::cout << "subjet branch is not here, change the config file otherwise segfaults will come for ye" << std::endl;
    exit(EXIT_FAILURE);
    return;
  }

  // subjet kinematics
  tree->SetBranchAddress(std::string(subjetType+"E").c_str(),&var_subjets_E_vec);
  tree->SetBranchAddress(std::string(subjetType+"pt").c_str(),&var_subjets_pt_vec);
  tree->SetBranchAddress(std::string(subjetType+"m").c_str(),&var_subjets_m_vec);
  tree->SetBranchAddress(std::string(subjetType+"eta").c_str(),&var_subjets_eta_vec);
  tree->SetBranchAddress(std::string(subjetType+"phi").c_str(),&var_subjets_phi_vec);
  
}//setJetsBranches



/*
 * Add the sub-jet branches to the output file.
 * @param tree A TChain pointer for the output file.
 * @param groomalgo A shortened version of the algorithm name.
 * @param groomIdx The full algorithm name. It is used as a key in the algorithms maps.
 */
void addSubJets(TTree * tree, std::string & groomalgo, std::string &  groomIdx)
{

  std::string samplePrefix = ""; // AntiKt10/ CamKt12 for example
  bool addLC = false; // some algorithms have "LC" in their name
  samplePrefix = algorithms.AlgoPrefix[groomIdx];

  if (algorithms.AlgoType[groomIdx].find("recluster") == std::string::npos) // we're doing reclustering
    addLC = true; // just add the LC to the name

  std::string jetString = returnJetType( samplePrefix, groomalgo, addLC, 2); //set to truth/ topo/ groomed
  // subjet string
  std::string subjetType = returnSubJetType(samplePrefix, groomIdx, addLC);
  if (subjetscalc)
    {
      tree->Branch(std::string(subjetType+"E").c_str(),&var_subjets_E,std::string(subjetType+"E/F").c_str());
      tree->Branch(std::string(subjetType+"pt").c_str(),&var_subjets_pt,std::string(subjetType+"pt/F").c_str());
      tree->Branch(std::string(subjetType+"m").c_str(),&var_subjets_m,std::string(subjetType+"m/F").c_str());
      tree->Branch(std::string(subjetType+"eta").c_str(),&var_subjets_eta,std::string(subjetType+"eta/F").c_str());
      tree->Branch(std::string(subjetType+"phi").c_str(),&var_subjets_phi,std::string(subjetType+"phi/F").c_str());
    }

  tree->Branch(std::string(jetString+"massdrop").c_str(),&var_massdrop, std::string(jetString+"massdrop/F").c_str());
  tree->Branch(std::string(jetString+"yt").c_str(),&var_yt,std::string(jetString+"yt/F").c_str());

} //addSubJets()


/*
 * Add the branches for the leptons in the output file.
 *
 * @param tree A TChain pointer to the output file.
 */
void addLeptonBranches(string & jetString, TTree * tree)
{
  if (!hvtllqq) // if not running hvtllqq analysis
    return;
  tree->Branch("leptons",&var_leptons);
  tree->Branch("ptcone20", &var_ptcone20);
  tree->Branch("etcone20", &var_etcone20);
  tree->Branch("charge", &var_charge);
  tree->Branch(string(jetString+"mass_llj").c_str(), &var_mllj, "mass_mllj/F");
  tree->Branch("mll", &var_mll, "mll/F");
  tree->Branch("ptll", &var_ptll, "ptll/F");
  tree->Branch("isElectronEvent", &var_isElectronEvent, "isElectronEvent/I");

} // addLeptonBranches


/*
 * Initialise all of the vectors used for reading in the input TTree.
 *
 */
void initVectors()
{

  for (int i = 0; i < 3; i++) // for jetType::TRUTH/TOPO/GROOMED
    {

      var_E_vec.push_back(0);
      var_pt_vec.push_back(0);
      var_m_vec.push_back(0);
      var_eta_vec.push_back(0);
      var_phi_vec.push_back(0);
      var_emfrac_vec.push_back(0);
      var_constit_index.push_back(0);
      var_constit_n.push_back(0);
      var_Tau1_vec.push_back(0);
      var_Tau2_vec.push_back(0);
      var_SPLIT12_vec.push_back(0);
      var_Dip12_vec.push_back(0);
      var_PlanarFlow_vec.push_back(0);
      var_Angularity_vec.push_back(0);

      var_Aplanarity_vec.push_back(0);
      var_Sphericity_vec.push_back(0);
      var_ThrustMaj_vec.push_back(0);
      var_ThrustMin_vec.push_back(0);
  
      var_TauWTA1_vec.push_back(0);
      var_TauWTA2_vec.push_back(0);
      
      var_ZCUT12_vec.push_back(0);
      var_ECF1_vec.push_back(0);
      var_ECF2_vec.push_back(0);
      var_ECF3_vec.push_back(0);
      var_ECF1_2_vec.push_back(0);
      var_ECF2_2_vec.push_back(0);
      var_ECF3_2_vec.push_back(0);
      var_Mu12_vec.push_back(0);
      

      if (!calcFoxWolfram20 && !preCalcFoxWolfram20)
	{
	  var_FoxWolfram0_vec.push_back(0);
	  var_FoxWolfram2_vec.push_back(0);
	}
      else if (preCalcFoxWolfram20)
	var_FoxWolfram20_vec.push_back(0);
      
      if (!calcSoftDrop)
	var_SoftDropTag_vec.push_back(0);
    }  

  // truth bosons
  var_truthboson_eta_vec = 0;
  var_truthboson_phi_vec = 0;
  var_truthboson_ID_vec = 0;

  var_YFilt_vec = 0;
  var_massFraction_vec = 0;
  var_ktycut2_vec = 0;
  //subjets
  subjet_index = 0;    
  var_subjets_E_vec = 0;
  var_subjets_pt_vec = 0;
  var_subjets_m_vec = 0;
  var_subjets_eta_vec = 0;
  var_subjets_phi_vec = 0;
  // clusters
  var_cl_pt_vec = 0;
  var_cl_n = 0;
  var_cl_eta_vec = 0;
  var_cl_phi_vec = 0;
  // leptons
  var_electrons_vec = 0;

  var_electronPhi_vec = 0;
  var_electronEta_vec = 0;
  var_electronPt_vec = 0;
  var_el_ptcone20_vec = 0;
  var_el_etcone20_vec = 0;

  var_muons_vec = 0;

  var_muonPt_vec = 0;
  var_muonEta_vec = 0;
  var_muonPhi_vec = 0;
  var_mu_ptcone20_vec = 0;
  var_mu_etcone20_vec = 0;
  var_mu_charge_vec = 0;

} // initVEctors

/*
 * Setup the TLVs for the electrons and muons.  Pending an investigation by me (Tim), these will be created from the pt, eta and phi variables
 * for each of the leptons.
 */
void setLeptonVectors()
{
  if(!hvtllqq) // if not running hvtlqq exit
    return;

  int elsize = 0;

  if (var_electronPt_vec != NULL)
    elsize = var_electronPt_vec->size(); // number of electrons

  for (int x = 0; x < elsize; x++) // loop through all electrons
    {
      TLorentzVector e = TLorentzVector();
      // create electron 4 vector
      e.SetPtEtaPhiM((*var_electronPt_vec)[x], (*var_electronEta_vec)[x], (*var_electronPhi_vec)[x], ELMASS);
      electrons.push_back(e);
      
    }

  int musize = 0;
  // get number of muons
  if (var_muonPt_vec != NULL)
    musize = var_muonPt_vec->size();
  TLorentzVector m = TLorentzVector();
  for (int x = 0; x < musize; x++)
    {
      
      // create muon 4-vector
      m.SetPtEtaPhiM((*var_muonPt_vec)[x], (*var_muonEta_vec)[x], (*var_muonPhi_vec)[x], MUMASS);
      muons.push_back(m);
      
    }
  

} // setLeptonVectors()


/*
 * Set all of the output variables.  These are set to variables calculated per event or to variables read in. 
 * This also sets the weights for the output file.
 * 
 * @param jet_idx_truth The index of the truth jet being used.
 * @param jet_idx_topo The index of the topo jet being used.
 * @param jet_idx_groomed The index of the groomed jet being used.
 * @param subjet_idx The index of the groomed jet within the subjet collection.
 * @param groomalgo The abbreviated name of the algorithm.
 * @param samplePrefix If it has "LC" in the name of the algorithm
 */
void setOutputVariables( int jet_idx_truth, int jet_idx_topo, int jet_idx_groomed, int subjet_idx, std::string & groomalgo, std::string &  samplePrefix)
{
  // set to 0 for now
  int jet_idx = 0;
  // set some of the variables that are not collections/ vectors
  // mc_event_weight and avgIntPerXing are stored differently in xaod and d3pd
  if (xAOD)
    {
      mc_event_weight_out = mc_event_weight_xaod->at(0);// remove ->at(0) when running on D3PD
      avgIntpXingOut = avgIntpXingIn_xaod;
      actualIntPerXingOut = actualIntPerXingIn;
      evt_kfactor_out = evt_kfactor;
      evt_filtereff_out = evt_filtereff;
      evt_sumWeights_out = evt_sumWeights;
      evt_nEvts_out = evt_nEvts;
      evt_xsec_out = evt_xsec;
      
    }
  else
    {
      mc_event_weight_out = mc_event_weight_d3pd;// remove ->at(0) when running on D3PD
      avgIntpXingOut = avgIntpXingIn_d3pd;
    }

  mc_channel_number_out = mc_channel_number;
  runNumberOut = runNumberIn;
  nvtxOut = nvtxIn;

  scale1fbOut = scale1fb;

  bool addLC = false; // some algorithms have "LC" in their name


  if (!recluster) // if we're not doing reclustering
    addLC = true; // just add the LC to the name

  // default weights
  var_k_factor = 1.0;  
  var_filter_eff = 1.0;
  var_xs = 1.0;
  

  // keep a list of run numbers, check this list instead of checking each of the
  // k_factors, xs and filter_eff maps
  if (runNumber_map.find(mc_channel_number)!=runNumber_map.end())
    {
      var_k_factor = k_factors[mc_channel_number];
      var_filter_eff = filt_eff[mc_channel_number];
      var_xs = xs[mc_channel_number];
    }


  // loop through the different types of jets and set the output variables
  for (int x = 0; x < jetType::MAX ; x++)
    {
      switch (x)
	{
	case jetType::TRUTH:
	  jet_idx = jet_idx_truth;
	  break;
	case jetType::TOPO:
	  jet_idx = jet_idx_topo;
	  break;
	case jetType::GROOMED:
	  jet_idx = jet_idx_groomed;
	  break;
	}
      
      if (jet_idx == -99) // mostly just topo jets
	continue;
      
      // all of the _vec variables have a default value that is not null.
      var_E[x]=(*var_E_vec[x])[jet_idx];
      var_pt[x]=(*var_pt_vec[x])[jet_idx];
      var_m[x]=(*var_m_vec[x])[jet_idx];
      var_eta[x]=(*var_eta_vec[x])[jet_idx];
      var_phi[x]=(*var_phi_vec[x])[jet_idx];
      if (x!=0 && !xAOD && var_emfrac_vec[x] != NULL)//useBranch(string(jetString +"emfrac")))
	var_emfrac[x]=(*var_emfrac_vec[x])[jet_idx];
      var_Tau1[x]=(*var_Tau1_vec[x])[jet_idx];
      var_Tau2[x]=(*var_Tau2_vec[x])[jet_idx];

 
      var_SPLIT12[x]=sqrt((*var_SPLIT12_vec[x])[jet_idx]);

      var_Dip12[x]=(*var_Dip12_vec[x])[jet_idx];
      var_PlanarFlow[x]=(*var_PlanarFlow_vec[x])[jet_idx];
      var_Angularity[x]=(*var_Angularity_vec[x])[jet_idx];
      
      var_Aplanarity[x] = (*var_Aplanarity_vec[x])[jet_idx];
      var_Sphericity[x] = (*var_Sphericity_vec[x])[jet_idx];
      var_ThrustMaj[x] = (*var_ThrustMaj_vec[x])[jet_idx];
      var_ThrustMin[x] = (*var_ThrustMin_vec[x])[jet_idx];
      
      
      // tau21 and tauwta21 are set in the main loop, not here, because we have to calculate them
      var_TauWTA1[x]=(*var_TauWTA1_vec[x])[jet_idx];
      var_TauWTA2[x]=(*var_TauWTA2_vec[x])[jet_idx];
      var_ZCUT12[x]=(*var_ZCUT12_vec[x])[jet_idx];

      var_Mu12[x] = (*var_Mu12_vec[x])[jet_idx];
      // if the FoxWolfram20 and SoftDropTag variables are not being calculated, but exist in the input file, set them here
      // foxwolfram20 = foxwolfram2/foxwolfram0      
      if (!calcFoxWolfram20 && !preCalcFoxWolfram20)
	var_FoxWolfram20[x] = (*var_FoxWolfram2_vec[x])[jet_idx]/(*var_FoxWolfram0_vec[x])[jet_idx];
      else if (preCalcFoxWolfram20)
	var_FoxWolfram20[x] = (*var_FoxWolfram20_vec[x])[jet_idx];


      if (!calcSoftDrop)
	var_softdrop[x] = (*var_SoftDropTag_vec[x])[jet_idx];


    } // end for loop

      // yfilt only exists for groomed jets
  if (var_YFilt_vec != NULL)
    {

      // yfilt doesn't exist in the split/filtered samples, so it needs to be derived
      float calcyfilt = (*var_SPLIT12_vec[jetType::GROOMED])[jet_idx_groomed]/(*var_m_vec[jetType::GROOMED])[jet_idx_groomed];
      var_YFilt=sqrt(calcyfilt);

    }
 
  // only store this for groomed jets
  if (subjetscalc)
    {
      var_subjets_E = (*var_subjets_E_vec)[subjet_idx];
      var_subjets_pt = (*var_subjets_pt_vec)[subjet_idx];
      var_subjets_m = (*var_subjets_m_vec)[subjet_idx];
      var_subjets_eta = (*var_subjets_eta_vec)[subjet_idx];
      var_subjets_phi = (*var_subjets_phi_vec)[subjet_idx];
    }

} //setOutputVariables


/*
 * Calculate the response values for the different variables.
 */
void calculateResponseValues()
{
  // all of the _vec variables have a default value that is not null.
  int g_idx = jetType::GROOMED;
  int t_idx = jetType::TRUTH;

  response_E = var_E[g_idx]/var_E[t_idx];
  response_pt = var_pt[g_idx]/var_pt[t_idx];
  response_m = var_m[g_idx]/var_m[t_idx];
  response_eta = var_eta[g_idx]/var_eta[t_idx];
  response_phi = var_phi[g_idx]/var_phi[t_idx];
  response_Tau1 = var_Tau1[g_idx]/var_Tau1[t_idx];
  response_Tau2 = var_Tau2[g_idx]/var_Tau2[t_idx];
  response_SPLIT12 = var_SPLIT12[g_idx]/var_SPLIT12[t_idx];
  response_Dip12 = var_Dip12[g_idx]/var_Dip12[t_idx];
  response_PlanarFlow = var_PlanarFlow[g_idx]/var_PlanarFlow[t_idx];
  response_Angularity = var_Angularity[g_idx]/var_Angularity[t_idx];
  response_Aplanarity = var_Aplanarity[g_idx]/var_Aplanarity[t_idx];
  response_Sphericity = var_Sphericity[g_idx]/var_Sphericity[t_idx];
  response_ThrustMaj = var_ThrustMaj[g_idx]/var_ThrustMaj[t_idx];
  response_ThrustMin = var_ThrustMin[g_idx]/var_ThrustMin[t_idx];
  // tau21 and tauwta21 are set in the main loop, because they are calculated there
  response_TauWTA1 = var_TauWTA1[g_idx]/var_TauWTA1[t_idx];
  response_TauWTA2 = var_TauWTA2[g_idx]/var_TauWTA2[t_idx];
  response_ZCUT12 = var_ZCUT12[g_idx]/var_ZCUT12[t_idx];
  response_Mu12 = var_Mu12[g_idx]/var_Mu12[t_idx];
  response_FoxWolfram20 = var_FoxWolfram20[g_idx]/var_FoxWolfram20[t_idx];
  response_softdrop = var_softdrop[g_idx]/var_softdrop[t_idx];

} // calculateResponseValues



/*
 * Clear all of the variables used for the output file.
 */
void clearOutputVariables()
{

  var_massdrop = 0;
  var_yt = 0;
  var_isElectronEvent = 0;
  var_mllj = 0;
  var_mll = 0;
  var_ptll = 0;

  var_leadingJetPt = 0;
  var_YFilt = 0;

  // response values
  response_E = 0;
  response_pt = 0;
  response_m = 0;
  response_eta = 0;
  response_phi = 0;
  response_Tau1 = 0;
  response_Tau2 = 0;
  response_SPLIT12 = 0;
  response_Dip12 = 0;
  response_PlanarFlow = 0;
  response_Angularity = 0;
  response_Tau21 = 0;
  response_Mu12 = 0;
  
  // extra output responseiables
  response_TauWTA1 = 0; 
  response_TauWTA2 = 0; 
  response_TauWTA2TauWTA1 = 0; 
  response_ZCUT12 = 0;
  
  response_Aplanarity = 0;
  response_Sphericity = 0;
  response_ThrustMaj = 0;
  response_ThrustMin = 0;
  
  response_QjetVol = 0;
  response_FoxWolfram20 = 0;
  response_softdrop = 0;
  response_EEC_C2_1 = 0;
  response_EEC_C2_2 = 0;
  response_EEC_D2_1 = 0;
  response_EEC_D2_2 = 0;
  
  electrons.clear();
  muons.clear();
  var_E.clear();
  var_pt.clear();
  var_m.clear();
  var_eta.clear();
  var_phi.clear();
  var_emfrac.clear();
  var_Tau1.clear();
  var_Tau2.clear();

  var_SPLIT12.clear();

  var_Dip12.clear();
  var_PlanarFlow.clear();
  var_Angularity.clear();

  var_Aplanarity.clear();
  var_Sphericity.clear();
  var_ThrustMaj.clear();
  var_ThrustMin.clear();


  var_Tau21.clear();
  var_TauWTA2TauWTA1.clear();
  var_TauWTA1.clear();
  var_TauWTA2.clear();
  var_ZCUT12.clear();
  var_Mu12.clear();
  var_FoxWolfram20.clear();
  var_QjetVol.clear();
  var_softdrop.clear();
  var_EEC_C2_1.clear();
  var_EEC_C2_2.clear();
  var_EEC_D2_1.clear();
  var_EEC_D2_2.clear();

  var_leptons.clear();
  var_ptcone20.clear();
  var_etcone20.clear();
  var_charge.clear();
  

} // clearOutputVariables


/*
 * Reset all of the variables that are written to the output file.  All variables are set to a default value.
 */
void resetOutputVariables()
{
  clearOutputVariables();

  var_leptons.push_back(TLorentzVector());
  var_ptcone20.push_back(-999);
  var_etcone20.push_back(-999);
  var_charge.push_back(-999);
  
  for (int i = 0; i < jetType::MAX; i++)
    {
      var_E.push_back(-999);
      var_pt.push_back(-999);
      var_m.push_back(-999);
      var_eta.push_back(-999);
      var_phi.push_back(-999);
      var_emfrac.push_back(-999);
      var_Tau1.push_back(-999);
      var_Tau2.push_back(-999);

      var_SPLIT12.push_back(-999);
      var_Dip12.push_back(-999);
      var_PlanarFlow.push_back(-999);
      var_Angularity.push_back(-999);

      var_Aplanarity.push_back(-999);
      var_Sphericity.push_back(-999);
      var_ThrustMaj.push_back(-999);
      var_ThrustMin.push_back(-999);

      var_Tau21.push_back(-999);
      var_TauWTA2TauWTA1.push_back(-999);
      var_TauWTA1.push_back(-999);
      var_TauWTA2.push_back(-999);
      var_ZCUT12.push_back(-999);
      var_Mu12.push_back(-999);

      var_FoxWolfram20.push_back(-999);
      var_QjetVol.push_back(-999);
      var_softdrop.push_back(-999);
      var_EEC_C2_1.push_back(-999);
      var_EEC_C2_2.push_back(-999);
      var_EEC_D2_1.push_back(-999);
      var_EEC_D2_2.push_back(-999);


    }
} //resetOutputVariables

/*
 * Setup all of the branches for the jets and leptons in the output tree.
 *
 * @param tree A pointer to the output TTree.
 * @param groomalgo The shortened version of the grooming algorithm.
 * @param groomIdx The full grooming algorithm name, used as a key in the algorithms maps.
e */
void setOutputBranches(TTree * tree, std::string & groomalgo, std::string & groomIdx)
{

  std::string samplePrefix = ""; // AntiKt10 for example
  bool addLC = false; // some algorithms have the string "LC" in them.
  samplePrefix = algorithms.AlgoPrefix[groomIdx];

  if (algorithms.AlgoType[groomIdx].find("recluster") == std::string::npos) // we're doing reclustering
    addLC = true; // just add the LC to the name

  // these branches will always exist
  tree->Branch("mc_event_weight",&mc_event_weight_out,"mc_event_weight/F");
  tree->Branch("mc_channel_number", &mc_channel_number_out,"mc_channel_number/I");
  tree->Branch("vxp_n", &nvtxOut, "vxp_n/I");
  tree->Branch("averageIntPerXing",&avgIntpXingOut,"averageIntPerXing/F");
  tree->Branch("scale1fb",&scale1fbOut, "scale1fb/F");

  if (xAOD)
    {
      tree->Branch("evt_kfactor",&evt_kfactor_out,"evt_kfactor/F");
      tree->Branch("evt_nEvts",&evt_nEvts_out,"evt_nEvts/F");
      tree->Branch("evt_filtereff",&evt_filtereff_out,"evt_filtereff/F");
      tree->Branch("evt_sumWeights",&evt_sumWeights_out,"evt_sumWeights/F");
      tree->Branch("evt_xsec",&evt_xsec_out,"evt_xsec/F");
      tree->Branch("actualIntPerXing",&actualIntPerXingOut,"actualIntPerXing/F");
    }

  bool addResponse = true;

  for (int i = 0; i < jetType::MAX; i++) // truth, topo, groomed
    {
      if (i!=0)
	addResponse = false;
     
      std::string jetString = returnJetType(samplePrefix, groomalgo, addLC,i); //set to truth/ topo/ groomed

      // check that we actually want to use this branch/ output this branch
      if (useBranch(string(jetString+"E")))
	{
	  tree->Branch(std::string(jetString+"E").c_str(),&var_E.at(i),std::string(jetString+"E/F").c_str());
	  // add the response branch, but only need one
	  if (addResponse)
	    {
	      tree->Branch("response_E",&response_E,"response_E/F");
	    }
	}
      if (useBranch(string(jetString+"pt")))
	{
	  tree->Branch(std::string(jetString+"pt").c_str(),&var_pt.at(i),std::string(jetString+"pt/F").c_str());
	  // add the response branch, but only need one
	  if (addResponse)
	    {
	      tree->Branch("response_pt",&response_pt,"response_pt/F");
	    }
	}
      if (useBranch(string(jetString+"m")))
	{
	  tree->Branch(std::string(jetString+"m").c_str(),&var_m.at(i),std::string(jetString+"m/F").c_str());
	  // add the response branch, but only need one
	  if (addResponse)
	    {
	      tree->Branch("response_m",&response_m,"response_m/F");
	    }
	}
      if (useBranch(string(jetString+"eta")))
	{
	  tree->Branch(std::string(jetString+"eta").c_str(),&var_eta.at(i),std::string(jetString+"eta/F").c_str());
	  // add the response branch, but only need one
	  if (addResponse)
	    {
	      tree->Branch("response_eta",&response_eta,"response_eta/F");
	    }
	}
      if (useBranch(string(jetString+"phi")))
	{
	  tree->Branch(std::string(jetString+"phi").c_str(),&var_phi.at(i),std::string(jetString+"phi/F").c_str());
	  // add the response branch, but only need one
	  if (addResponse)
	    {
	      tree->Branch("response_phi",&response_phi,"response_phi/F");
	    }
	}
      if (i != jetType::TRUTH && useBranch(string(jetString+"emfrac"))) // emfrac doesn't exist for truth jets
	{
	  tree->Branch(std::string(jetString+"emfrac").c_str(),&var_emfrac.at(i),std::string(jetString+"emfrac"+"/F").c_str());
	}
      if (useBranch(string(jetString+"Tau1")))
	{
	  tree->Branch(std::string(jetString+"Tau1").c_str(),&var_Tau1.at(i),std::string(jetString+"Tau1/F").c_str());
	  // add the response branch, but only need one
	  if (addResponse)
	    {
	      tree->Branch("response_Tau1",&response_Tau1,"response_Tau1/F");
	    }
	}
      if (useBranch(string(jetString+"Tau2")))
<<<<<<< HEAD
      tree->Branch(std::string(jetString+"Tau2").c_str(),&var_Tau2.at(i),std::string(jetString+"Tau2/F").c_str());
<<<<<<< HEAD
      if (useBranch(string(jetString+"Tau3")))
      tree->Branch(std::string(jetString+"Tau3").c_str(),&var_Tau3.at(i),std::string(jetString+"Tau3/F").c_str());
      if (useBranch(string(jetString+"WIDTH")))
      tree->Branch(std::string(jetString+"WIDTH").c_str(),&var_WIDTH.at(i),std::string(jetString+"WIDTH/F").c_str());
      if (useBranch(string(jetString+"SPLIT12")))
      tree->Branch(std::string(jetString+"SPLIT12").c_str(),&var_SPLIT12.at(i),std::string(jetString+"SPLIT12/F").c_str());
      if (useBranch(string(jetString+"SPLIT23")))
      tree->Branch(std::string(jetString+"SPLIT23").c_str(),&var_SPLIT23.at(i),std::string(jetString+"SPLIT23/F").c_str());
      if (useBranch(string(jetString+"SPLIT34")))
      tree->Branch(std::string(jetString+"SPLIT34").c_str(),&var_SPLIT34.at(i),std::string(jetString+"SPLIT34/F").c_str());
      if (useBranch(string(jetString+"Dip12")))
      tree->Branch(std::string(jetString+"Dip12").c_str(),&var_Dip12.at(i),std::string(jetString+"Dip12/F").c_str());
      if (useBranch(string(jetString+"Dip13")))
      tree->Branch(std::string(jetString+"Dip13").c_str(),&var_Dip13.at(i),std::string(jetString+"Dip13/F").c_str());
      if (useBranch(string(jetString+"Dip23")))
      tree->Branch(std::string(jetString+"Dip23").c_str(),&var_Dip23.at(i),std::string(jetString+"Dip23/F").c_str());
      if (useBranch(string(jetString+"DipExcl12")))
      tree->Branch(std::string(jetString+"DipExcl12").c_str(),&var_DipExcl12.at(i),std::string(jetString+"DipExcl12/F").c_str());
      if (useBranch(string(jetString+"PlanarFlow")))
=======

<<<<<<< HEAD
      tree->Branch(std::string(jetString+"SPLIT12").c_str(),&var_SPLIT12.at(i),std::string(jetString+"SPLIT12/F").c_str());
      tree->Branch(std::string(jetString+"Dip12").c_str(),&var_Dip12.at(i),std::string(jetString+"Dip12/F").c_str());
>>>>>>> removed variables
      tree->Branch(std::string(jetString+"PlanarFlow").c_str(),&var_PlanarFlow.at(i),std::string(jetString+"PlanarFlow/F").c_str());
      if (useBranch(string(jetString+"Angularity")))
      tree->Branch(std::string(jetString+"Angularity").c_str(),&var_Angularity.at(i),std::string(jetString+"Angularity/F").c_str());
<<<<<<< HEAD
      if (useBranch(string(jetString+"QW")))
      tree->Branch(std::string(jetString+"QW").c_str(),&var_QW.at(i),std::string(jetString+"QW/F").c_str());
      if (useBranch(string(jetString+"PullMag")))
      tree->Branch(std::string(jetString+"PullMag").c_str(),&var_PullMag.at(i),std::string(jetString+"PullMag/F").c_str());
      if (useBranch(string(jetString+"PullPhi")))
      tree->Branch(std::string(jetString+"PullPhi").c_str(),&var_PullPhi.at(i),std::string(jetString+"PullPhi/F").c_str());
      if (useBranch(string(jetString+"Pull_C00")))
      tree->Branch(std::string(jetString+"Pull_C00").c_str(),&var_Pull_C00.at(i),std::string(jetString+"Pull_C00/F").c_str());
      if (useBranch(string(jetString+"Pull_C01")))
      tree->Branch(std::string(jetString+"Pull_C01").c_str(),&var_Pull_C01.at(i),std::string(jetString+"Pull_C01/F").c_str());
      if (useBranch(string(jetString+"Pull_C10")))
      tree->Branch(std::string(jetString+"Pull_C10").c_str(),&var_Pull_C10.at(i),std::string(jetString+"Pull_C10/F").c_str());
      if (useBranch(string(jetString+"Pull_C11")))
      tree->Branch(std::string(jetString+"Pull_C11").c_str(),&var_Pull_C11.at(i),std::string(jetString+"Pull_C11/F").c_str());
      if (useBranch(string(jetString+"YFilt")))
      tree->Branch(std::string(jetString+"YFilt").c_str(),&var_YFilt,std::string(jetString+"YFilt/F").c_str());
=======

>>>>>>> Fixed the emfrac bug for llqq samples

      if (useBranch(string(jetString+"ActiveArea")))
	tree->Branch(std::string(jetString+"ActiveArea").c_str(),&var_ActiveArea,std::string(jetString+"ActiveArea/F").c_str());
      if (useBranch(string(jetString+"VoronoiArea")))
	tree->Branch(std::string(jetString+"VoronoiArea").c_str(),&var_VoronoiArea,std::string(jetString+"VoronoiArea/F").c_str());
=======
	{
	  tree->Branch(std::string(jetString+"Tau2").c_str(),&var_Tau2.at(i),std::string(jetString+"Tau2/F").c_str());
	  // add the response branch, but only need one
	  if (addResponse)
	    {
	      tree->Branch("response_Tau2",&response_Tau2,"response_Tau2/F");
	    }
	}

      if (useBranch(string(jetString+"SPLIT12")))
	{
	  tree->Branch(std::string(jetString+"SPLIT12").c_str(),&var_SPLIT12.at(i),std::string(jetString+"SPLIT12/F").c_str());
	  // add the response branch, but only need one
	  if (addResponse)
	    {
	      tree->Branch("response_SPLIT12",&response_SPLIT12,"response_SPLIT12/F");
	    }
	}
      if (useBranch(string(jetString+"Dip12")))
	{
	  tree->Branch(std::string(jetString+"Dip12").c_str(),&var_Dip12.at(i),std::string(jetString+"Dip12/F").c_str());
	  // add the response branch, but only need one
	  if (addResponse)
	    {
	      tree->Branch("response_Dip12",&response_Dip12,"response_Dip12/F");
	    }
	}

      if (useBranch(string(jetString+"PlanarFlow")))
	{
	  tree->Branch(std::string(jetString+"PlanarFlow").c_str(),&var_PlanarFlow.at(i),std::string(jetString+"PlanarFlow/F").c_str());
	  // add the response branch, but only need one
	  if (addResponse)
	    {
	      tree->Branch("response_PlanarFlow",&response_PlanarFlow,"response_PlanarFlow/F");
	    }
	}
      if (useBranch(string(jetString+"Angularity")))
	{
	  tree->Branch(std::string(jetString+"Angularity").c_str(),&var_Angularity.at(i),std::string(jetString+"Angularity/F").c_str());
	  // add the response branch, but only need one
	  if (addResponse)
	    {
	      tree->Branch("response_Angularity",&response_Angularity,"response_Angularity/F");
	    }
	}
	  
>>>>>>> Added the scale factors that are calculated in the post processing of the xAODs.
      if (useBranch(string(jetString+"Aplanarity")))
	{
	  tree->Branch(std::string(jetString+"Aplanarity").c_str(),&var_Aplanarity.at(i),std::string(jetString+"Aplanarity/F").c_str());
	  // add the response branch, but only need one
	  if (addResponse)
	    {
	      tree->Branch("response_Aplanarity",&response_Aplanarity,"response_Aplanarity/F");
	    }
	}
      if (useBranch(string(jetString+"Sphericity")))
	{
	  tree->Branch(std::string(jetString+"Sphericity").c_str(),&var_Sphericity.at(i),std::string(jetString+"Sphericity/F").c_str());
		  // add the response branch, but only need one
	  if (addResponse)
	    {
	      tree->Branch("response_Sphericity",&response_Sphericity,"response_Sphericity/F");
	    }
	}
      if (useBranch(string(jetString+"ThrustMaj")))
	{
	  tree->Branch(std::string(jetString+"ThrustMaj").c_str(),&var_ThrustMaj.at(i),std::string(jetString+"ThrustMaj/F").c_str());
	  // add the response branch, but only need one
	  if (addResponse)
	    {
	      tree->Branch("response_ThrustMaj",&response_ThrustMaj,"response_ThrustMaj/F");
	    }
	}
      if (useBranch(string(jetString+"ThrustMin")))
	{
	  tree->Branch(std::string(jetString+"ThrustMin").c_str(),&var_ThrustMin.at(i),std::string(jetString+"ThrustMin/F").c_str());
	  // add the response branch, but only need one
	  if (addResponse)
	    {
	      tree->Branch("response_ThrustMin",&response_ThrustMin,"response_ThrustMin/F");
	    }
	}
      
      if (useBranch(string(jetString+"TauWTA1")))
	{
	  tree->Branch(std::string(jetString+"TauWTA1").c_str(),&var_TauWTA1.at(i),std::string(jetString+"TauWTA1/F").c_str());
	  // add the response branch, but only need one
	  if (addResponse)
	    {
	      tree->Branch("response_TauWTA1",&response_TauWTA1,"response_TauWTA1/F");
	    }
	}
      if (useBranch(string(jetString+"TauWTA2")))
<<<<<<< HEAD
	tree->Branch(std::string(jetString+"TauWTA2").c_str(),&var_TauWTA2.at(i),std::string(jetString+"TauWTA2/F").c_str());
      if (useBranch(string(jetString+"TauWTA3")))
	tree->Branch(std::string(jetString+"TauWTA3").c_str(),&var_TauWTA3.at(i),std::string(jetString+"TauWTA3/F").c_str());
      if (useBranch(string(jetString+"ZCUT12")))
	tree->Branch(std::string(jetString+"ZCUT12").c_str(),&var_ZCUT12.at(i),std::string(jetString+"ZCUT12/F").c_str());
<<<<<<< HEAD
      if (useBranch(string(jetString+"ZCUT23")))
	tree->Branch(std::string(jetString+"ZCUT23").c_str(),&var_ZCUT23.at(i),std::string(jetString+"ZCUT23/F").c_str());
      if (useBranch(string(jetString+"ZCUT34")))
	tree->Branch(std::string(jetString+"ZCUT34").c_str(),&var_ZCUT34.at(i),std::string(jetString+"ZCUT34/F").c_str());
=======
      if (useBranch(string(jetString+"Mu12")))
	tree->Branch(std::string(jetString+"Mu12").c_str(),&var_Mu12.at(i),std::string(jetString+"Mu12/F").c_str());
>>>>>>> Added new variables to run ECF on xAOD.  Some things are broken for 8 tev now like mc_event_weight
      
<<<<<<< HEAD
<<<<<<< HEAD
      if (useBranch(string(returnJetType(samplePrefix, groomalgo, addLC,0)+"TauWTA2/TauWTA1")))
	tree->Branch(std::string(returnJetType(samplePrefix, groomalgo, addLC,0)+"TauWTA2/TauWTA1").c_str(),&var_TauWTA2TauWTA1.at(0),std::string(jetString+"TauWTA2TauWTA1/F").c_str());
      if (useBranch(string(returnJetType(samplePrefix, groomalgo, addLC,1)+"TauWTA2/TauWTA1")))
	tree->Branch(std::string(returnJetType(samplePrefix, groomalgo, addLC,1)+"TauWTA2/TauWTA1").c_str(),&var_TauWTA2TauWTA1.at(1),std::string(jetString+"TauWTA2TauWTA1/F").c_str());
      if (useBranch(string(jetString+"TauWTA2/TauWTA1")))
	tree->Branch(std::string(jetString+"TauWTA2/TauWTA1").c_str(),&var_TauWTA2TauWTA1.at(2),std::string(jetString+"TauWTA2TauWTA1/F").c_str());  
<<<<<<< HEAD
      
=======
      tree->Branch(std::string(jetString+"YFilt").c_str(),&var_YFilt,std::string(jetString+"YFilt/F").c_str());


      tree->Branch(std::string(jetString+"Aplanarity").c_str(),&var_Aplanarity,std::string(jetString+"Aplanarity/F").c_str());
      tree->Branch(std::string(jetString+"Sphericity").c_str(),&var_Sphericity,std::string(jetString+"Sphericity/F").c_str());
      tree->Branch(std::string(jetString+"ThrustMaj").c_str(),&var_ThrustMaj,std::string(jetString+"ThrustMaj/F").c_str());
      tree->Branch(std::string(jetString+"ThrustMin").c_str(),&var_ThrustMin,std::string(jetString+"ThrustMin/F").c_str());
      
=======
=======
      if (useBranch(string(returnJetType(samplePrefix, groomalgo, addLC,0)+"TauWTA2TauWTA1")))
	tree->Branch(std::string(returnJetType(samplePrefix, groomalgo, addLC,0)+"TauWTA2TauWTA1").c_str(),&var_TauWTA2TauWTA1.at(0),std::string(jetString+"TauWTA2TauWTA1/F").c_str());
      if (useBranch(string(returnJetType(samplePrefix, groomalgo, addLC,1)+"TauWTA2TauWTA1")))
	tree->Branch(std::string(returnJetType(samplePrefix, groomalgo, addLC,1)+"TauWTA2TauWTA1").c_str(),&var_TauWTA2TauWTA1.at(1),std::string(jetString+"TauWTA2TauWTA1/F").c_str());
      if (useBranch(string(jetString+"TauWTA2TauWTA1")))
	tree->Branch(std::string(jetString+"TauWTA2TauWTA1").c_str(),&var_TauWTA2TauWTA1.at(2),std::string(jetString+"TauWTA2TauWTA1/F").c_str());  
>>>>>>> Added FoxWolfram and SoftDropTag from the D3PD as another option instead of calculating them by hand.
=======
>>>>>>> Fixed the emfrac bug for llqq samples
=======
	{
	  tree->Branch(std::string(jetString+"TauWTA2").c_str(),&var_TauWTA2.at(i),std::string(jetString+"TauWTA2/F").c_str());
	  // add the response branch, but only need one
	  if (addResponse)
	    {
	      tree->Branch("response_TauWTA2",&response_TauWTA2,"response_TauWTA2/F");
	    }
	}

      if (useBranch(string(jetString+"ZCUT12")))
	{
	  tree->Branch(std::string(jetString+"ZCUT12").c_str(),&var_ZCUT12.at(i),std::string(jetString+"ZCUT12/F").c_str());
	  // add the response branch, but only need one
	  if (addResponse)
	    {
	      tree->Branch("response_ZCUT12",&response_ZCUT12,"response_ZCUT12/F");
	    }
	}
      if (useBranch(string(jetString+"Mu12")))
	{
	  tree->Branch(std::string(jetString+"Mu12").c_str(),&var_Mu12.at(i),std::string(jetString+"Mu12/F").c_str());
      	  // add the response branch, but only need one
	  if (addResponse)
	    {
	      tree->Branch("response_Mu12",&response_Mu12,"response_Mu12/F");
	    }
	}
>>>>>>> Added the scale factors that are calculated in the post processing of the xAODs.


      if (useBranch(string(jetString+"QJetsVol")))
	{
	  tree->Branch(std::string(jetString+"QJetsVol").c_str(), &var_QjetVol.at(i), std::string(jetString+"QJetsVol").c_str());
	  // add the response branch, but only need one
	  if (addResponse)
	    {
	      tree->Branch("response_QJetsVol",&response_QJetsVol,"response_QJetsVol/F");
	    }
	}
      if (useBranch(string(jetString+"FoxWolfram20")))
	{
	  tree->Branch(std::string(jetString+"FoxWolfram20").c_str(), &var_FoxWolfram20.at(i), std::string(jetString+"FoxWolfram20").c_str());
	  // add the response branch, but only need one
	  if (addResponse)
	    {
	      tree->Branch("response_FoxWolfram20",&response_FoxWolfram20,"response_FoxWolfram20/F");
	    }
	}
	
      if (useBranch(string(jetString+"SoftDrop")))
<<<<<<< HEAD
	tree->Branch(std::string(jetString+"SoftDrop").c_str(), &var_softdrop.at(i), std::string(jetString+"SoftDrop").c_str());
<<<<<<< HEAD
      if (useBranch(string(jetString+"EEC_C1")))
	tree->Branch(std::string(jetString+"EEC_C1").c_str(), &var_EEC_C1.at(i), std::string(jetString+"EEC_C1").c_str());
      if (useBranch(string(jetString+"EEC_C2")))
	tree->Branch(std::string(jetString+"EEC_C2").c_str(), &var_EEC_C2.at(i), std::string(jetString+"EEC_C2").c_str());
      if (useBranch(string(jetString+"EEC_D1")))
	tree->Branch(std::string(jetString+"EEC_D1").c_str(), &var_EEC_D1.at(i), std::string(jetString+"EEC_D1").c_str());
      if (useBranch(string(jetString+"EEC_D2")))
	tree->Branch(std::string(jetString+"EEC_D2").c_str(), &var_EEC_D2.at(i), std::string(jetString+"EEC_D2").c_str());
>>>>>>> Added qjets vol, softdrop tag, eec, fox wolfram.  Has not yet been tested, but it does compile.  QjetsPlugin is in a sep folder and should be added to the git repo.

      tree->Branch(std::string(jetString+"TauWTA1").c_str(),&var_TauWTA1.at(i),std::string(jetString+"TauWTA1/F").c_str());
      tree->Branch(std::string(jetString+"TauWTA2").c_str(),&var_TauWTA2.at(i),std::string(jetString+"TauWTA2/F").c_str());
      tree->Branch(std::string(jetString+"ZCUT12").c_str(),&var_ZCUT12.at(i),std::string(jetString+"ZCUT12/F").c_str());

	  
      tree->Branch(std::string(returnJetType(samplePrefix, groomalgo, addLC,0)+"TauWTA2/TauWTA1").c_str(),&var_TauWTA2TauWTA1.at(0),std::string(jetString+"TauWTA2TauWTA1/F").c_str());
      tree->Branch(std::string(returnJetType(samplePrefix, groomalgo, addLC,1)+"TauWTA2/TauWTA1").c_str(),&var_TauWTA2TauWTA1.at(1),std::string(jetString+"TauWTA2TauWTA1/F").c_str());
      tree->Branch(std::string(jetString+"TauWTA2/TauWTA1").c_str(),&var_TauWTA2TauWTA1.at(2),std::string(jetString+"TauWTA2TauWTA1/F").c_str());  
	  
>>>>>>> removed variables
      // add a calculated variable Tau2/Tau1
      if (useBranch(string(returnJetType(samplePrefix, groomalgo, addLC,0)+"Tau21")))
	tree->Branch(std::string(returnJetType(samplePrefix, groomalgo, addLC,0)+"Tau21").c_str(),&var_Tau21.at(0),std::string(jetString+"Tau21/F").c_str());
      if (useBranch(string(returnJetType(samplePrefix, groomalgo, addLC,0)+"Tau21")))
	tree->Branch(std::string(returnJetType(samplePrefix, groomalgo, addLC,1)+"Tau21").c_str(),&var_Tau21.at(1),std::string(jetString+"Tau21/F").c_str());
      if (useBranch(string(jetString+"Tau21")))
	tree->Branch(std::string(jetString+"Tau21").c_str(),&var_Tau21.at(2),std::string(jetString+"Tau21/F").c_str());  
=======
=======
	{
	  tree->Branch(std::string(jetString+"SoftDrop").c_str(), &var_softdrop.at(i), std::string(jetString+"SoftDrop").c_str());
	  // add the response branch, but only need one
	  if (addResponse)
	    {
	      tree->Branch("response_SoftDrop",&response_SoftDrop,"response_SoftDrop/F");
	    }
	}
>>>>>>> Added the scale factors that are calculated in the post processing of the xAODs.
      if (useBranch(string(jetString+"EEC_C2_1")))
	{
	  tree->Branch(std::string(jetString+"EEC_C2_1").c_str(), &var_EEC_C2_1.at(i), std::string(jetString+"EEC_C2_1").c_str());
	  // add the response branch, but only need one
	  if (addResponse)
	    {
	      tree->Branch("response_EEC_C2_1",&response_EEC_C2_1,"response_EEC_C2_1/F");
	    }
	}
      if (useBranch(string(jetString+"EEC_C2_2")))
	{
	  tree->Branch(std::string(jetString+"EEC_C2_2").c_str(), &var_EEC_C2_2.at(i), std::string(jetString+"EEC_C2_2").c_str());
	  // add the response branch, but only need one
	  if (addResponse)
	    {
	      tree->Branch("response_EEC_C2_2",&response_EEC_C2_2,"response_EEC_C2_2/F");
	    }
	}
      if (useBranch(string(jetString+"EEC_D2_1")))
	{
	  tree->Branch(std::string(jetString+"EEC_D2_1").c_str(), &var_EEC_D2_1.at(i), std::string(jetString+"EEC_D2_1").c_str());
	  // add the response branch, but only need one
	  if (addResponse)
	    {
	      tree->Branch("response_EEC_D2_1",&response_EEC_D2_1,"response_EEC_D2_1/F");
	    }
	}
      if (useBranch(string(jetString+"EEC_D2_2")))
	{
	  tree->Branch(std::string(jetString+"EEC_D2_2").c_str(), &var_EEC_D2_2.at(i), std::string(jetString+"EEC_D2_2").c_str());
	  // add the response branch, but only need one
	  if (addResponse)
	    {
	      tree->Branch("response_EEC_D2_2",&response_EEC_D2_2,"response_EEC_D2_2/F");
	    }
	}

>>>>>>> Fixed the emfrac bug for llqq samples

    }
  for (int wta_idx = 0; wta_idx < jetType::MAX; wta_idx++)
    {
      string jetstring = returnJetType(samplePrefix, groomalgo, addLC, wta_idx);
      bool addResponse = wta_idx == 0 ? true : false;
      if (useBranch(string(jetstring+"TauWTA2TauWTA1")))
	{
	  tree->Branch(std::string(jetstring+"TauWTA2TauWTA1").c_str(),&var_TauWTA2TauWTA1.at(wta_idx),std::string(jetstring+"TauWTA2TauWTA1/F").c_str());
	  // add the response branch, but only need one
	  if (addResponse)
	    {
	      tree->Branch("response_TauWTA2TauWTA1",&response_TauWTA2TauWTA1,"response_TauWTA2TauWTA1/F");
	    }
	}
    }

  // add a calculated variable Tau2/Tau1
  for (int tau_idx = 0 ; tau_idx < jetType::MAX; tau_idx++)
    {
      string jetstring = returnJetType(samplePrefix, groomalgo, addLC,tau_idx);
      bool addResponse = wta_idx == 0 ? true : false;
      if (useBranch(string(jetstring+"Tau21")))
	{
	  tree->Branch(std::string(jetstring+"Tau21").c_str(),&var_Tau21.at(tau_idx),std::string(jetstring+"Tau21/F").c_str());
	  if (addResponse)
	    {
	      tree->Branch("response_Tau21",&response_Tau21,"response_Tau21/F");
	    }
	}

    }

  string jetstring = returnJetType(samplePrefix, groomalgo, addLC, jetType::GROOMED);
  // this only exists for the groomed jets
  if (useBranch(string(jetstring+"YFilt")))
    tree->Branch(std::string(jetstring+"YFilt").c_str(),&var_YFilt,std::string(jetstring+"YFilt/F").c_str());
  

  std::string jetString = returnJetType(samplePrefix, groomalgo, addLC,jetType::GROOMED); 
  // add the lepton branches if running hvtllqq analysis
  if (hvtllqq)
    addLeptonBranches(jetString, tree);

} // setOutputBranches



/*
 * Check if the branch exists in a tree and then set the address if it does.
 *
 * @param tree A TChain pointer to the input file.
 * @param name Name of the branch.
 * @param var_vec A pointer to a vector of Float_t that is the variable used to read in from the tree.
 */
void setAddress(TChain * tree, std::string name, std::vector<Float_t> * var_vec)
{
  if(tree->GetListOfBranches()->FindObject(name.c_str()) )
    tree->SetBranchAddress(name.c_str(),&var_vec);
  else
    std::cout << "branch not found: " << name << std::endl;
  
} //setAddress


/*
 * read in the weights for the different samples.  Weights are based on RunNumbers, the file should have
 * RunNumber,SampleName,xs,k-factor,filter-efficiency
 * weights come from these two twiki pages: https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/BosonTaggingPaper2012#7th_July_2014_Mixture_with_small and https://twiki.cern.ch/twiki/bin/view/AtlasProtected/BoostedBosonTaggingD3PD#Properly_normalizing_MC_samples
 */
void readWeights()
{
  // open input file
  std::cout << "using weights file: " << weightsfile << std::endl;
  ifstream wf (weightsfile);
  string line;
  std::cout << "reading weights" << std::endl;
  // read all lines from file
  while (getline( wf, line))
    {
      trim(line); // remove trailing whitespace
      std::vector<std::string> strs; // use this to hold the vector obtained by splitting the input line
      boost::split(strs, line, boost::is_any_of(","));

      int cntr = 0;
      long runNumber = 0;
      // iterate through the elements in strs vector/ tokens in the input line
      for (std::vector<std::string>::iterator it = strs.begin(); it != strs.end(); it ++)
	{
	  switch (cntr)
	    {
	      // runnumber
	    case 0:
	      try
		{
		  runNumber = std::stol(*it);
		  runNumber_map[runNumber] = 1;
		}
	      catch (exception &e)
		{
		  std::cout << "error converting to run number:" << *it << std::endl;
		}
	      k_factors[runNumber] = 1;
	      filt_eff[runNumber] = 1;
	      xs[runNumber] = 1;
	      break;
	      // This will be the sample name, but not using it right now.
	    case 1:
	      // This will be the sample name, but not using it right now.
	      break;
	      // cross section
	    case 2:
	      try
		{
		  xs[runNumber] = std::stof(*it);
		}
	      catch (exception &e)
		{
		  std::cout << "error converting to xs:" << *it << std::endl;
		}
	      break;
	      // k factor
	    case 3:
	      try
		{
		  k_factors[runNumber] = std::stof(*it);
		}
	      catch (exception &e)
		{
		  std::cout << "error converting to k-factor:" << *it << std::endl;
		}
	      break;
	      // filter efficiency
	    case 4:
	      try
		{
		  filt_eff[runNumber] = std::stof(*it);
		}
	      catch (exception &e)
		{
		  std::cout << "error converting to eff:" << *it << std::endl;
		}
	      break;
	    default:
	      // this is probably an extra , on the end or something and is probably not needed!
	      break;
	      
	    }
	  cntr++;
	}
    }

  wf.close();
  
} // readWeights


/* 
 * Create an output file containing a number of datapoints for each algorithm + pt range.
 * each datapoint contains bkg/signal for each bin in bkg and sig arguments
 * this output file will be read in by the plotting code later and used for reweighting the signal sample
 * 
 * @param bkg Pointer to background TH1F.
 * @param sig Pointer to signal TH1F.
 * @param fname Name of output reweight file.
 */
void createPtReweightFile(TH1F * bkg, TH1F * sig, std::string & fname)
{
  // number of bins
  int bins = bkg->GetNbinsX();
  // open output file
  ofstream out(fname);
  double weight;

  for (int b = 1; b <= bins; b++) // loop through all bins
    {
      weight = bkg->GetBinContent(b);
      // don't divide by 0
      if (sig->GetBinContent(b) != 0)
	weight/=sig->GetBinContent(b);
      else
	{
	  std::cout << "oh no, we've been rumbled!  We have no signal in this bin and now we have to reweight by over 9000! But actually we're just going to make it 0." << std::endl;
	  weight = 0;
	}
      // write out low edge, bkg value, signal value
      out << bkg->GetXaxis()->GetBinLowEdge(b) << "," << bkg->GetBinContent(b) << "," << sig->GetBinContent(b) << endl;     
    }
  out.close();  
} // end createPtReweightHistogram


/*
 * Loads algorithm config from algorithms.xml.  It sets up maps for each algorithm for the abbreviation, the plot label, type of 
 * algorithm (pruned/trimmed/splitfiltered), the subjet grooming algorithm name (sometimes different), subjet index branch and bin
 * label.
 *
 * @param filename The filename of the configuration file.
 */
void Algorithms::load(const std::string & filename)
{
  // Create an empty property tree object
  using boost::property_tree::ptree;
  ptree pt;

  // read in xml
  read_xml(filename,pt);
  // loop through all of the xml entries
  BOOST_FOREACH( ptree::value_type & v, pt.get_child("config.algorithms"))
    {
      // name of algorithm - CamKt12topolcBLAH
      std::string name =  v.second.get<std::string>("Algorithm","");
      if (name != "")
	{
	  // The name without CamKtX/ AntiKtX
	  std::string groomalgo = v.second.get<std::string>("GroomingAlgorithm","");
	  if (groomalgo == "")
	    {
	      groomalgo = name.substr(name.find("Topo"),name.length()-1); // shortened version of name, - AntiKt10 for eg.
	    }
	  AlgoNames[name] = groomalgo;
	  // the jet algorithm being used - CamKtX for example
	  std::string jetalgo = v.second.get<std::string>("JetAlgorithm","");
	  if (jetalgo == "")
	    {
	      // up until Topo in string is the algorithm
	      jetalgo = name.substr(0,name.find("Topo"));
	      if (jetalgo.find("LC") != std::string::npos)
		jetalgo.erase(jetalgo.end()-2,jetalgo.end()); // erase LC from the name
	    }

	  // set all of the remaining properties
	  AlgoPrefix[name] = jetalgo;
	  AlgoList[name] = v.second.get<std::string>("Abbreviation","");
	  AlgoListN[name] = v.second.get<std::string>("PlotLabel","");
	  AlgoType[name] = v.second.get<std::string>("Type","NONE");
	  // subjet string
	  subjetMap[name] = v.second.get<std::string>("SubjetGroomingAlgorithm");
	  subjetIndex[name] = v.second.get<std::string>("SubjetIndexBranch");
	  binLabel[name] = v.second.get<std::string>("BinLabel","");
	}
    };

} // algorithms::load()



/*
 * Calculate the Energy correction using formula EEC = sum1*sum2/(sum3^2)
 * Sum1 = Sum(ijk) : pt(i)*pt(j)*pt(k) * [dR(ij)dR(ik)dR(jk) ]^beta
 * Sum2 = Sum(i) : pt(i)
 * Sum3 = Sum(ij) : pt(i)*pt(j)* [dR(ij)]^beta
 *
 * @param jettype The jet type - truth, topo or groomed
 * @param beta The beta value which is the power that dR is raised to in Sum1.
 * @param exp The power that Sum3 is raised to.
 * @return double containing eec
 */
double calculateEEC(int jettype, float beta, float exp)
{
  double Sum1=0;
  
  int size = (*var_pt_vec[jettype]).size();
  // calculate sum1
  for(int i=0; i<size; i++){
    for(int j=i+1; j<size; j++){
      for(int k=j+1; k<size; k++){
	// calculate deltaR values
	float dr_ij = DeltaR((*var_eta_vec[jettype])[i],(*var_phi_vec[jettype])[i],(*var_eta_vec[jettype])[j],(*var_phi_vec[jettype])[j]);
	float dr_ik = DeltaR((*var_eta_vec[jettype])[i],(*var_phi_vec[jettype])[i],(*var_eta_vec[jettype])[k],(*var_phi_vec[jettype])[k]);
	float dr_jk = DeltaR((*var_eta_vec[jettype])[j],(*var_phi_vec[jettype])[j],(*var_eta_vec[jettype])[k],(*var_phi_vec[jettype])[k]);
	// increment sum1
	Sum1 += ((*var_pt_vec[jettype])[i]/1000.)*((*var_pt_vec[jettype])[j]/1000.)*((*var_pt_vec[jettype])[k]/1000.) * pow( dr_ij * dr_ik * dr_jk, beta );
      }
    }
  }

  // sum2
  double Sum2=0;
  for(int i=0; i<size; i++){
    Sum2 += (*var_pt_vec[jettype])[i]/1000.;
  }

  //sum3
  double Sum3=0;
  for(int i=0; i<size; i++){
    for(int j=i+1; j<size; j++){
      float dr_ij = DeltaR((*var_eta_vec[jettype])[i],(*var_phi_vec[jettype])[i],(*var_eta_vec[jettype])[j],(*var_phi_vec[jettype])[j]);
      Sum3 += ((*var_pt_vec[jettype])[i]/1000.)*((*var_pt_vec[jettype])[j]/1000.) * pow( dr_ij, beta );
    }
  }

  if (Sum3 == 0) // don't want to divide by zero
    return 0;
  return Sum1*Sum2/pow(Sum3, exp);

}//calculateEEC


/*
 * Calculate the Energy correction factors
 * ECF1 = Sum(i) : pt(i)
 * ECF2 = Sum(ij) : pt(i)*pt(j)* [dR(ij)]^beta
 * ECF3 = Sum(ijk) : pt(i)*pt(j)*pt(k) * [dR(ij)dR(ik)dR(jk) ]^beta
 *
 * @param cons A vector<TLV> of the clusters for each jet
 * @param jettype The jet type - truth, topo or groomed
 * @param beta The beta value which is the power that dR is raised to in Sum1.
 * @param exp The power that Sum3 is raised to.
 * @return double containing eec
 */
void calculateECF(vector<TLorentzVector> & cons, int jettype, float beta)
{
  //double beta=0.3;
  // we store two beta values, so just getting the idx right to store it in the vector
  int beta_idx = beta <= 1 ? 0 : 1;

  double ECF3=0;
  for(int i=0; i<(int)cons.size(); i++){
        for(int j=i+1; j<(int)cons.size(); j++){
            for(int k=j+1; k<(int)cons.size(); k++){
                ECF3 += cons.at(i).Pt()*cons.at(j).Pt()*cons.at(k).Pt() * pow( cons.at(i).DeltaR(cons.at(j)) * cons.at(i).DeltaR(cons.at(k)) * cons.at(j).DeltaR(cons.at(k)), beta );
            }
        }
    }

    double ECF1=0;
    for(int i=0; i<(int)cons.size(); i++){
        ECF1 += cons.at(i).Pt();
    }

    double ECF2=0;
    for(int i=0; i<(int)cons.size(); i++){
        for(int j=i+1; j<(int)cons.size(); j++){
            ECF2 += cons.at(i).Pt()*cons.at(j).Pt() * pow( cons.at(i).DeltaR(cons.at(j)), beta );
        }
    }

    if( ECF3==0 || ECF2==0 || ECF1==0 ){
         cout<<"Error in eecorr calculation: something is 0"<<endl;

        return;
    }
    else{
      // TODO: dirty, should have the vector size sorted earlier
      if ((*var_ECF1_vec[jettype]).size() == 1)
	{
	  (*var_ECF1_vec[jettype])[beta_idx] = ECF1;
	  (*var_ECF2_vec[jettype])[beta_idx] = ECF2;
	  (*var_ECF3_vec[jettype])[beta_idx] = ECF3;
	}
      else
	{
	  (*var_ECF1_vec[jettype]).push_back(ECF1);
	  (*var_ECF2_vec[jettype]).push_back(ECF2);
	  (*var_ECF3_vec[jettype]).push_back(ECF3);
	}
        return;
    }

}//calculateECF

/*
 * In some samples the EEC variables are available already.  So little calculation needs to happen, but they do need to be set still.  This is a method that does this.
 * http://arxiv.org/abs/1411.0665  See equations 2 and 3.
 *
 * @params jettype Gives the index of the jet - truth/ topo or groomed
 * @params jetidx Is the index in the jet collection
 */
void setEEC(int jettype, int jetidx)
{
  // this is not obvious.  If we calculate EEC then we set the ecf variables explicitly, with index 0 for
  // beta 1 and index 1 for beta 2.
  // However, if we read ECF from the ntuple then the index should relate to the jet's index in the collection.
  int beta_1_idx = calcEEC ? 0 : jetidx;
  int beta_2_idx = calcEEC ? 1 : jetidx;
  // get ecf2 and 3 from the input ntuple/ calculated variables
  float ecf1_1 = (*var_ECF1_vec[jettype])[beta_1_idx];
  float ecf2_1 = (*var_ECF2_vec[jettype])[beta_1_idx];
  float ecf3_1 = (*var_ECF3_vec[jettype])[beta_1_idx];
  float e2_1 = ecf2_1/(ecf1_1*ecf1_1);
  float e3_1 = ecf3_1/(pow(ecf1_1,3));
  // calculate the EEC values using the formulae given in the paper linked above
  var_EEC_C2_1[jettype] = e3_1/pow(e2_1,2);
  var_EEC_D2_1[jettype] = e3_1/pow(e2_1,3);


  // get ecf2 and 3 from the input ntuple
  if (!xAOD) // the beta=2 ones do not exist in the xaods and there is no cluster info
    {
      float ecf1_2 = (*var_ECF1_vec[jettype])[beta_2_idx];
      float ecf2_2 = (*var_ECF2_vec[jettype])[beta_2_idx];
      float ecf3_2 = (*var_ECF3_vec[jettype])[beta_2_idx];
      float e2_2 = ecf2_2/(ecf1_2*ecf1_2);
      float e3_2 = ecf3_2/(pow(ecf1_2,3));
      var_EEC_C2_2[jettype] = e3_2/pow(e2_2,2);
      var_EEC_D2_2[jettype] = e3_2/pow(e2_2,3);
    }
  else if (beta2available) // some xAODs now have the beta 2 info
    {
      float ecf1_2 = (*var_ECF1_2_vec[jettype])[beta_2_idx];
      float ecf2_2 = (*var_ECF2_2_vec[jettype])[beta_2_idx];
      float ecf3_2 = (*var_ECF3_2_vec[jettype])[beta_2_idx];
      float e2_2 = ecf2_2/(ecf1_2*ecf1_2);
      float e3_2 = ecf3_2/(pow(ecf1_2,3));
      var_EEC_C2_2[jettype] = e3_2/pow(e2_2,2);
      var_EEC_D2_2[jettype] = e3_2/pow(e2_2,3);
    }
  else
    {
      var_EEC_C2_2[jettype] = -999;
      var_EEC_D2_2[jettype] = -999;
    }

  
}//setEEC


/*
 * Create TLVs containing the clusters for a given jet.
 *
 * @param jettype The type of jet - truth, topo or groomed.
 * @param jetidx The index of the jet in the truth/ topo or groomed collection.
 * @param cluster A vector<TLV> that will have the clusters added to it.
 *
 * @return A bool indicating if it was a success or not.  If the index ==-1 (which happens for truth jets) it will return false, for example.
 */
bool createClusters(int jettype, int jetidx, std::vector<TLorentzVector> & cluster)
{
  // number of constituents
  int size = (*var_constit_index[jettype])[jetidx].size();
  int cl_size = (*var_cl_pt_vec).size();
  cluster.reserve(size);

  // create TLV that gets reused each iteration
  TLorentzVector constit(0., 0., 0., 0.);

  //for (iterator::vector i = 0; i < size; i ++)
  for (vector<int>::iterator it = (*var_constit_index[jettype])[jetidx].begin(); it != (*var_constit_index[jettype])[jetidx].end(); it ++)
  //for (vector<int>::iterator it = tmpvec.begin(); it != tmpvec.end(); it ++)
    {
      // get the index in the cl_lc collection
      int idx = (*it);

      // if index is out of bounds there is no cluster info for this jet
      if (idx < 0 || idx >= cl_size)
	{
	  if (DEBUG)
	    std::cout << "warning: the index for the cluster is out of bounds" << std::endl;
	  return false;
	}

      if (DEBUG == 1)
	{
	  std::cout << "jet type: " << jettype << std::endl;
	  std::cout << "cluster idx: " << idx << std::endl;
	  std::cout << "cluster pt: " << (*var_cl_pt_vec)[idx] << std::endl;
	  std::cout << "cluster eta: " << (*var_cl_eta_vec)[idx] << std::endl;
	  std::cout << "cluster phi: " << (*var_cl_phi_vec)[idx] << std::endl;
	}
      // set values
      constit.SetPtEtaPhiM((*var_cl_pt_vec)[idx], (*var_cl_eta_vec)[idx], (*var_cl_phi_vec)[idx], 0.);
      // add to vector
      cluster.push_back(constit);
    }
  return true;
} // createClusters


/*
 * Calculate FoxWolfram20.
 * Formula used:
 * H_l = Sum_ij{frac{pi*pj}{E^2}*P_l(cos\theta_{ij}}, where P_l are legendre polynomials
 * Using the ratio between the 2nd and 0th moments. 
 *
 * @param clusters vector<TLV> of constituents of jet being analysed.
 * @return The calculated FoxWolfram20 value
 */
double calculateFoxWolfram20(vector<TLorentzVector> & clusters)
{
  // don't forget to assign the cl_ variables
  // The cluster variables should be for a single jet...
  float esum = 0;
  vector<double> FoxWolframMoments (5,0.0);
  // loop through i and j
  for (int i = 0; i < clusters.size(); i++)
    {
      //TLorentzVector cli;
      //cli.SetPtEtaPhiM((*var_cl_pt_vec)[i], (*var_cl_eta_vec)[i], (*var_cl_phi_vec)[i], 0);
      double pti = clusters[i].Pt();
      for (int j = i+1; j < clusters.size(); j++)
	{
	  //TLorentzVector clj;
	  //clj.SetPtEtaPhiM((*var_cl_pt_vec)[j], (*var_cl_eta_vec)[j], (*var_cl_phi_vec)[j], 0);
	  double ptj = clusters[j].Pt();
	  double costheta12 = TMath::Cos(clusters[j].Angle(clusters[i].Vect()));

	  // calculate the legendre polynomials
	  double p0 = 1.;
	  double p1 = costheta12;
	  double p2 = 0.5*(3.*costheta12*costheta12 - 1.);
	  double p3 = 0.5*(5.*costheta12*costheta12*costheta12 - 3.*costheta12);
	  double p4 = 0.125*(35.*costheta12*costheta12*costheta12*costheta12 - 30.*costheta12*costheta12 + 3.);
	  
	  // add to the moments array
	  FoxWolframMoments[0] += pti*ptj*p0;
	  FoxWolframMoments[1] += pti*ptj*p1;
	  FoxWolframMoments[2] += pti*ptj*p2;
	  FoxWolframMoments[3] += pti*ptj*p3;
	  FoxWolframMoments[4] += pti*ptj*p4;
			     
	} // loop over j
      esum+=clusters[i].E();

    } // loop i

  // divide by E^2
  for (int x = 0 ; x < 5; x++)
    {
      FoxWolframMoments[x]/= esum*esum;

    }
  //ratio of 2nd and 0th moments
  return double(FoxWolframMoments[2]/FoxWolframMoments[0]);

} // calculateFoxWolfram20

/*
 * Calculate the soft drop tag for a given jet.  This will be at a certain eff and fake rate. Definition found here: http://arxiv.org/pdf/1402.2657.pdf. Soft drop: highest pT, parameters: Beta=-1, zcut=0.04,0.06 and 0.08.
 * Code based on code found here: https://svnweb.cern.ch/trac/atlasoff/browser/PhysicsAnalysis/JetPhys/JetSubStructure/tags/JetSubStructure-00-00-21-10/src/SubstructureMethods/TaggingMethods.h
 *
 * @param cluster A vector<TLV> containing constituents of jet.
 *
 * @return Return the pt fraction for the highest softdrop with given eff and fake rate.
 */
int calculateSoftDropTag(vector<TLorentzVector> & cluster)
{
  // recluster jet with CA - don't be confused by the 1.5 radius- if a jet with R=0.4 is read in, then only the constituents from the jet will be used
  // in the rebuilt jet
  double radius=1.0;   
  fastjet::JetDefinition jd(fastjet::cambridge_algorithm, radius);
  boost::shared_ptr<const fastjet::ClusterSequence> cs(new const fastjet::ClusterSequence(cluster, jd));
  std::vector<fastjet::PseudoJet> jetVect = fastjet::sorted_by_pt(cs->inclusive_jets(0.));
  
  double beta = -1.0;
  int softdrop_level = 0;
  int highest_softdrop_level = 0;
  
  // set up pseudojets for parents and current jet
  fastjet::PseudoJet parent1, parent2;
  fastjet::PseudoJet currjet(jetVect[0]);

  float highest_pt_fraction = 0;
  // loop until current jet has no more parents.
  while (cs->has_parents(currjet, parent1, parent2))
    {
      // looking for leading parent
      if(parent1.perp2() < parent2.perp())
	{
	  std::swap(parent1,parent2);
	}
      // pt fraction of smaller sub-jet 
      double ptfraction = parent2.perp()/(parent1.perp() + parent2.perp());
      double power = pow((sqrt(parent1.plain_distance(parent2)))/radius, beta);
      // check pt fraction of current jet of pt of parents.
      if (ptfraction > 0.04*power) softdrop_level = 1; //eff 55%, fake rake 5%
      if (ptfraction > 0.06*power) softdrop_level = 2; //eff 30%, fake rate 2%
      if (ptfraction > 0.08*power) softdrop_level = 3; //eff 20%, fake rate 1%
      // new highest softdrop level
      if (softdrop_level > highest_softdrop_level)
	{
	  highest_softdrop_level = softdrop_level;
	  highest_pt_fraction = ptfraction;
	}
      // new current jet
      currjet = parent1;
    }
  
  return highest_softdrop_level; // highest_pt_fraction;
} //calculateSoftDropTag


/*
 * This method is based on the code from here:  http://jets.physics.harvard.edu/Qjets/html/Welcome.html
 * Take cells from an event, cluster them into jets, and perform Q-jet pruning on the hardest one

 * The main parameters governing the behavior of the algorithm are 
 * zcut: this is the z-cut used in the pruning algorithm (typically 0.1 or 0.15)
 * dcut_fctr: this is used in determining the angular dcut used in pruning. 
 as in the pruning algorithm, dcut = dcut_fctr * 2 m / pT (dcut_fctr is typically 0.5)
 * exp_min and exp_max determine the form of d_ij used in clustering the jet we have d_ij = min(pTi,pTj)^exp_min * max(pTi,pTj)^exp_max * R_ij^2, so for kT clustering, (exp_min,exp_max) = (2,0), while for C/A clustering they are (0,0)
 * rigidity: this determines how close the algorithm is to "classical" C/A or kT. This is denoted as \alpha in the Qjets paper
 *
 * @param clusters a vector<TLV> containing the input clusters for the jet
 * @return The q-jet volatility
*/
double calculateQJetsVol_v2(vector<TLorentzVector> & clusters)
{
  // sort input clusters by pt
  std::cout << "calc qjets" << std::endl;
  sort(clusters.begin(), clusters.end(), [] (TLorentzVector a, TLorentzVector b){
	  return a.Pt() > b.Pt();
	});//ComparePt);

  // create vector<PseudoJets> called input_particles, using input clusters for input
  vector<fastjet::PseudoJet> input_particles;
  for (std::vector<TLorentzVector>::iterator it = clusters.begin(); it != clusters.end(); it++)
    {
      input_particles.push_back(fastjet::PseudoJet((*it).Px(), (*it).Py(), (*it).Pz(), (*it).E()));
    }

  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm,radius,fastjet::E_scheme, fastjet::Best);
  fastjet::ClusterSequence clust_seq(input_particles, jet_def);  
  vector<fastjet::PseudoJet> inclusive_jets = clust_seq.inclusive_jets();
  
  // set up the parameters for the reclustering.
  // dcut = m/pt
  double dcut_fctr = sorted_by_pt(inclusive_jets)[0].m()/sorted_by_pt(inclusive_jets)[0].perp();
  double zcut(0.1), exp_min(0.), exp_max(0.), rigidity(0.1);
  double truncation_fctr(0.0);
 
  // initialise the plugin
  QjetsPlugin qjet_plugin(zcut, dcut_fctr, exp_min, exp_max, rigidity, truncation_fctr);
  fastjet::JetDefinition qjet_def(&qjet_plugin);

  vector<fastjet::PseudoJet> constits = clust_seq.constituents(sorted_by_pt(inclusive_jets)[0]);

  vector<double> masses;

  for(int i = 0 ; i < nqjets ; i++){
    fastjet::ClusterSequence qjet_seq(constits, qjet_def);
    vector<fastjet::PseudoJet> inclusive_jets2 = sorted_by_pt(qjet_seq.inclusive_jets());
    masses.push_back(inclusive_jets2[0].m());
  }

  // calculate variance and mean of mass distribution
  double massmean = mean(masses);
  double variance = var(masses, massmean);
  return (double)sqrt(variance)/massmean;
} // calculateQJetsVol_v2

/*
 * Set the radius for the jet algorithm using the prefix of the algorithm.
 *
 * @param prefix The prefix of the jet algorithm, which contains the radius.
 */
void setRadius(std::string & prefix)
{

  // set up the regex pattern, finds numbers at the end of the string.
  const char * pattern = "\\d+";
  boost::regex re(pattern);
  boost::sregex_iterator it(prefix.begin(), prefix.end(), re);
  boost::sregex_iterator end;
  // get each number returned by regex and append to a string
  string number = "";
  for ( ; it!=end; it++)
    {
      number.append((*it).str());
    }
  // convert the string to a float
  radius = stof(number);
  // divide by 10 to get the correct value
  radius/=10;

} //setRadius

/*
 * Print the entries in a vector<TLV>.
 *
 * @param tlv The vector of TLVs to print.
 */
void printTLV(vector<TLorentzVector> & tlv)
{
  for (int i = 0; i < tlv.size(); i++)
    {
      std::cout << "Entry in tlv number: " << i << std::endl;
      std::cout << "Pt: " << tlv[i].Pt() << std::endl;
      std::cout << "M: " << tlv[i].M() << std::endl;
      std::cout << "Eta: " << tlv[i].Eta() << std::endl;
      std::cout << "Phi: " << tlv[i].Phi() << std::endl;
    }
  
} //printTLV

/*
 * Create a map with the names of all of the branches within a tree from a TObjArray.
 *
 * @param arr The input TObjArray from which the map will be created.
 * @return unordered_map<string,bool> with list of branches in array.  All are set to false
 */
std::unordered_map<std::string, bool> createBranchMap(TObjArray *& arr)
{
  std::unordered_map<std::string, bool> br_map;
  // create object to iterate through the TObjArray
  TIter it(arr);
  TBranch * tb;
  while (tb = (TBranch *) it.Next())
    {
      std::string name = tb->GetName();
      // set map value
      br_map[name] = false;
    }
  
  return br_map;

} // createBranchMap


/*
 * create a float vector and add a variable to it so it is not empty
 */
void floatvec(vector<float> *& tmp)
{
  //vector<float> * 
  tmp = new vector<float>();
  tmp->push_back(-999);
  //return tmp;
} //floatvec

/*
 * create an int vector and add a variable to it so it is not empty
 */
//vector<int> * intvec()
void intvec(vector<int> *& tmp)
{
  //vector<int> * 
  tmp = new vector<int>();
  tmp->push_back(-999);
  //return tmp;
} //intvec


/*
 * create a tlv vector and add a variable to it so it is not empty
 */
//vector<TLorentzVector> * tlvvec()
void tlvvec(vector<TLorentzVector> *& tmp)
{
  //vector<TLorentzVector> * 
  tmp = new vector<TLorentzVector>();
  TLorentzVector tlv (0.,0.,0.,0.);
  tmp->push_back(tlv);
  //return tmp;
} //tlvvec

/*
 * create a vector<vector<int> > and add a variable to it so it is not empty
 */
//vector< vector<int> > * vecintvec()
void vecintvec(vector< vector<int> > *& tmp)
{
  vector<int> vec;
  vec.push_back(-999);
  //vector<vector<int> > * 
    tmp = new vector<vector<int> > ();
  tmp->push_back(vec);
  //return tmp;
} //vecintvec


/*
 * Used to catch exceptions where an attempted map element access is done where the element doesn't exist
 */
void SignalHandlerMapAccess(int signal)
{

  std::cout << "Signal: " << signal << std::endl;
  throw "!Access Violation!";

}
