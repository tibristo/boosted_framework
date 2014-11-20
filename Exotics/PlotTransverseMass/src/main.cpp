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
//#include "LinkDef.h"
using namespace boost::algorithm;
using namespace std;
namespace po = boost::program_options;

bool ComparePt(TLorentzVector a, TLorentzVector b) { return a.Pt() > b.Pt(); }

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
      ("branches-file", po::value<string>(&branchesFile)->default_value(""),"Name of file containing branches, otherwise Alg_branches.txt is used. Becareful with this, because if the branches are not read in then when any of them are used later on a segfault will occur, so make sure there is one that it will use.")
      ("xAOD-jets", po::value<bool>(&xAODJets)->default_value(false),"Indicate if we are running over xAOD output and there is the word Jets appended to the algorithm name.")
      ("xAOD-emfrac", po::value<bool>(&xAODemfrac)->default_value(false),"Indicate if we are running over xAOD output and emfrac is not available.")
      ("hvtllqq-selection", po::value<bool>(&hvtllqq)->default_value(false),"Indicate if we are running the HVTllqq analysis selection.")
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
  // create histograms that will be filled later on
  createHistos(alg_in_orig);

  // run the code to fill the mass histograms
  if (massHistos)
    {
      runAlgorithm(inputTChain[sampleType::BACKGROUND], inputTChain[sampleType::SIGNAL], algorithms.AlgoNames[alg_in_orig], alg_in_orig, massHistos);
  
      //Normalise all histograms of interest
      if (scaleHistsFlag)
	{
	  scaleHists();
	} // end scaleHists if
 
      //Plot things
      if (makePlotsFlag)
	{
	  makePlots(alg_in_orig);
	}
      
      //1. Get the MPV
      if (getMPVFlag)
	{
	  getMPV();
	}
  
      //2. Get the mass window which gives 68% W mass efficiency 
      if (applyMassWindowFlag)
	{
	  for (int j=0; j<nPtBins; j++){
	    Qw(WidthMassWindow[j],TopEdgeMassWindow[j],Wprime_Lead_CA12_mass[1][j], 0.68);
	    BottomEdgeMassWindow[j]=TopEdgeMassWindow[j]-WidthMassWindow[j];    
	  }
	  
	  for (int j=0; j<nFineBins; j++){
	    Qw(WidthMassWindow_finePt[j],TopEdgeMassWindow_finePt[j],Wprime_finePtBin_mass[1][j], 0.68);
	    BottomEdgeMassWindow_finePt[j]=TopEdgeMassWindow_finePt[j]-WidthMassWindow_finePt[j];    
	  }
	  
	}
      
      //3. Check background fraction in the 68% mass window  
      if (checkBkgFrac)
	{
	  for (int j=0; j<nPtBins; j++){
	    QCDfrac[j]=qcd_Lead_CA12_mass[histType::GROOMEDJET][j]->Integral(qcd_Lead_CA12_mass[histType::GROOMEDJET][j]->FindBin(BottomEdgeMassWindow[j]),qcd_Lead_CA12_mass[histType::GROOMEDJET][j]->FindBin(TopEdgeMassWindow[j]));
	    
	  }
	
	  for (int j=0; j<nFineBins; j++){
	    QCDfrac_finePt[j]=qcd_finePtBin_mass[histType::GROOMEDJET][j]->Integral(qcd_finePtBin_mass[histType::GROOMEDJET][j]->FindBin(BottomEdgeMassWindow_finePt[j]),qcd_finePtBin_mass[histType::GROOMEDJET][j]->FindBin(TopEdgeMassWindow_finePt[j]));
	  }
	  
	}

      if (makePtPlotsFlag)
	makePtPlots(alg_in_orig);
    } // if massHists

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
    
    while(tempFrac<frac && imax != Nbins){
      //fraction in bin imax=0,1,2,...                                                                                       
      tempFrac+=histo->GetBinContent(imax)/integral;
      imax+=1;
    }                                          
    
    double width = histo->GetBinCenter(imax) - histo->GetBinCenter(i); 
    
    double top_edge = histo->GetBinCenter(imax);
    
    if(imax != Nbins and width<minWidth){ 
      minWidth = width;
      topEdge = top_edge;
    }
  }  
}

/*
 *
 */
void runAlgorithm(TChain *inputTree, TChain *inputTree1, TString groomAlgo, std::string & groomAlgoIndex, bool massHistos)//, int fileidx)
{
  if (!massHistos)
    {
      nAlgos++;
      return;
    }
  initializeVariables();

  getMassHistograms(inputTree->GetTree(), inputTree1->GetTree(), groomAlgo, groomAlgoIndex);
  nAlgos++;
  deleteVectors();
}

/*
 *
 */
void getMassHistograms(TTree *inputTree, TTree *inputTree1, TString groomAlgo, std::string & groomAlgoIndex){
  
  cout<<"getBranches() for " << groomAlgo <<endl;
  getBranches(inputTree, inputTree1, groomAlgo, groomAlgoIndex);

  /////________________________________________________________________________
  /////  START LOOP : QCD
  /////________________________________________________________________________
  
  int nentries0=(int)inputTree->GetEntries();
  cout<<"start " << groomAlgo << " loop"<<endl;
  
  for (int ientry=0;ientry<nentries0;ientry++) {
    inputTree->GetEntry(ientry);
    
    //count total events
    nEvt_0=nEvt_0+1;
    
    //variable for QCD event weight
    Float_t weight = (qcd_mc_event_weight);
    //cout << "weight from tree " << weight << endl;
    UInt_t RunNumber = (qcd_mc_channel_number);
    //cout << "QCD channel number " << RunNumber << endl;
    
    //add weight to the QCD event from here https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/BosonTaggingPaper2012#7th_July_2014_Mixture_with_small
    if (RunNumber == 147913) weight*=1.6664E+03*1.9139E-03; //JZ3W
    if (RunNumber == 147914) weight*=2.7646E+01*1.4296E-03; //JZ4W
    if (RunNumber == 147915) weight*=3.0317E-01*5.5040E-03; //JZ5W
    if (RunNumber == 147916) weight*=7.5078E-03*1.5252E-02; //JZ6W   
    if (RunNumber == 147917) weight*=1.3760E-03*7.6369E-02; //JZ7W
    //cout << "weight after re-weight " << weight << endl;
    

    //NOW I WILL RECLUSTER MY JETS AND FEED THE APPROPRIATE VARIABLES 
    //FOR THE RECLUSTERING OPTIONS

    if (algorithms.AlgoType[groomAlgoIndex].find("recluster") != std::string::npos ){

      vector<TLorentzVector> small_jets;
      for (int n=0; n<(*qcd_CA12_topo_pt).size(); n++){
	TLorentzVector tempJet;
	float emfractmp = 0.5;
	if (xAODemfrac)
	  emfractmp = (*qcd_CA12_topo_emfrac)[n];
	if (emfractmp<0.99 && (*qcd_CA12_topo_pt)[n]/1000>15.0){
	  tempJet.SetPtEtaPhiE((*qcd_CA12_topo_pt)[n], (*qcd_CA12_topo_eta)[n], (*qcd_CA12_topo_phi)[n], (*qcd_CA12_topo_E)[n]);
	  small_jets.push_back(tempJet);
	}
      }

      vector<TLorentzVector> reclus_jets;
      

      reclus_jets = Recluster(small_jets);
      
      sort(reclus_jets.begin(), reclus_jets.end(), ComparePt);


      qcd_CA12_groomed_pt= new vector<float>();
      qcd_CA12_groomed_eta= new vector<float>();
      qcd_CA12_groomed_phi= new vector<float>();
      qcd_CA12_groomed_mass= new vector<float>();
      qcd_CA12_groomed_E= new vector<float>();
      qcd_CA12_groomed_emfrac= new vector<float>();
      
      
      //for (int n=0; n<reclus_jets.size(); n++){
      qcd_CA12_groomed_pt->push_back(reclus_jets[0].Pt());
      qcd_CA12_groomed_eta->push_back(reclus_jets[0].Eta());
      qcd_CA12_groomed_phi->push_back(reclus_jets[0].Phi());
      qcd_CA12_groomed_mass->push_back(reclus_jets[0].M());
      qcd_CA12_groomed_E->push_back(reclus_jets[0].E());
      qcd_CA12_groomed_emfrac->push_back(0.5);
      
      //}

    }

    /////________________________________________________________________________
    ///// I AM NOT GOING TO MATCH TRUTH JETS TO RECO ANYMORE
    ////--------------------------------------------------------------------------
   
    
    //check which is my CA12 ungroomed reco jet with EMfrac<0.99 and |eta|<1.2
    //loop over all topo jets and get the leading with cuts
    
    if (algorithms.AlgoType[groomAlgoIndex].find("truthmatch") != std::string::npos){
      
      bool hasTopoJet=false;
      int chosenTopoJetIndex=-99;
      float emfractmp = 0.5;
      for (int i=0; i<(*qcd_CA12_topo_pt).size(); i++){
	if (xAODemfrac)
	  emfractmp = (*qcd_CA12_topo_emfrac)[i];
	if (!hasTopoJet && emfractmp<0.99 && fabs((*qcd_CA12_topo_eta)[i])<1.2) {
	  chosenTopoJetIndex=i;
	  hasTopoJet=true;
	  nEvt_1=nEvt_1+1;
	} 
      }
      //cout << "leading topo jet with good selection is jet #" << chosenTopoJetIndex << " with pT " << (*qcdSF_Y9_CA12_topo_pt)[chosenTopoJetIndex]/1000 << endl;
      
      //now match my truth jet with the chosen topo jet
      int chosenLeadTruthJetIndex=-99;
      if (hasTopoJet){
	for (int i=0; i<(*qcd_CA12_truth_pt).size(); i++){
	  if (chosenLeadTruthJetIndex<0 && DeltaR((*qcd_CA12_topo_eta)[chosenTopoJetIndex],(*qcd_CA12_topo_phi)[chosenTopoJetIndex],(*qcd_CA12_truth_eta)[i],(*qcd_CA12_truth_phi)[i])<0.9){ 
	    chosenLeadTruthJetIndex=i;
	    nEvt_2=nEvt_2+1;
	  }	  
	}	
      }
      //qcd_Lead_CA12_pt[1]->Fill((*qcdSF_Y9_CA12_truth_pt)[0]/1000, weight);    
      
      //now I fill histograms for truth jet only
      
      if (chosenLeadTruthJetIndex>=0){
	qcd_Lead_CA12_pt[histType::TRUTHJET]->Fill((*qcd_CA12_truth_pt)[chosenLeadTruthJetIndex]/1000, weight);
	qcd_Lead_CA12_mass[histType::TRUTHJET][0]->Fill((*qcd_CA12_truth_mass)[chosenLeadTruthJetIndex]/1000, weight);
	
	if (((*qcd_CA12_truth_pt)[0]/1000)>1000.0){
	  qcd_Lead_CA12_mass[histType::TRUTHJET][1]->Fill((*qcd_CA12_truth_mass)[chosenLeadTruthJetIndex]/1000, weight);
	}
	
	if (((*qcd_CA12_truth_pt)[0]/1000)>500.0 && ((*qcd_CA12_truth_pt)[0]/1000)<1000.0){
	  qcd_Lead_CA12_mass[histType::TRUTHJET][2]->Fill((*qcd_CA12_truth_mass)[chosenLeadTruthJetIndex]/1000, weight);
	}
	if (((*qcd_CA12_truth_pt)[0]/1000)>1000.0 && ((*qcd_CA12_truth_pt)[0]/1000)<1500.0){
	  qcd_Lead_CA12_mass[histType::TRUTHJET][3]->Fill((*qcd_CA12_truth_mass)[chosenLeadTruthJetIndex]/1000, weight);
	}
	if (((*qcd_CA12_truth_pt)[0]/1000)>1500.0 && ((*qcd_CA12_truth_pt)[0]/1000)<2000.0){
	  qcd_Lead_CA12_mass[histType::TRUTHJET][4]->Fill((*qcd_CA12_truth_mass)[chosenLeadTruthJetIndex]/1000, weight);
	}
	if (((*qcd_CA12_truth_pt)[0]/1000)>2000.0){
	  qcd_Lead_CA12_mass[histType::TRUTHJET][5]->Fill((*qcd_CA12_truth_mass)[chosenLeadTruthJetIndex]/1000, weight);
	}
      }
    }
    
    
    //NOW BACK TO GROOMING AND BENS SELECTION

    int chosenLeadTruthJetIndex=0;
        
    //Now I have which events to make my pt reweight with, and to match to, etc
    //Start plotting things with the other algos for QCD
    int chosenLeadGroomedIndex=-99;

    for (int i=0; i<(*qcd_CA12_groomed_pt).size(); i++){
      float emfractmp = 0.5;
      if (xAODemfrac)
	emfractmp = (*qcd_CA12_groomed_emfrac)[i];
      if (chosenLeadTruthJetIndex>=0 && chosenLeadGroomedIndex<0 && DeltaR((*qcd_CA12_truth_eta)[chosenLeadTruthJetIndex],(*qcd_CA12_truth_phi)[chosenLeadTruthJetIndex],(*qcd_CA12_groomed_eta)[i],(*qcd_CA12_groomed_phi)[i])<0.9 && emfractmp<0.99 && fabs((*qcd_CA12_groomed_eta)[i])<1.2){
	// && (*qcd_CA12_groomed_pt)[i]/1000>100.0 ){ 
	
	chosenLeadGroomedIndex=i;
	//nEvt_3=nEvt_3+1;
      }     
    }
    

    // Now I will still ask for both selections of the same grooming algorithm to be in the sample
    // To simplify the coding below, I will say the truth jet is always 0
    chosenLeadTruthJetIndex=0; // wasn't this calculated earlier?
    
    if (chosenLeadGroomedIndex>=0 && chosenLeadTruthJetIndex>=0) {
      // arg why so much dividing by 1000???  Do this once!
      nEvt_3=nEvt_3+1;

      qcd_Lead_CA12_pt[histType::GROOMEDJET]->Fill((*qcd_CA12_groomed_pt)[chosenLeadGroomedIndex]/1000, weight);
      qcd_Lead_CA12_mass[histType::GROOMEDJET][0]->Fill((*qcd_CA12_groomed_mass)[chosenLeadGroomedIndex]/1000, weight);
      
      if (((*qcd_CA12_truth_pt)[0]/1000)>1000.0){
	qcd_Lead_CA12_mass[histType::GROOMEDJET][1]->Fill((*qcd_CA12_groomed_mass)[chosenLeadGroomedIndex]/1000, weight);
      }
      
      if (((*qcd_CA12_truth_pt)[0]/1000)>500.0 && ((*qcd_CA12_truth_pt)[0]/1000)<1000.0){
	qcd_Lead_CA12_mass[histType::GROOMEDJET][2]->Fill((*qcd_CA12_groomed_mass)[chosenLeadGroomedIndex]/1000, weight);
      }
      if (((*qcd_CA12_truth_pt)[0]/1000)>1000.0 && ((*qcd_CA12_truth_pt)[0]/1000)<1500.0){
	qcd_Lead_CA12_mass[histType::GROOMEDJET][3]->Fill((*qcd_CA12_groomed_mass)[chosenLeadGroomedIndex]/1000, weight);
      }
      if (((*qcd_CA12_truth_pt)[0]/1000)>1500.0 && ((*qcd_CA12_truth_pt)[0]/1000)<2000.0){
	qcd_Lead_CA12_mass[histType::GROOMEDJET][4]->Fill((*qcd_CA12_groomed_mass)[chosenLeadGroomedIndex]/1000, weight);
      }
      if (((*qcd_CA12_truth_pt)[0]/1000)>2000.0){
	qcd_Lead_CA12_mass[histType::GROOMEDJET][5]->Fill((*qcd_CA12_groomed_mass)[chosenLeadGroomedIndex]/1000, weight);
      }
      
      //now for fine pT bins
      if (((*qcd_CA12_truth_pt)[0]/1000)>0.0 && ((*qcd_CA12_truth_pt)[0]/1000)<250.0 ){
	qcd_finePtBin_mass[histType::GROOMEDJET][0]->Fill((*qcd_CA12_groomed_mass)[chosenLeadGroomedIndex]/1000, weight);
      }

      for (int k=0; k<11; k++){
	if (k==10){
	  if (((*qcd_CA12_truth_pt)[0]/1000)>250.0*k){
	    qcd_finePtBin_mass[histType::GROOMEDJET][k]->Fill((*qcd_CA12_groomed_mass)[chosenLeadGroomedIndex]/1000, weight);
	  }
	}
	else	if (((*qcd_CA12_truth_pt)[0]/1000)>250.0*k && ((*qcd_CA12_truth_pt)[0]/1000)<250.0*(k+1)){
	  qcd_finePtBin_mass[histType::GROOMEDJET][k]->Fill((*qcd_CA12_groomed_mass)[chosenLeadGroomedIndex]/1000, weight);
	}
      }


      
    }
    
    //for pt reweighing purposes
    qcd_PtReweight->Fill((*qcd_CA12_truth_pt)[0]/1000, weight);
    
  } // end loop over QCD events
  //end loop over QCD events
  
  //cout << "QCD efficiencies for " << groomAlgo << ": " << nEvt_0 << " " << nEvt_1 << " " << nEvt_2  << " " << nEvt_3  << " " << nEvt_4 << endl;
  
  
  //NOW I NEED A LOOP OVER Wp SAMPLES TO GET THEIR WEIGHT
  // I want to put this into another method

  int nentriesW=(int)inputTree1->GetEntries();
  cout<<"start Wp weight loop"<<endl;
  
  for (int ientry=0;ientry<nentriesW;ientry++) {
    inputTree1->GetEntry(ientry);

    Wp_PtReweight->Fill((*Wp_CA12_truth_pt)[0]/1000);

  }


  //NOW I HAVE THE HISTOGRAMS TO DIVIDE AND GET THE WEIGHT HISTOGRAM

  for (int j=1; j<21; j++){
    //loop bins 1-20
    float weight=0.0;
    if (Wp_PtReweight->GetBinContent(j)>0){
      weight = qcd_PtReweight->GetBinContent(j)/Wp_PtReweight->GetBinContent(j);
      pTweights->SetBinContent(j,weight);
    }
    else pTweights->SetBinContent(j,0.0);
    //cout << weight << endl;
  }


  
  //START LOOP: over Wp samples
  int nentries1=(int)inputTree1->GetEntries();
  cout<<"start Wp " << groomAlgo << " loop"<<endl;
  
  for (int ientry=0;ientry<nentries1;ientry++) {
    inputTree1->GetEntry(ientry);
    
    nEvt1_0++;;
    

    double weight = 1.0;
    
    //Now I have to scale ever 175 GeV
    // this is the pT I need to use for the bin and rescale
    for (int j=0; j<20; j++){
      if ((*Wp_CA12_truth_pt)[0]/1000>j*175.0 && (*Wp_CA12_truth_pt)[0]/1000<=175.0*(j+1)){
	weight = pTweights->GetBinContent(j+1); // nothing is being scaled here?!  This doesn't do anything
      }
    }
    
    

    //NOW I WILL RECLUSTER MY JETS AND FEED THE APPROPRIATE VARIABLES 
    //FOR THE RECLUSTERING OPTIONS

    if (algorithms.AlgoType[groomAlgoIndex].find("recluster") != std::string::npos){

      vector<TLorentzVector> small_jets;
      for (int n=0; n<(*Wp_CA12_topo_pt).size(); n++){
	float emfractmp = 0.5;
	if (xAODemfrac)
	  emfractmp = (*Wp_CA12_topo_emfrac)[n];
	TLorentzVector tempJet;
	if (emfractmp<0.99){
	  tempJet.SetPtEtaPhiE((*Wp_CA12_topo_pt)[n], (*Wp_CA12_topo_eta)[n], (*Wp_CA12_topo_phi)[n], (*Wp_CA12_topo_E)[n]);
	  small_jets.push_back(tempJet);
	}
      }

      vector<TLorentzVector> reclus_jets;
      
      reclus_jets = Recluster(small_jets);
      
      sort(reclus_jets.begin(), reclus_jets.end(), ComparePt);

      Wp_CA12_groomed_pt= new vector<float>();
      Wp_CA12_groomed_eta= new vector<float>();
      Wp_CA12_groomed_phi= new vector<float>();
      Wp_CA12_groomed_mass= new vector<float>();
      Wp_CA12_groomed_E= new vector<float>();
      Wp_CA12_groomed_emfrac= new vector<float>();
      

      //for (int n=0; n<reclus_jets.size(); n++){
      Wp_CA12_groomed_pt->push_back(reclus_jets[0].Pt());
      Wp_CA12_groomed_eta->push_back(reclus_jets[0].Eta());
      Wp_CA12_groomed_phi->push_back(reclus_jets[0].Phi());
      Wp_CA12_groomed_mass->push_back(reclus_jets[0].M());
      Wp_CA12_groomed_E->push_back(reclus_jets[0].E());
      Wp_CA12_groomed_emfrac->push_back(0.5);
	
      // }
      
    }


    
    //now for ungroomed truth:
      
    if (algorithms.AlgoType[groomAlgoIndex].find("truthmatch") != std::string::npos){
      //check which is my CA12 ungroomed reco jet with EMfrac<0.99 and |eta|<1.2
      //loop over all topo jets and get the leading with cuts
      bool hasTopoJet=false;
      int chosenTopoJetIndex=-99;
      for (int i=0; i<(*Wp_CA12_topo_pt).size(); i++){
	float emfractmp = 0.5;
	if (xAODemfrac)
	  emfractmp = (*Wp_CA12_topo_emfrac)[i];
	if (!hasTopoJet && emfractmp<0.99 && fabs((*Wp_CA12_topo_eta)[i])<1.2) {
	  chosenTopoJetIndex=i;
	  hasTopoJet=true;
	  nEvt1_1=nEvt1_1+1;
	} 
      }
      //cout << "leading topo jet with good selection is jet #" << chosenTopoJetIndex << " with pT " << (*WpSF_CA12_topo_pt)[chosenTopoJetIndex]/1000 << endl;
      
      //now match my truth jet with the chosen topo jet
      int chosenLeadTruthJetIndex=-99;
      if (hasTopoJet){
	for (int i=0; i<(*Wp_CA12_truth_pt).size(); i++){
	  if (chosenLeadTruthJetIndex<0 && DeltaR((*Wp_CA12_topo_eta)[chosenTopoJetIndex],(*Wp_CA12_topo_phi)[chosenTopoJetIndex],(*Wp_CA12_truth_eta)[i],(*Wp_CA12_truth_phi)[i])<0.9){ 
	    chosenLeadTruthJetIndex=i;
	    nEvt1_2=nEvt1_2+1;
	  }	  
	}	
      }

      if (chosenLeadTruthJetIndex>=0){
	Wprime_Lead_CA12_pt[histType::GROOMEDJET]->Fill((*Wp_CA12_truth_pt)[chosenLeadTruthJetIndex]/1000, weight);
	Wprime_Lead_CA12_mass[histType::GROOMEDJET][0]->Fill((*Wp_CA12_truth_mass)[chosenLeadTruthJetIndex]/1000, weight);
	
	if (((*Wp_CA12_truth_pt)[0]/1000)>1000.0){
	  Wprime_Lead_CA12_mass[histType::GROOMEDJET][1]->Fill((*Wp_CA12_truth_mass)[chosenLeadTruthJetIndex]/1000, weight);
	}
	
	if (((*Wp_CA12_truth_pt)[0]/1000)>500.0 && ((*Wp_CA12_truth_pt)[0]/1000)<1000.0){
	  Wprime_Lead_CA12_mass[histType::GROOMEDJET][2]->Fill((*Wp_CA12_truth_mass)[chosenLeadTruthJetIndex]/1000, weight);
	}
	if (((*Wp_CA12_truth_pt)[0]/1000)>1000.0 && ((*Wp_CA12_truth_pt)[0]/1000)<1500.0){
	  Wprime_Lead_CA12_mass[histType::GROOMEDJET][3]->Fill((*Wp_CA12_truth_mass)[chosenLeadTruthJetIndex]/1000, weight);
	}
	if (((*Wp_CA12_truth_pt)[0]/1000)>1500.0 && ((*Wp_CA12_truth_pt)[0]/1000)<2000.0){
	  Wprime_Lead_CA12_mass[histType::GROOMEDJET][4]->Fill((*Wp_CA12_truth_mass)[chosenLeadTruthJetIndex]/1000, weight);
	}
	if (((*Wp_CA12_truth_pt)[0]/1000)>2000.0){
	  Wprime_Lead_CA12_mass[histType::GROOMEDJET][5]->Fill((*Wp_CA12_truth_mass)[chosenLeadTruthJetIndex]/1000, weight);
	}	
      }

    } //end of groomedAlgo=0 if

    
    

    int chosenLeadTruthJetIndex=0;
    
    //Now I have which events to make my pt reweight with, and to match to, etc
    //Start plotting things with the other algos for QCD
    int chosenLeadGroomedIndex=-99;
    for (int i=0; i<(*Wp_CA12_groomed_pt).size(); i++){
      float emfractmp = 0.5;
      if (xAODemfrac)
	emfractmp = (*Wp_CA12_groomed_emfrac)[i]; 
      if (chosenLeadTruthJetIndex>=0 && chosenLeadGroomedIndex<0 && DeltaR((*Wp_CA12_truth_eta)[chosenLeadTruthJetIndex],(*Wp_CA12_truth_phi)[chosenLeadTruthJetIndex],(*Wp_CA12_groomed_eta)[i],(*Wp_CA12_groomed_phi)[i])<0.9 && emfractmp<0.99 && fabs((*Wp_CA12_groomed_eta)[i])<1.2){
	//  && (*Wp_CA12_groomed_pt)[i]/1000>100.0 ){ 	
	chosenLeadGroomedIndex=i;
	//nEvt1_3=nEvt1_3+1;
      }     
    }


       
    if (chosenLeadGroomedIndex>=0 && chosenLeadTruthJetIndex>=0) {

      nEvt1_3=nEvt1_3+1;

      Wprime_Lead_CA12_pt[histType::GROOMEDJET]->Fill((*Wp_CA12_groomed_pt)[chosenLeadGroomedIndex]/1000, weight);
      Wprime_Lead_CA12_mass[histType::GROOMEDJET][0]->Fill((*Wp_CA12_groomed_mass)[chosenLeadGroomedIndex]/1000, weight);
      
      if (((*Wp_CA12_truth_pt)[0]/1000)>1000.0){
	Wprime_Lead_CA12_mass[histType::GROOMEDJET][1]->Fill((*Wp_CA12_groomed_mass)[chosenLeadGroomedIndex]/1000, weight);
      }
      
      if (((*Wp_CA12_truth_pt)[0]/1000)>500.0 && ((*Wp_CA12_truth_pt)[0]/1000)<1000.0){
	Wprime_Lead_CA12_mass[histType::GROOMEDJET][2]->Fill((*Wp_CA12_groomed_mass)[chosenLeadGroomedIndex]/1000, weight);
      }
      if (((*Wp_CA12_truth_pt)[0]/1000)>1000.0 && ((*Wp_CA12_truth_pt)[0]/1000)<1500.0){
	Wprime_Lead_CA12_mass[histType::GROOMEDJET][3]->Fill((*Wp_CA12_groomed_mass)[chosenLeadGroomedIndex]/1000, weight);
      }
      if (((*Wp_CA12_truth_pt)[0]/1000)>1500.0 && ((*Wp_CA12_truth_pt)[0]/1000)<2000.0){
	Wprime_Lead_CA12_mass[histType::GROOMEDJET][4]->Fill((*Wp_CA12_groomed_mass)[chosenLeadGroomedIndex]/1000, weight);
      }
      if (((*Wp_CA12_truth_pt)[0]/1000)>2000.0){
	//if ((*Wp_CA12_groomed_pt)[chosenLeadGroomedIndex]/1000<500.00)
	//cout << (*Wp_CA12_groomed_pt)[chosenLeadGroomedIndex]/1000 << endl;
	Wprime_Lead_CA12_mass[histType::GROOMEDJET][5]->Fill((*Wp_CA12_groomed_mass)[chosenLeadGroomedIndex]/1000, weight);
      }

      //fine pt bins      
      for (int k=0; k<12; k++){
	if (((*Wp_CA12_truth_pt)[0]/1000)>250.0*k && ((*Wp_CA12_truth_pt)[0]/1000)<250.0*(k+1)){
	  Wprime_finePtBin_mass[histType::GROOMEDJET][k]->Fill((*Wp_CA12_groomed_mass)[chosenLeadGroomedIndex]/1000, weight);
	}
      }
      
    }
    //for rescaling check purposes
    
    Wprime_Lead_CA12_scaled_pt->Fill((*Wp_CA12_truth_pt)[0]/1000, weight);
    /*
      Wprime_Lead_CA12_scaled_pt[8]->Fill((*Wp_CA12_truth_pt)[0]/1000, weight);
      Wprime_Lead_CA12_mass[8][0]->Fill((*Wp_CA12_truth_mass)[0]/1000, weight);
    */
  }
  //end loop over Wprime events
  
  cout << "Wprime efficiencies for " << groomAlgo << ": " << nEvt1_0 << " " << nEvt1_1 << " " << nEvt1_2  << " " << nEvt1_3  << " " << nEvt1_4 << endl;
  

}


/*
 *
 */
void createHistos(std::string & alg){

  //create my histograms
  // make an array with histograms for all sorts of algorithms used
  int M = 100;
  for (int i=0; i<2; i++){ // truth and groomed
    //std::cout << nAlgos << std::endl;
    //for (int ii=0; ii<nAlgos; ii++){
    //int i = algoMap[ii];
    //std::cout << i << std::endl;
    std::string alg_idx = i == 0 ? "CamKt12Truth" : alg;
    for (int j=0; j<nFineBins; j++){
      qcd_finePtBin_mass[i][j] = new TH1F(std::string("qcd_finePtBin_mass_"+algorithms.AlgoList[alg_idx]+finePtBins[j]).c_str(), std::string("qcd_finePtBin_mass_"+algorithms.AlgoList[alg_idx]+finePtBins[j]).c_str(), 150, 0.0, 300.0);
      //qcd_finePtBin_mass_curve[i][j] = new TH1D("qcd_finePtBin_mass_"+algorithms.AlgoList[alg_idx]+finePtBins[j]+"_ROC","qcd_finePtBin_mass_"+algorithms.AlgoList[alg_idx]+finePtBins[j]+"_ROC",M,0,1);
      Wprime_finePtBin_mass[i][j] = new TH1F(std::string("Wprime_finePtBin_mass_"+algorithms.AlgoList[alg_idx]+finePtBins[j]).c_str(), std::string("Wprime_finePtBin_mass_"+algorithms.AlgoList[alg_idx]+finePtBins[j]).c_str(), 150, 0.0, 300.0);

      qcd_finePtBin_mass[i][j]->Sumw2();
      Wprime_finePtBin_mass[i][j]->Sumw2();
    }


    for (int j=0; j<nPtBins; j++){
      qcd_Lead_CA12_mass[i][j] = new TH1F(std::string("qcd_Lead_CA12_mass_"+algorithms.AlgoList[alg_idx]+pTbins[j]).c_str(), std::string("qcd_Lead_CA12_mass_"+algorithms.AlgoList[alg_idx]+pTbins[j]).c_str(), 300, 0.0, 300.0);

      qcd_Lead_CA12_mass[i][j]->Sumw2();
      Wprime_Lead_CA12_mass[i][j] = new TH1F(std::string("Wprime_Lead_CA12_mass_"+algorithms.AlgoList[alg_idx]+pTbins[j]).c_str(), std::string("Wprime_Lead_CA12_mass_"+algorithms.AlgoList[alg_idx]+pTbins[j]).c_str(), 300, 0.0, 300.0);

      Wprime_Lead_CA12_mass[i][j]->Sumw2();
    }

    qcd_Lead_CA12_pt[i] = new TH1F(std::string("qcd_Lead_CA12_pt_"+algorithms.AlgoList[alg_idx]).c_str(), std::string("qcd_Lead_CA12_pt_"+algorithms.AlgoList[alg_idx]).c_str(), 20, 0.0, 3500.0);
    Wprime_Lead_CA12_pt[i] = new TH1F(std::string("Wprime_Lead_CA12_pt_"+algorithms.AlgoList[alg_idx]).c_str(), std::string("Wprime_Lead_CA12_pt_"+algorithms.AlgoList[alg_idx]).c_str(), 20, 0.0, 3500.0);
    
    

      
  }
  
  qcd_Lead_CA12_mass[2][0] = new TH1F(std::string("qcd_Lead_CA12_mass_"+algorithms.AlgoList["LeadTruthJet"]).c_str(), std::string("qcd_Lead_CA12_mass_"+algorithms.AlgoList["LeadTruthJet"]).c_str(), 100, 0.0, 300.0);
  qcd_Lead_CA12_mass[2][0]->Sumw2();
  Wprime_Lead_CA12_mass[2][0] = new TH1F(std::string("Wprime_Lead_CA12_mass_"+algorithms.AlgoList["LeadTruthJet"]).c_str(), std::string("Wprime_Lead_CA12_mass_"+algorithms.AlgoList["LeadTruthJet"]).c_str(), 100, 0.0, 300.0);
  Wprime_Lead_CA12_mass[2][0]->Sumw2();
  
  qcd_Lead_CA12_pt[2] = new TH1F(std::string("qcd_Lead_CA12_pt_"+algorithms.AlgoList["LeadTruthJet"]).c_str(), std::string("qcd_Lead_CA12_pt_"+algorithms.AlgoList["LeadTruthJet"]).c_str(), 350, 0.0, 3500.0);
  Wprime_Lead_CA12_pt[2] = new TH1F(std::string("Wprime_Lead_CA12_pt_"+algorithms.AlgoList["LeadTruthJet"]).c_str(), std::string("Wprime_Lead_CA12_pt_"+algorithms.AlgoList["LeadTruthJet"]).c_str(), 350, 0.0, 3500.0);
  //Wprime_Lead_CA12_scaled_pt[2] = new TH1F(std::string("Wprime_Lead_CA12_scaled_pt_"+algorithms.AlgoList["LeadTruthJet"]).c_str(), std::string("Wprime_Lead_CA12_scaled_pt_"+algorithms.AlgoList["LeadTruthJet"]).c_str(), 350, 0.0, 3500.0);
  
  Wprime_Lead_CA12_scaled_pt = new TH1F(std::string("Wprime_Lead_CA12_scaled_pt_"+algorithms.AlgoList[alg]).c_str(), std::string("Wprime_Lead_CA12_scaled_pt_"+algorithms.AlgoList[alg]).c_str(), 20, 0.0, 3500.0);

  qcd_PtReweight = new TH1F(std::string("qcd_ptreweight_"+algorithms.AlgoNames[alg]).c_str(), std::string("qcd_ptreweight_"+algorithms.AlgoNames[alg]).c_str(), 20, 0.0, 3500.0);
  Wp_PtReweight = new TH1F(std::string("Wp_ptreweight_"+algorithms.AlgoNames[alg]).c_str(), std::string("Wp_ptreweight_"+algorithms.AlgoNames[alg]).c_str(), 20, 0.0, 3500.0);
  pTweights = new TH1F(std::string("pT_weights_"+algorithms.AlgoNames[alg]).c_str(), std::string("pt weights_"+algorithms.AlgoNames[alg]).c_str(), 20, 0.0, 3500.0);

} // createHistos

/*
 * Initialise all of the qcd and WPrime variables that are used for reading in from the input files.
 */
void initializeVariables(){

  qcd_CA12_truth_pt = 0;
  qcd_CA12_truth_eta = 0;
  qcd_CA12_truth_phi = 0;
  qcd_CA12_truth_mass = 0;
  qcd_CA12_truth_E = 0;
  
  qcd_CA12_topo_pt = 0;
  qcd_CA12_topo_eta = 0;
  qcd_CA12_topo_phi = 0;
  qcd_CA12_topo_mass = 0;
  qcd_CA12_topo_E = 0;
  qcd_CA12_topo_emfrac = 0;
  
  qcd_CA12_groomed_pt = 0;
  qcd_CA12_groomed_eta = 0;
  qcd_CA12_groomed_phi = 0;
  qcd_CA12_groomed_mass = 0;
  qcd_CA12_groomed_E = 0;
  qcd_CA12_groomed_emfrac = 0;
  
  qcd_mc_channel_number = 0; 
  qcd_mc_event_weight = 0; 

  //Wprime split filtering both cuts

  Wp_CA12_truth_pt = 0;
  Wp_CA12_truth_eta = 0;
  Wp_CA12_truth_phi = 0;
  Wp_CA12_truth_mass = 0;
  Wp_CA12_truth_E = 0;
  Wp_mc_channel_number = 0; 
  //Wp_mc_event_weight = 0; 

  Wp_CA12_topo_pt = 0;
  Wp_CA12_topo_eta = 0;
  Wp_CA12_topo_phi = 0;
  Wp_CA12_topo_mass = 0;
  Wp_CA12_topo_E = 0;
  Wp_CA12_topo_emfrac = 0;
  
  Wp_CA12_groomed_pt = 0;
  Wp_CA12_groomed_eta = 0;
  Wp_CA12_groomed_phi = 0;
  Wp_CA12_groomed_mass = 0;
  Wp_CA12_groomed_E = 0;
  Wp_CA12_groomed_emfrac = 0;
 
  nEvt_0 = 0;
  nEvt_1 = 0;
  nEvt_2 = 0;
  nEvt_3 = 0;
  nEvt_4 = 0;
  nEvt1_0 = 0;
  nEvt1_1 = 0;
  nEvt1_2 = 0;
  nEvt1_3 = 0;
  nEvt1_4 = 0;

   
} //initialiseVariables


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
 * Delete the wprime and qcd vectors used for reading in from the file.
 */
void deleteVectors(){


  delete qcd_CA12_groomed_pt;
  delete qcd_CA12_groomed_eta;
  delete qcd_CA12_groomed_phi;
  delete qcd_CA12_groomed_mass;
  delete qcd_CA12_groomed_E;
  delete qcd_CA12_groomed_emfrac;
      
  delete Wp_CA12_groomed_pt;
  delete Wp_CA12_groomed_eta;
  delete Wp_CA12_groomed_phi;
  delete Wp_CA12_groomed_mass;
  delete Wp_CA12_groomed_E;
  delete Wp_CA12_groomed_emfrac;


}

/*
 * Set up the input branches for the wprime vs qcd comparison for truth and groomed algorithm.
 *
 * @param inputTree Pointer to input background TTree.
 * @param inputTree1 Pointer to input signal TTree.
 * @param groomAlgo Abbreviated name of the grooming algorithm.
 * @param groomAlgoIndex The full name of the algorithm that is used to access the algorithms maps.
 */
void getBranches(TTree *inputTree, TTree *inputTree1, TString groomAlgo, std::string & groomAlgoIndex){

  //if not reclustering
  if (algorithms.AlgoType[groomAlgoIndex].find("recluster") == std::string::npos){
    
    inputTree->SetBranchAddress("mc_channel_number", &qcd_mc_channel_number); // segfaults for groomalgo = <incomplete type>, groomalgoindex = 3, TopoSplitFilteredMu100SmallR30YCut4
    inputTree->SetBranchAddress("mc_event_weight", &qcd_mc_event_weight);
    inputTree->SetBranchAddress(std::string("jet_CamKt12Truth_pt").c_str(), &qcd_CA12_truth_pt);
    inputTree->SetBranchAddress(std::string("jet_CamKt12Truth_eta").c_str(), &qcd_CA12_truth_eta);
    inputTree->SetBranchAddress(std::string("jet_CamKt12Truth_phi").c_str(), &qcd_CA12_truth_phi);
    inputTree->SetBranchAddress(std::string("jet_CamKt12Truth_m").c_str(), &qcd_CA12_truth_mass);
    inputTree->SetBranchAddress(std::string("jet_CamKt12Truth_E").c_str(), &qcd_CA12_truth_E);
    
    inputTree->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LCTopo_pt").c_str(), &qcd_CA12_topo_pt);
    inputTree->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LCTopo_eta").c_str(), &qcd_CA12_topo_eta);
    inputTree->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LCTopo_phi").c_str(), &qcd_CA12_topo_phi);
    inputTree->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LCTopo_m").c_str(), &qcd_CA12_topo_mass);
    inputTree->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LCTopo_E").c_str(), &qcd_CA12_topo_E);
    if (inputTree->GetListOfBranches()->FindObject(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LCTopo_emfrac").c_str()) )
      {
	inputTree->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LCTopo_emfrac").c_str(), &qcd_CA12_topo_emfrac);
      }
    inputTree->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LC"+groomAlgo+"_pt").c_str(), &qcd_CA12_groomed_pt);
    inputTree->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LC"+groomAlgo+"_eta").c_str(), &qcd_CA12_groomed_eta);
    inputTree->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LC"+groomAlgo+"_phi").c_str(), &qcd_CA12_groomed_phi);
    inputTree->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LC"+groomAlgo+"_m").c_str(), &qcd_CA12_groomed_mass);
    inputTree->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LC"+groomAlgo+"_E").c_str(), &qcd_CA12_groomed_E);
    if (inputTree->GetListOfBranches()->FindObject(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LC" + groomAlgo+"_emfrac").c_str()) )
      {
	inputTree->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LC"+groomAlgo+"_emfrac").c_str(), &qcd_CA12_groomed_emfrac);
      }
    
    inputTree1->SetBranchAddress("mc_channel_number", &Wp_mc_channel_number);
    inputTree1->SetBranchAddress("mc_event_weight", &Wp_mc_event_weight);
    inputTree1->SetBranchAddress(std::string("jet_CamKt12Truth_pt").c_str(), &Wp_CA12_truth_pt);
    inputTree1->SetBranchAddress(std::string("jet_CamKt12Truth_eta").c_str(), &Wp_CA12_truth_eta);
    inputTree1->SetBranchAddress(std::string("jet_CamKt12Truth_phi").c_str(), &Wp_CA12_truth_phi);
    inputTree1->SetBranchAddress(std::string("jet_CamKt12Truth_m").c_str(), &Wp_CA12_truth_mass);
    inputTree1->SetBranchAddress(std::string("jet_CamKt12Truth_E").c_str(), &Wp_CA12_truth_E);
    
    inputTree1->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LCTopo_pt").c_str(), &Wp_CA12_topo_pt);
    inputTree1->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LCTopo_eta").c_str(), &Wp_CA12_topo_eta);
    inputTree1->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LCTopo_phi").c_str(), &Wp_CA12_topo_phi);
    inputTree1->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LCTopo_m").c_str(), &Wp_CA12_topo_mass);
    inputTree1->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LCTopo_E").c_str(), &Wp_CA12_topo_E);
    if (inputTree1->GetListOfBranches()->FindObject(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LCTopo_emfrac").c_str()) )
      {
	inputTree1->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LCTopo_emfrac").c_str(), &Wp_CA12_topo_emfrac);
      }
    
    inputTree1->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LC"+groomAlgo+"_pt").c_str(), &Wp_CA12_groomed_pt);
    inputTree1->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LC"+groomAlgo+"_eta").c_str(), &Wp_CA12_groomed_eta);
    inputTree1->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LC"+groomAlgo+"_phi").c_str(), &Wp_CA12_groomed_phi);
    inputTree1->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LC"+groomAlgo+"_m").c_str(), &Wp_CA12_groomed_mass);
    inputTree1->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LC"+groomAlgo+"_E").c_str(), &Wp_CA12_groomed_E);
    if (inputTree1->GetListOfBranches()->FindObject(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LC"+groomAlgo+"_emfrac").c_str()) )
      {
	inputTree1->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LC"+groomAlgo+"_emfrac").c_str(), &Wp_CA12_groomed_emfrac);
      }


  }
 



  //if reclustering algorithm
  else if (algorithms.AlgoType[groomAlgoIndex].find("recluster") != std::string::npos){

    inputTree->SetBranchAddress("mc_channel_number", &qcd_mc_channel_number);
    inputTree->SetBranchAddress("mc_event_weight", &qcd_mc_event_weight);
    inputTree->SetBranchAddress(std::string("jet_CamKt12Truth_pt").c_str(), &qcd_CA12_truth_pt);
    inputTree->SetBranchAddress(std::string("jet_CamKt12Truth_eta").c_str(), &qcd_CA12_truth_eta);
    inputTree->SetBranchAddress(std::string("jet_CamKt12Truth_phi").c_str(), &qcd_CA12_truth_phi);
    inputTree->SetBranchAddress(std::string("jet_CamKt12Truth_m").c_str(), &qcd_CA12_truth_mass);
    inputTree->SetBranchAddress(std::string("jet_CamKt12Truth_E").c_str(), &qcd_CA12_truth_E);
    
    //AntiKt2LCTopo

    inputTree->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+""+groomAlgo+"_pt").c_str(), &qcd_CA12_groomed_pt);
    inputTree->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+""+groomAlgo+"_eta").c_str(), &qcd_CA12_groomed_eta);
    inputTree->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+""+groomAlgo+"_phi").c_str(), &qcd_CA12_groomed_phi);
    inputTree->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+""+groomAlgo+"_m").c_str(), &qcd_CA12_groomed_mass);
    inputTree->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+""+groomAlgo+"_E").c_str(), &qcd_CA12_groomed_E);
    if (inputTree->GetListOfBranches()->FindObject(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+groomAlgo+"_emfrac").c_str()) )
      {
	inputTree->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+""+groomAlgo+"_emfrac").c_str(), &qcd_CA12_groomed_emfrac);
      }
    
    
    
    inputTree1->SetBranchAddress("mc_channel_number", &Wp_mc_channel_number);
    //inputTree1->SetBranchAddress("mc_event_weight", &Wp_mc_event_weight);
    inputTree1->SetBranchAddress(std::string("jet_CamKt12Truth_pt").c_str(), &Wp_CA12_truth_pt);
    inputTree1->SetBranchAddress(std::string("jet_CamKt12Truth_eta").c_str(), &Wp_CA12_truth_eta);
    inputTree1->SetBranchAddress(std::string("jet_CamKt12Truth_phi").c_str(), &Wp_CA12_truth_phi);
    inputTree1->SetBranchAddress(std::string("jet_CamKt12Truth_m").c_str(), &Wp_CA12_truth_mass);
    inputTree1->SetBranchAddress(std::string("jet_CamKt12Truth_E").c_str(), &Wp_CA12_truth_E);
      
    inputTree1->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LCTopo_pt").c_str(), &Wp_CA12_topo_pt);
    inputTree1->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LCTopo_eta").c_str(), &Wp_CA12_topo_eta);
    inputTree1->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LCTopo_phi").c_str(), &Wp_CA12_topo_phi);
    inputTree1->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LCTopo_m").c_str(), &Wp_CA12_topo_mass);
    inputTree1->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LCTopo_E").c_str(), &Wp_CA12_topo_E);
    if (inputTree1->GetListOfBranches()->FindObject(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LCTopo_emfrac").c_str()) )
      {
	inputTree1->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LCTopo_emfrac").c_str(), &Wp_CA12_topo_emfrac);
      }
    
    inputTree1->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+""+groomAlgo+"_pt").c_str(), &Wp_CA12_groomed_pt);
    inputTree1->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+""+groomAlgo+"_eta").c_str(), &Wp_CA12_groomed_eta);
    inputTree1->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+""+groomAlgo+"_phi").c_str(), &Wp_CA12_groomed_phi);
    inputTree1->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+""+groomAlgo+"_m").c_str(), &Wp_CA12_groomed_mass);
    inputTree1->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+""+groomAlgo+"_E").c_str(), &Wp_CA12_groomed_E);
    if (inputTree1->GetListOfBranches()->FindObject(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+groomAlgo+"_emfrac").c_str()) )
      {
	inputTree1->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+""+groomAlgo+"_emfrac").c_str(), &Wp_CA12_groomed_emfrac);
      }
    

  }
} //getBranches()

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
 * Create mass plots.
 * 
 * @param algo The name of the grooming algorithm being run.
 */
void makePlots(std::string & algo){
  ////////PLOT THINGS -> This should go into it's own class or at the very least it's own method, TODO

  // CANVAS FOR THE MASS PLOTS
  
  for (int i=0; i<3; i++){//-2; i++){
    std::string alg_idx;
    switch (i)
      {
      case 0:
	alg_idx = "TruthJet";
	break;
      case 1:
	alg_idx = algo;
	break;
      case 2:
	alg_idx = "LeadTruthJet";
	break;
      }
    for (int j=0; j<nPtBins; j++){

      c1[i][j] = new TCanvas();
      
      float maximo = 0.0;
      //Maximum of histograms to set the canvas Y axis properly
      float max1 = qcd_Lead_CA12_mass[i][j]->GetMaximum();
      float max2 = Wprime_Lead_CA12_mass[i][j]->GetMaximum();
      if (max1>max2)
	maximo = max1;
      else
	maximo = max2; 

      qcd_Lead_CA12_mass[i][j]->GetXaxis()->SetTitle("Jet mass [GeV]");
      qcd_Lead_CA12_mass[i][j]->GetYaxis()->SetTitle("Normalised events");
      
      Wprime_Lead_CA12_mass[i][j]->GetXaxis()->SetTitle("Jet mass [GeV]");
      Wprime_Lead_CA12_mass[i][j]->GetYaxis()->SetTitle("Normalised events");
      
      qcd_Lead_CA12_mass[i][j]->SetLineColor(4);
      Wprime_Lead_CA12_mass[i][j]->SetLineColor(2);
      Wprime_Lead_CA12_mass[i][j]->SetMaximum(maximo*1.40);
      
      //Plotting 
      Wprime_Lead_CA12_mass[i][j]->SetMarkerSize(0.0);
      //Wprime_Lead_CA12_mass[j]->DrawCopy("E");
      Wprime_Lead_CA12_mass[i][j]->DrawCopy("hist");
      TH1F * Wprime_mass_err = (TH1F*) Wprime_Lead_CA12_mass[i][j]->Clone();
      Wprime_mass_err->SetFillStyle(3001);
      Wprime_mass_err->SetMarkerSize(0.0);
      Wprime_mass_err->SetFillColor(2);
      Wprime_mass_err->DrawCopy("E2same");
      qcd_Lead_CA12_mass[i][j]->SetMarkerSize(0.0);
      qcd_Lead_CA12_mass[i][j]->Draw("histsame");
      TH1F * QCD_mass_err = (TH1F*) qcd_Lead_CA12_mass[i][j]->Clone();
      QCD_mass_err->SetFillStyle(3001);
      QCD_mass_err->SetMarkerSize(0.0);
      QCD_mass_err->SetFillColor(4);
      QCD_mass_err->DrawCopy("E2same");

      //Now for all labels and shit
      TLatex texw;
      texw.SetNDC();
      texw.SetTextSize(0.045);
      texw.SetTextFont(72);
      texw.DrawLatex(0.20,0.88,"ATLAS");
      
      TLatex p;
      p.SetNDC();
      p.SetTextFont(42);
      p.SetTextSize(0.045);
      p.SetTextColor(kBlack);
      p.DrawLatex(0.31,0.88,"Internal Simulation");
      
      TLatex p2;
      p2.SetNDC();
      p2.SetTextFont(42);
      p2.SetTextSize(0.035);
      p2.SetTextColor(kBlack);
      p2.DrawLatex(0.20,0.82,std::string(algorithms.AlgoListN[alg_idx]).c_str());

    
      TLatex p3;
      p3.SetNDC();
      p3.SetTextFont(42);
      p3.SetTextSize(0.035);
      p3.SetTextColor(kBlack);
      p3.DrawLatex(0.20,0.77,pTbinsN[j]);
     
      TLegend *leg = new TLegend(0.82,0.93,0.93,0.76);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextFont(42);
      leg->SetTextSize(0.040);
      leg->AddEntry( qcd_Lead_CA12_mass[i][j] , "QCD" , "l" );
      leg->AddEntry( Wprime_Lead_CA12_mass[i][j] , "Signal" , "l" );
      leg->Draw();

      c1[i][j]->SaveAs(std::string("jetMass_"+algorithms.AlgoList[alg_idx]+pTbins[j]+".eps").c_str());
      
    }
  }
  
  /////////END OF PLOTTING THINGS
  
  ///// Fine pT bins plots
  
  for (int i=0; i<3; i++){//-2; i++){
  //for (int ii=0; ii<nAlgos-1; ii++){//-2; i++){
  //int i = algoMap[ii];
    std::string alg_idx;
    switch (i)
      {
      case 0:
	alg_idx = "TruthJet";
	break;
      case 1:
	alg_idx = algo;
	break;
      case 2:
	alg_idx = "LeadTruthJet";
	break;
      }
    for (int j=0; j<nFineBins; j++){
      
      c3[i][j] = new TCanvas();
      
      float maximo = 0.0;
      //Maximum of histograms to set the canvas Y axis properly
      float max1 = qcd_finePtBin_mass[i][j]->GetMaximum();
      float max2 = Wprime_finePtBin_mass[i][j]->GetMaximum();
      if (max1>max2)
	maximo = max1;
      else
	maximo = max2; 

      qcd_finePtBin_mass[i][j]->GetXaxis()->SetTitle("Jet mass [GeV]");
      qcd_finePtBin_mass[i][j]->GetYaxis()->SetTitle("Normalised events");
      
      Wprime_finePtBin_mass[i][j]->GetXaxis()->SetTitle("Jet mass [GeV]");
      Wprime_finePtBin_mass[i][j]->GetYaxis()->SetTitle("Normalised events");
      
      qcd_finePtBin_mass[i][j]->SetLineColor(4);
      Wprime_finePtBin_mass[i][j]->SetLineColor(2);
      Wprime_finePtBin_mass[i][j]->SetMaximum(maximo*1.40);
      
      //Plotting 
      Wprime_finePtBin_mass[i][j]->SetMarkerSize(0.0);
      //Wprime_Lead_CA12_mass[j]->DrawCopy("E");
      Wprime_finePtBin_mass[i][j]->DrawCopy("hist");
      TH1F * Wprime_mass_err = (TH1F*) Wprime_finePtBin_mass[i][j]->Clone();
      Wprime_mass_err->SetFillStyle(3001);
      Wprime_mass_err->SetMarkerSize(0.0);
      Wprime_mass_err->SetFillColor(2);
      Wprime_mass_err->DrawCopy("E2same");
      qcd_finePtBin_mass[i][j]->SetMarkerSize(0.0);
      qcd_finePtBin_mass[i][j]->Draw("histsame");
      TH1F * QCD_mass_err = (TH1F*) qcd_finePtBin_mass[i][j]->Clone();
      QCD_mass_err->SetFillStyle(3001);
      QCD_mass_err->SetMarkerSize(0.0);
      QCD_mass_err->SetFillColor(4);
      QCD_mass_err->DrawCopy("E2same");

      //Now for all labels and shit
      TLatex texw;
      texw.SetNDC();
      texw.SetTextSize(0.045);
      texw.SetTextFont(72);
      texw.DrawLatex(0.20,0.88,"ATLAS");
      
      TLatex p;
      p.SetNDC();
      p.SetTextFont(42);
      p.SetTextSize(0.045);
      p.SetTextColor(kBlack);
      p.DrawLatex(0.31,0.88,"Internal Simulation");
      
      TLatex p2;
      p2.SetNDC();
      p2.SetTextFont(42);
      p2.SetTextSize(0.035);
      p2.SetTextColor(kBlack);
      p2.DrawLatex(0.20,0.82,std::string(algorithms.AlgoListN[alg_idx]).c_str());

    
      TLatex p3;
      p3.SetNDC();
      p3.SetTextFont(42);
      p3.SetTextSize(0.035);
      p3.SetTextColor(kBlack);
      p3.DrawLatex(0.20,0.77,pTFinebinsN[j]);
     
      TLegend *leg = new TLegend(0.82,0.93,0.93,0.76);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextFont(42);
      leg->SetTextSize(0.040);
      leg->AddEntry( qcd_finePtBin_mass[i][j] , "QCD" , "l" );
      leg->AddEntry( Wprime_finePtBin_mass[i][j] , "Signal" , "l" );
      leg->Draw();

      c3[i][j]->SaveAs(std::string("jetMass_fineBins_"+algorithms.AlgoList[alg_idx]+finePtBins[j]+".eps").c_str());

    }
    }
} // end makePlots()

/*
 * Create plots in pT bins.
 * @param algo The name of the grooming algorithm being run.
 */
void makePtPlots(std::string & algo){

  // Now I want to make plots which enclose all those information, for each pT bin

  for (int i=0; i<nPtBins; i++){
    hMassLow[i] = new TH1F("MassLow_"+pTbins[i],"MassLow_"+pTbins[i],1, 0, nAlgos-2);
    hMassHigh[i] = new TH1F("MassHigh_"+pTbins[i],"MassHigh_"+pTbins[i],1, 0, nAlgos-2);
    hMPV[i] = new TH1F("MPV_"+pTbins[i],"MPV_"+pTbins[i],1, 0, nAlgos-2);
    hWmass[i] = new TH1F("Wmass_"+pTbins[i],"Wmass_"+pTbins[i],1, 0, nAlgos-2);
    hQCDeff[i] = new TH1F("QCDeff_"+pTbins[i],"QCDeff_"+pTbins[i], 1, 0, nAlgos-2);

  }

  // for (int i=0; i<nAlgos-2; i++){
  //for now only the algos I already have: 0,1,2,3,4,5,6,7,8,9 
  //for (int ii=0; ii<nAlgos; ii++){  //-2; ii++){  
  //int i = algoMap[ii];
    for (int j=0; j<nPtBins; j++){
      hMassLow[j]->SetBinContent(1,BottomEdgeMassWindow[j]);
      hMassHigh[j]->SetBinContent(1,TopEdgeMassWindow[j]);
      hMPV[j]->SetBinContent(1,myMPV[j]);
      hQCDeff[j]->SetBinContent(1,QCDfrac[j]);
      hWmass[j]->SetBinContent(1,80.0);
    }
    //}
  
  //NOW PLOT WINDOW WIDTH VS PT

    windowsVsPt = new TH1F(std::string("windowVsPt_"+algorithms.AlgoList[algo]).c_str(), std::string("windowVsPt_"+algorithms.AlgoList[algo]).c_str(), 10, 0.0, 2500.);
    windowsVsPt->Sumw2();
    windowsVsPt->SetBinContent(1,WidthMassWindow_finePt[0]);
    windowsVsPt->SetBinContent(2,WidthMassWindow_finePt[1]);
    windowsVsPt->SetBinContent(3,WidthMassWindow_finePt[2]);
    windowsVsPt->SetBinContent(4,WidthMassWindow_finePt[3]);
    windowsVsPt->SetBinContent(5,WidthMassWindow_finePt[4]);
    windowsVsPt->SetBinContent(6,WidthMassWindow_finePt[5]);
    windowsVsPt->SetBinContent(7,WidthMassWindow_finePt[6]);
    windowsVsPt->SetBinContent(8,WidthMassWindow_finePt[7]);
    windowsVsPt->SetBinContent(9,WidthMassWindow_finePt[8]);
    windowsVsPt->SetBinContent(10,WidthMassWindow_finePt[9]);

    
    hQCDeff_finePt = new TH1F(std::string("QCD_fineBin_"+algorithms.AlgoList[algo]).c_str(), std::string("QCD_fineBin_"+algorithms.AlgoList[algo]).c_str(),10,0.0, 2500.);
    
    for (int k=1; k<nFineBins-1; k++){
      
      hQCDeff_finePt->SetBinContent(k, QCDfrac_finePt[k-1]);
      
    }

    for (int k=0; k<nFineBins; k++){
      windowsVsPt->SetBinError(k+1,2.0);

    }
    //}
  

    TCanvas *wVsPt;//[nAlgos]; // was -2
    TPad *pad11;//[nAlgos];
    TPad *pad22;//[nAlgos];
	     
    //for (int ii=0; ii<nAlgos; ii++){//ii = 1; ...... -2; ii++){
    //int i = algoMap[ii];
    wVsPt = new TCanvas("wVsPt", "wVsPt", 700, 550);
    
    pad11 = new TPad("pad1","pad1",0,0.4,1,1);
    pad11->SetBottomMargin(0.009);
    pad11->Draw();
    pad22 = new TPad("pad2","pad2",0,0,1,0.4);
    pad22->SetTopMargin(0.009);
    pad22->SetBottomMargin(0.25);
    pad22->Draw();

    pad11->cd();

    windowsVsPt->GetYaxis()->SetTitle("Mass width (68% W eff) [GeV]");  
    windowsVsPt->GetYaxis()->SetTitleSize(0.06);  
    windowsVsPt->GetYaxis()->SetTitleOffset(0.9);

    windowsVsPt->GetXaxis()->SetTitle("Jet CA12 Truth p_{T} [GeV]");  
    windowsVsPt->SetMaximum(windowsVsPt->GetMaximum()*1.3);
    windowsVsPt->Draw("P");

    TLatex texw;
    texw.SetNDC();
    texw.SetTextSize(0.065);
    texw.SetTextFont(72);
    texw.DrawLatex(0.20,0.88,"ATLAS");
    
    TLatex p;
    p.SetNDC();
    p.SetTextFont(42);
    p.SetTextSize(0.065);
    p.SetTextColor(kBlack);
    p.DrawLatex(0.31,0.88,"Internal Simulation");
    
    TLatex p2;
    p2.SetNDC();
    p2.SetTextFont(42);
    p2.SetTextSize(0.055);
    p2.SetTextColor(kBlack);
    p2.DrawLatex(0.20,0.80,std::string(algorithms.AlgoListN[algo]).c_str());
    

    

    pad22->cd();
    hQCDeff_finePt->SetLineColor(2);
    //hQCDeff_finePt->SetMaximum(1.1*hQCDeff_finePt->GetMaximum());
    hQCDeff_finePt->SetMinimum(0.05);
    hQCDeff_finePt->SetMaximum(0.5);
    hQCDeff_finePt->Draw();


    hQCDeff_finePt->GetYaxis()->SetTitle("QCD fraction");
    hQCDeff_finePt->GetXaxis()->SetTitle("Jet CA12 Truth p_{T} [GeV]");  
    hQCDeff_finePt->GetXaxis()->SetTitleSize(0.09);
    hQCDeff_finePt->GetXaxis()->SetLabelSize(0.08);
    hQCDeff_finePt->GetYaxis()->SetLabelSize(0.08);
    hQCDeff_finePt->GetYaxis()->SetTitleSize(0.09);
    hQCDeff_finePt->GetYaxis()->SetTitleOffset(0.55);
    

    wVsPt->SaveAs(std::string("windowsVsPt_"+algorithms.AlgoList[algo]+".eps").c_str());
    
    //}

  
  /////PLOT THOSE INTO APPROPRIATE FORMAT

  for (int i=0; i<nPtBins; i++){
    
    c2[i] = new TCanvas("canvas", "canvas", 700,550);
    
    pad1[i] = new TPad("pad1","pad1",0,0.55,1,1);
    pad1[i]->SetBottomMargin(0.009);
    pad1[i]->Draw();
    pad2[i] = new TPad("pad2","pad2",0,0,1,0.55);
    pad2[i]->SetTopMargin(0.009);
    pad2[i]->SetBottomMargin(0.5);
    pad2[i]->Draw();
    
    pad1[i]->cd();
    hMassHigh[i]->SetLineColor(4);
    hMassHigh[i]->GetYaxis()->SetTitle("Jet mass [GeV]");
    hMassHigh[i]->GetYaxis()->SetTitleSize(0.09);
    hMassHigh[i]->GetYaxis()->SetTitleOffset(0.5);
    hMassHigh[i]->SetMaximum((1.3*hMassHigh[i]->GetMaximum()));
    hMassHigh[i]->SetMinimum(0.0);
    hMassHigh[i]->Draw();
    hWmass[i]->SetLineStyle(2);
    hWmass[i]->Draw("same");
    hMPV[i]->Draw("same");
    hMassLow[i]->SetLineColor(2);
    hMassLow[i]->Draw("same");

    //Now for all labels and shit
    TLatex texw;
    texw.SetNDC();
    texw.SetTextSize(0.1);
    texw.SetTextFont(72);
    texw.DrawLatex(0.52,0.85,"ATLAS");
    
    TLatex p;
    p.SetNDC();
    p.SetTextFont(42);
    p.SetTextSize(0.1);
    p.SetTextColor(kBlack);
    p.DrawLatex(0.65,0.85,"Internal Simulation");
    
    //   p2.DrawLatex(0.20,0.82,AlgoListN[i]);
    
    TLatex p3;
    p3.SetNDC();
    p3.SetTextFont(42);
    p3.SetTextSize(0.07);
    p3.SetTextColor(kBlack);
    p3.DrawLatex(0.65,0.72,pTbinsN[i]);
    
    
    pad2[i]->cd();
    hQCDeff[i]->SetLineColor(3);
    hQCDeff[i]->SetMaximum(1.3*hQCDeff[i]->GetMaximum());
    hQCDeff[i]->Draw();


    hQCDeff[i]->GetYaxis()->SetTitle("QCD fraction");
    hQCDeff[i]->GetYaxis()->SetTitleSize(0.07);
    hQCDeff[i]->GetYaxis()->SetTitleOffset(0.68);

    for (int k=1; k<=hQCDeff[i]->GetXaxis()->GetNbins(); k++){
      hQCDeff[i]->GetXaxis()->SetBinLabel(k,std::string(algorithms.binLabel[algo/*k-1*/]).c_str());
      hQCDeff[i]->LabelsOption("v");
      
    }

    TLegend *leg = new TLegend(0.65,0.88,0.93,0.97);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.040);
    leg->SetNColumns(2);
    leg->AddEntry( hMassHigh[i], "Upper bound" , "l" );
    leg->AddEntry( hMPV[i] , "MPV" , "l" );
    leg->AddEntry( hMassLow[i] , "Lower bound" , "l" );
    leg->AddEntry( hQCDeff[i] , "QCD fraction" , "l" );
    leg->Draw();
    
    
    c2[i]->SaveAs("massWindows_"+pTbins[i]+".eps");
    
    
    
  }
  
  
  ///end of fine pT bins plots
  

  ///// SAVE THINGS TO HISTOGRAM FILE
  
  TFile f("MunichPlots.root","recreate");

  //Usage:
  //TH1D *S = new TH1D();
  //TH1D *B = new TH1D();
  //...Fill S and B with signal (S) and background (B)
  //Need S and B to have the same number of bins!
  //TH1D *curve = new TH1D("","",M,0,1); where M is however fine you want the ROC curve binning to be.
  //ROC(S,B,curve);
  int type = 1;
  
  //for (int ii=0; ii<nAlgos; ii++){//-2; ii++){
  //int i = algoMap[ii];
    windowsVsPt->Write();
    //}


    pTweights->Write();
    qcd_PtReweight->Write();
    Wp_PtReweight->Write();
    makeROC(type, Wp_PtReweight, qcd_PtReweight, (*PtReweight_curve), "PTReweight", true);
    //PtReweight_curve->Write();
    Wprime_Lead_CA12_scaled_pt->Write();

    //Lead_CA12_scaled_pt_curve->Write();



    for (int i=0; i<3; i++){//-1; ii++){
    //int i = algoMap[ii];
    
    qcd_Lead_CA12_pt[i]->Write();
    Wprime_Lead_CA12_pt[i]->Write();
    makeROC(type, Wprime_Lead_CA12_pt[i],qcd_Lead_CA12_pt[i],(*Lead_CA12_pt_curve)[i], Wprime_Lead_CA12_pt[i]->GetName(), true);
    makeROC(type, Wprime_Lead_CA12_scaled_pt, qcd_Lead_CA12_pt[i], (*Lead_CA12_scaled_pt_curve)[i], Wprime_Lead_CA12_scaled_pt->GetName(), true);// Is his right?!
    //Lead_CA12_pt_curve->Write();


    for (int j=0; j<nPtBins; j++){
      qcd_Lead_CA12_mass[i][j]->Write();
      Wprime_Lead_CA12_mass[i][j]->Write();
      makeROC(type, Wprime_Lead_CA12_mass[i][j], qcd_Lead_CA12_mass[i][j], (*Lead_CA12_mass_curve[i])[j], Wprime_Lead_CA12_mass[i][j]->GetName(), true);
      //Lead_CA12_mass_curve[j]->Write();
    }

    for (int j=0; j<nFineBins; j++){
      qcd_finePtBin_mass[i][j]->Write();
      Wprime_finePtBin_mass[i][j]->Write();
      makeROC(type, Wprime_finePtBin_mass[i][j], qcd_finePtBin_mass[i][j], (*finePtBin_mass_curve[i])[j], Wprime_finePtBin_mass[i][j]->GetName(), true);
      //finePtBin_mass_curve[j]->Write();
    }
    }
  
  for (int i=0; i<nPtBins; i++){
    hMassLow[i]->Write();
    hMassHigh[i]->Write();
    hMPV[i]->Write();
    hQCDeff[i]->Write();    
    hWmass[i]->Write();
  }

  //Wprime_Lead_CA12_scaled_pt[7]->Write();
  //Wprime_Lead_CA12_scaled_pt[8]->Write();
  //pTweights->Write();

  f.Close();

} // end makePtPlots()

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
  double mass_max = 0.0;
  double mass_min = 0.0;

  bool signal = false;

  // need branches for specific algorithms because otherwise we end up with a million branches that we don't need
  vector<std::pair<std::string,int> > branches;

  std::string i = algorithm;//algoMap[ii];
    std::cout << "Doing mass window plots for " << algorithms.AlgoNames[i] << std::endl;

    std::string groomAlgoIndex = i;
    std::string prefix = "";
    TObjArray * brancharray = inputTChain[0]->GetListOfBranches();
    if (branchesFile != "")
      {
	std::cout << "branchesFile: " << branchesFile << std::endl;
	branches = getListOfJetBranches(branchesFile, brancharray);
      }
    else
      {
	std::cout << "branchesFile not defined" << std::endl;
	branches = getListOfJetBranches(algorithms.AlgoNames[i], brancharray);
      }
    // loop through the different pt bins. j == 0 is the inclusive one
    for (int j=0; j<1; j++){//nPtBins; j++){

      signal = false;
      std::stringstream ss;
      ss << algorithms.AlgoNames[i] << "_" << pTbins[j];
      // loop through background and signal
      for (int k = 0; k < 2; k++)
	{
	  signal = k == 1? false : true;
	  // Initialise all the vectors to.. something 
	  initVectors();

	  // set up the mass window
	  mass_max = TopEdgeMassWindow[j];
	  mass_min = BottomEdgeMassWindow[j];
	  
	  // set the tchain pointer to the signal or background sample
	  int tchainIdx = signal ? sampleType::SIGNAL : sampleType::BACKGROUND;
	  std::stringstream ss_fname; // store the name of the output file and include the i and j indices!
	  std::string bkg = signal ? "sig": "bkg";
	  ss_fname << ss.str() << "_" << bkg << ".root";
	  boost::filesystem::path dir(std::string(algorithms.AlgoNames[i]+fileid_global));
	  boost::filesystem::create_directory(dir);

	  inputTChain[tchainIdx]->GetEntries();

	  // turn off the branches we're not interested in
	  for (std::vector<std::pair<string,int > >::iterator it = branches.begin(); it != branches.end(); it++)
	    {

	      if ((*it).second == 0)
		{
		  inputTChain[tchainIdx]->SetBranchStatus((*it).first.c_str(),0);
		}
	      
	    }

	  setJetsBranches(inputTChain[tchainIdx], algorithms.AlgoNames[i], i); //set all of the branches for the output tree for the jets	  

	  // output file
	  TFile * outfile = new TFile(std::string(algorithms.AlgoNames[i]+fileid_global+"/"+ss_fname.str()).c_str(),"RECREATE");   
	  TTree * outTree = new TTree(treeName.c_str(),treeName.c_str());

	  // reset all of the output variables so that they have default values
	  resetOutputVariables();
	  // set up all of the branches for the output tree
	  setOutputBranches(outTree, algorithms.AlgoNames[i], i);

	  // if we are going to add the subjet branches we set them up here
	  if (subjetscalc || subjetspre)
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
	  // add branches for runnumber, mc_event_number etc.
>>>>>>> Updated weights file.  TLV works. Lots of protection put in against reading branches that are not valid.  Still need to add the option to turn off event and lepton selection tlection to do jet grooming algorithm analysis.
	  addInfoBranches(outTree);
>>>>>>> Fully commented.

<<<<<<< HEAD
>>>>>>> Some general cleanup. Set vars to -999 as default so that we can see when something doesn't get set. Particularly important for Topo jets which are not often found.
=======
	  // get number of entries to run over in the tree
>>>>>>> Updated weights file.  TLV works. Lots of protection put in against reading branches that are not valid.  Still need to add the option to turn off event and lepton selection tlection to do jet grooming algorithm analysis.
	  long entries = (long)inputTChain[tchainIdx]->GetEntries();

	  // Setting up this with a high limit of 3.5 TeV so we don't miss anything.  Lots of bins - 200, so we can
	  // do a lot of tuning of the scale factor regions later on!
	  pt_reweight = new TH1F(std::string("pt_reweight"+bkg).c_str(),std::string("pt_reweight_"+bkg).c_str(), 200, 0, 3500);
	  
	  // set up the number of events in total
	  NEvents = entries;
	  // we keep track of the weighted (by mc_event_weight) number of events too
	  NEvents_weighted.clear();

	  std::cout << "total entries: " << entries << std::endl;
	  
	  // mass variable
	  double mass = 0;
	  for (long n = 0; n < entries; n++)
	    {
	      inputTChain[tchainIdx]->GetEntry(n);

	      setSelectionVectors(); // again reading from tree, but just setting jet_pt_x whatever
	      // which makes it easier to use later on.
	      resetOutputVariables(); 
	     

	      if (n%1000==0)
		std::cout << "Entry: "<< n << " / " << entries <<  std::endl;

	      if (NEvents_weighted.find(mc_channel_number) != NEvents_weighted.end())
		NEvents_weighted[mc_channel_number] += mc_event_weight;
	      else
		NEvents_weighted[mc_channel_number] = mc_event_weight;


	      // what about reclustered jets?! argggg
	      int chosenLeadTruthJetIndex=-99;
	      int chosenLeadTopoJetIndex=-99;

	       if ((*jet_pt_truth)[0] / 1000.0 < 100)
		{
       		  continue;
		}
	      
	      if (algorithms.AlgoType[groomAlgoIndex].find("truthmatch") == std::string::npos) // check the pt is in the correct bin
		{
		  if ((*jet_pt_truth)[0]/1000.0 < ptrange[j].first || (*jet_pt_truth)[0]/1000.0 > ptrange[j].second)
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
	      
	      if (algorithms.AlgoType[groomAlgoIndex].find("truthmatch") != std::string::npos)
		{
		  //check which is my CA12 ungroomed reco jet with EMfrac<0.99 and |eta|<1.2
		  //loop over all topo jets and get the leading with cuts
		  bool hasTopoJet=false;
		
		  for (int jet_i=0; jet_i<(*jet_pt_topo).size(); jet_i++)
		    {
		      float emfractmp = 0.5;
		      if (xAODemfrac)
			emfractmp = (*jet_emfrac_topo)[jet_i];
		      if (!hasTopoJet && emfractmp<0.99 && fabs((*jet_eta_topo)[jet_i])<1.2) 
			{
			  chosenLeadTopoJetIndex=jet_i;
			  hasTopoJet=true;
			} 
		    } // end loop over jet_pt_topo
		  
		  //now match truth jet with the chosen topo jet
		  if (hasTopoJet)
		    {
		      for (int jet_i=0; jet_i<(*jet_pt_truth).size(); jet_i++)
			{
			  if (chosenLeadTruthJetIndex<0 && DeltaR((*jet_eta_topo)[chosenLeadTopoJetIndex],(*jet_phi_topo)[chosenLeadTopoJetIndex],(*jet_eta_truth)[jet_i],(*jet_phi_truth)[jet_i])<0.9)
			    { 
			      chosenLeadTruthJetIndex=jet_i;
			    }	  
			}	// end loop over jet_pt_truth
		    } // end if(hasTopoJet)
		} // end if(groomAlgoIndex==0)
	      
	      chosenLeadTruthJetIndex = algorithms.AlgoType[groomAlgoIndex].find("truthmatch") != std::string::npos ? chosenLeadTruthJetIndex : 0;

	      // veto the event if we have no good truth jets.
	      if (chosenLeadTruthJetIndex < 0)
		continue;
	      //Now I have which events to make my pt reweight with, and to match to, etc
	      int chosenLeadGroomedIndex=-99;
	      for (int jet_i=0; jet_i<(*jet_pt_groomed).size(); jet_i++)
		{
		  float emfractmp = 0.5;
		  if (xAODemfrac)
		    emfractmp = (*jet_emfrac_groomed)[jet_i];
		  if (chosenLeadGroomedIndex<0 && DeltaR((*jet_eta_truth)[chosenLeadTruthJetIndex],(*jet_phi_truth)[chosenLeadTruthJetIndex],(*jet_eta_groomed)[jet_i],(*jet_phi_groomed)[jet_i])<0.9 && emfractmp<0.99 && fabs((*jet_eta_groomed)[jet_i])<1.2)
		    {
		      chosenLeadGroomedIndex=jet_i;
		      continue;
		    }     
		} // end loop over jet_pt_groomed
	      
	      if (chosenLeadGroomedIndex == -99) // failed selection
		{
		  continue;	      
		}
	      
	      leadGroomedIndex = chosenLeadGroomedIndex;
	      leadTruthIndex = chosenLeadTruthJetIndex;
	      leadTopoIndex = chosenLeadTopoJetIndex;

	      // set the mass to the leading groomed jet
	      mass = (*var_m_vec[2])[leadGroomedIndex]/1000.0 ;

	      // if we are not applying a mass window we do not apply any mass cuts
	      if (applyMassWindow && (mass > mass_max && mass < mass_min))
		{
		  continue;
		}
	      // set the mass of the two leptons and leading groomed jet
	      setLLJMass(leadGroomedIndex);

	      // leading subjet index
	      int lead_subjet = 0;

	      // if we are calculating the subjet variables: massdrop and momentum balance
	      if (subjetscalc)
		{

		  std::vector<int> subjet_idx = (*subjet_index).at(leadGroomedIndex); // only groomed ones.....
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
		    
		      double subjet_mass =  (*var_subjets_m_vec)[subjet_leading.first];
		    
		      // calculate the mass drop - leading subjet/mass
		      double mu12 = subjet_mass/(mass*1000);
		      var_massdrop = mu12;
		    
		      // momentum balance
		      double dRsub12 = DeltaR (eta_1, phi_1, eta_2, phi_2);
		      Float_t yt = (pt2*dRsub12)/(mass*1000);		    
		      yt*=yt;
		      var_yt = yt;
		    }
		} // end if(subjets)
	      // if we have the mass fraction and momentum balance in the input file already
	      else if (subjetspre)
		{
		  var_massdrop = var_massFraction_vec;
		  var_yt = var_ktycut2_vec;
		}

	      // tau21
	      var_Tau21[0]=(*var_Tau2_vec[0])[leadTruthIndex]/(*var_Tau1_vec[0])[leadTruthIndex];
	      var_Tau21[1]=(*var_Tau2_vec[1])[leadTopoIndex]/(*var_Tau1_vec[1])[leadTopoIndex];
	      var_Tau21[2]=(*var_Tau2_vec[2])[leadGroomedIndex]/(*var_Tau1_vec[2])[leadGroomedIndex];
	    
	      // set up tauwta variables and zcut12
	      if (useBranch(string("TauWTA2TauWTA1"),true) && useBranch(string("TauWTA2"), true) && useBranch(string("TauWTA1"), true) )
		{
		  // need to check if these branches are on....
		  if (var_TauWTA1_vec[0] != NULL && var_TauWTA2_vec[0] != NULL)
		    var_TauWTA2TauWTA1[0]=(*var_TauWTA2_vec[0])[leadTruthIndex]/(*var_TauWTA1_vec[0])[leadTruthIndex];
		  if (var_TauWTA1_vec[2] != NULL && var_TauWTA2_vec[2] != NULL)
		    var_TauWTA2TauWTA1[2]=(*var_TauWTA2_vec[2])[leadGroomedIndex]/(*var_TauWTA1_vec[2])[leadGroomedIndex];
		  if (leadTopoIndex == -99)
		    var_TauWTA2TauWTA1[1]=-99;
		  else if (var_TauWTA1_vec[1] != NULL && var_TauWTA2_vec[1] != NULL)
		    var_TauWTA2TauWTA1[1]=(*var_TauWTA2_vec[1])[leadTopoIndex]/(*var_TauWTA1_vec[1])[leadTopoIndex];
		}
		
	    
	      // make sure all of the other output variables have their values set
	      setOutputVariables(leadTruthIndex, leadTopoIndex, leadGroomedIndex, lead_subjet, algorithms.AlgoNames[i] , i);
	      outTree->Fill();
	      pt_reweight->Fill((*jet_pt_truth)[chosenLeadTruthJetIndex]/1000.0);

	    } // end loop over nentries

	  // write the rweight th1f things to the outfile...
	  // calculate the overall reweight with a new histogram!
	  outTree->GetCurrentFile()->Write();
	  pt_reweight->Write();
	  pt_reweight_arr[tchainIdx] = (TH1F*)pt_reweight->Clone();
	  // stupid clone method needs this so that it doesn't delete this histo when closing the file
	  // http://root.cern.ch/phpBB3/viewtopic.php?f=3&t=11486
	  pt_reweight_arr[tchainIdx]->SetDirectory(0);

	  outTree->GetCurrentFile()->Close();

	  std::stringstream ss2; // store the name of the output file and include the i and j indices!
	  std::string bkg2 = signal ? "sig": "bkg";
	  ss2 << algorithms.AlgoNames[i] << fileid_global << "/" << ss.str() << "_" << bkg2 << ".nevents";
	  ofstream ev_out(ss2.str());
	  // write out the weighted events to a text file
	  for (std::map<int,float>::iterator it = NEvents_weighted.begin(); it!= NEvents_weighted.end(); it++)
	    ev_out << it->first << "," << NEvents_weighted[it->first] << std::endl;
	  ev_out.close();
	  delete outfile;

	} // end loop of datatype

      // create the reweighting file
      std::string fname = algorithms.AlgoNames[i]+fileid_global+"/" + ss.str() + ".ptweights";
      createPtReweightFile(pt_reweight_arr[sampleType::BACKGROUND], pt_reweight_arr[sampleType::SIGNAL], fname);
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
std::pair<int,int> getTwoLeadingSubjets(std::vector<int> & jet_idx, std::vector<float> * subjet_pt)
{
  // jet_idx contains the indices in the subjet vector for the subjets of jet i
  double max_pt = 0;
  double sec_pt = 0;
  int max = -1;
  int sec = -1;

  for (int i = 0; i < jet_idx.size(); i++)
    {
      int idx = jet_idx[i];
      if (max_pt < (*subjet_pt)[idx])
	{
	  sec_pt = max_pt;
	  sec = max;
	  max_pt = (*subjet_pt)[idx];
	  max = idx;
	}
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
 * @param brancharray A TObjArray pointer containing all of the branches in the input tree
 * @return vector of pairs containing all of the branches and if they should be on or off.
 */
vector<std::pair<std::string,int> > getListOfJetBranches(std::string &algorithm, TObjArray * brancharray)
{
  // Ideally we want this stored in an XML file, but for now it'll have to be a standard text file because I'm short on time!
  vector<pair<string,int> > branches;
  
  // reset branchmap
  branchmap.clear();
  std::string filename = algorithm.find("branches.txt") == std::string::npos ? algorithm+"_branches.txt" : algorithm;

  // open input file
  ifstream in(filename);
  string line;
  // keep a record of all the branches in the config file in a map
  while (getline(in, line))
    {
      trim(line);
      branchmap[line] = 1;
    }

  // need to add some essential ones in case they get forgotten in that config file :)
  branchmap["RunNumber"] = 1;
  branchmap["mc_channel_number"] = 1;
  branchmap["vxp_n"] = 1;
  branchmap["nVertices"] = 1;
  branchmap["averageIntPerXing"] = 1;
  branchmap["mc_event_weight"] = 1;

  // need to add leptons
  branchmap["electrons"] = 1;

  branchmap["el_pt"] = 1;
  branchmap["el_eta"] = 1;
  branchmap["el_phi"] = 1;
  branchmap["el_ptcone20"] = 1;
  branchmap["el_etcone20"] = 1;

  branchmap["muons"] = 1;

  branchmap["mu_pt"] = 1;
  branchmap["mu_eta"] = 1;
  branchmap["mu_phi"] = 1;
  branchmap["mu_ptcone20"] = 1;
  branchmap["mu_etcone20"] = 1;
  branchmap["mu_charge"] = 1;
  in.close();

  // loop through all of the branches in the tree and check if they exist in the branch map
  TIter it(brancharray);
  TBranch * tb;
  while (tb = (TBranch *) it.Next())
    {
      std::string name = tb->GetName();
      trim(name);
      if (branchmap.find(name) != branchmap.end() )
	branches.push_back(make_pair(name, 1));
      else
	branches.push_back(make_pair(name, 0));
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
 * Check if a branch is in the branchmap so we don't try to set branch addresses we aren't using
 *
 * @param branch Name of the branch to check
 * @param partialmatch Only look for substrings of keys that match the string
 *
 * @return bool indicating if it is in the branchmap
 */
bool useBranch(std::string branch, bool partialmatch)
{
  // look for the key in the branch map
  if (!partialmatch)
    {
      if (branchmap.find(branch) != branchmap.end())
	return true;
    }
  // look for a partial match
  else
    {
      int brlen = branch.length();
      for (std::map<std::string, int>::iterator it = branchmap.begin(); it != branchmap.end(); it++)
	{
	  if (it->first.length() >= brlen && it->first.compare(it->first.length()-brlen, brlen, branch) == 0)
	    return true;
	}
    }

  return false;
}

/*
 * Adds essential output branches to the output file.
 *
 * @param tree Pointer to the output TTree.
 */
void addInfoBranches(TTree * tree)
{
  tree->Branch("leadTruthIndex",&leadTruthIndex, "leadTruthIndex/I");
  tree->Branch("leadTopoIndex",&leadTopoIndex, "leadTopoIndex/I");
  tree->Branch("leadGroomedIndex",&leadGroomedIndex, "leadGroomedIndex/I");

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
 * @param list Reference to a TObjArray pointer that contains a list of all branches in the tree.
 * @param vec Reference to a vector<TLV> pointer which will be used to read in from the file.
 * @param branch The name of the branch we are setting the address for.
 */
void setVector(TChain *& tree, TObjArray *& list, vector<TLorentzVector> *& vec, string branch)//const char * branch)
{
  // first check to see if the branch is actually in the tree, then check to see if we are interested in using it
  if (list->FindObject(branch.c_str()) && useBranch(branch))
    tree->SetBranchAddress(branch.c_str(), &vec);
  else
    {
      std::cout << "missing branch " << branch << ", might cause unexpected behaviour, removing from branchmap" << std::endl;
      branchmap.erase(branch);
    }
}

/*
 * Set up the branch address for a vector<Float_t>.  Overloaded version of the one for vector<TLV>.
 *
 * @param tree Reference to a TChain pointer.  This is the input TTree.
 * @param list Reference to a TObjArray pointer that contains a list of all branches in the tree.
 * @param vec Reference to a vector<Float_t> pointer which will be used to read in from the file.
 * @param branch The name of the branch we are setting the address for.
 */
void setVector(TChain *& tree, TObjArray *& list, vector<Float_t> *& vec, string branch)//const char * branch)
{
  // first check to see if the branch is actually in the tree, then check to see if we are interested in using it
  if (list->FindObject(branch.c_str()) && useBranch(branch))
    tree->SetBranchAddress(branch.c_str(), &vec);
  else
    {
      std::cout << "missing branch " << branch << ", might cause unexpected behaviour, removing from branchmap" << std::endl;
      branchmap.erase(branch);
    }
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
  if (var_Tau1_vec[i] != NULL && var_Tau1_vec[i]->size() > jet)
    var_Tau1_vec[i]->erase(var_Tau1_vec[i]->begin()+jet);
  if (var_Tau2_vec[i] != NULL && var_Tau2_vec[i]->size() > jet)
    var_Tau2_vec[i]->erase(var_Tau2_vec[i]->begin()+jet);
  if (var_Tau3_vec[i] != NULL && var_Tau3_vec[i]->size() > jet)
    var_Tau3_vec[i]->erase(var_Tau3_vec[i]->begin()+jet);
  if (var_WIDTH_vec[i] != NULL && var_WIDTH_vec[i]->size() > jet)
    var_WIDTH_vec[i]->erase(var_WIDTH_vec[i]->begin()+jet);
  if (var_SPLIT12_vec[i] != NULL && var_SPLIT12_vec[i]->size() > jet)
    var_SPLIT12_vec[i]->erase(var_SPLIT12_vec[i]->begin()+jet);
  if (var_SPLIT23_vec[i] != NULL && var_SPLIT23_vec[i]->size() > jet)
    var_SPLIT23_vec[i]->erase(var_SPLIT23_vec[i]->begin()+jet);
  if (var_SPLIT34_vec[i] != NULL && var_SPLIT34_vec[i]->size() > jet)
    var_SPLIT34_vec[i]->erase(var_SPLIT34_vec[i]->begin()+jet);
  if (var_Dip12_vec[i] != NULL && var_Dip12_vec[i]->size() > jet)
    var_Dip12_vec[i]->erase(var_Dip12_vec[i]->begin()+jet);
  if (var_Dip13_vec[i] != NULL && var_Dip13_vec[i]->size() > jet)
    var_Dip13_vec[i]->erase(var_Dip13_vec[i]->begin()+jet);
  if (var_Dip23_vec[i] != NULL && var_Dip23_vec[i]->size() > jet)
    var_Dip23_vec[i]->erase(var_Dip23_vec[i]->begin()+jet);
  if (var_DipExcl12_vec[i] != NULL && var_DipExcl12_vec[i]->size() > jet)
    var_DipExcl12_vec[i]->erase(var_DipExcl12_vec[i]->begin()+jet);
  if (var_PlanarFlow_vec[i] != NULL && var_PlanarFlow_vec[i]->size() > jet)
    var_PlanarFlow_vec[i]->erase(var_PlanarFlow_vec[i]->begin()+jet);
  if (var_Angularity_vec[i] != NULL && var_Angularity_vec[i]->size() > jet)
    var_Angularity_vec[i]->erase(var_Angularity_vec[i]->begin()+jet);
  if (var_QW_vec[i] != NULL && var_QW_vec[i]->size() > jet)
    var_QW_vec[i]->erase(var_QW_vec[i]->begin()+jet);
  if (var_PullMag_vec[i] != NULL && var_PullMag_vec[i]->size() > jet)
    var_PullMag_vec[i]->erase(var_PullMag_vec[i]->begin()+jet);
  if (var_PullPhi_vec[i] != NULL && var_PullPhi_vec[i]->size() > jet)
    var_PullPhi_vec[i]->erase(var_PullPhi_vec[i]->begin()+jet);
  if (var_Pull_C00_vec[i] != NULL && var_Pull_C00_vec[i]->size() > jet)
    var_Pull_C00_vec[i]->erase(var_Pull_C00_vec[i]->begin()+jet);
  if (var_Pull_C01_vec[i] != NULL && var_Pull_C01_vec[i]->size() > jet)
    var_Pull_C01_vec[i]->erase(var_Pull_C01_vec[i]->begin()+jet);
  if (var_Pull_C10_vec[i] != NULL && var_Pull_C10_vec[i]->size() > jet)
    var_Pull_C10_vec[i]->erase(var_Pull_C10_vec[i]->begin()+jet);
  if (var_Pull_C11_vec[i] != NULL && var_Pull_C11_vec[i]->size() > jet)
    var_Pull_C11_vec[i]->erase(var_Pull_C11_vec[i]->begin()+jet);
  if (var_YFilt_vec != NULL && var_YFilt_vec->size() > jet)
    var_YFilt_vec->erase(var_YFilt_vec->begin()+jet);
  
  if (var_ActiveArea_vec[i] != NULL && var_ActiveArea_vec[i]->size() >> jet)
    var_ActiveArea_vec[i]->erase(var_ActiveArea_vec[i]->begin()+jet);
  if (var_VoronoiArea_vec[i] != NULL && var_VoronoiArea_vec[i]->size() >> jet)
    var_VoronoiArea_vec[i]->erase(var_VoronoiArea_vec[i]->begin()+jet);
  if (var_Aplanarity_vec[i] != NULL && var_Aplanarity_vec[i]->size() >> jet)
    var_Aplanarity_vec[i]->erase(var_Aplanarity_vec[i]->begin()+jet);
  if (var_Sphericity_vec[i] != NULL && var_Sphericity_vec[i]->size() >> jet)
    var_Sphericity_vec[i]->erase(var_Sphericity_vec[i]->begin()+jet);
  if (var_ThrustMaj_vec[i] != NULL && var_ThrustMaj_vec[i]->size() >> jet)
    var_ThrustMaj_vec[i]->erase(var_ThrustMaj_vec[i]->begin()+jet);
  if (var_ThrustMin_vec[i] != NULL && var_ThrustMin_vec[i]->size() >> jet)
    var_ThrustMin_vec[i]->erase(var_ThrustMin_vec[i]->begin()+jet);

  if (var_TauWTA1_vec[i] != NULL && var_TauWTA1_vec[i]->size() > jet)
    var_TauWTA1_vec[i]->erase(var_TauWTA1_vec[i]->begin()+jet);
  if (var_TauWTA2_vec[i] != NULL && var_TauWTA2_vec[i]->size() > jet)
    var_TauWTA2_vec[i]->erase(var_TauWTA2_vec[i]->begin()+jet);
  if (var_TauWTA3_vec[i] != NULL && var_TauWTA3_vec[i]->size() > jet)
    var_TauWTA3_vec[i]->erase(var_TauWTA3_vec[i]->begin()+jet);
  if (var_ZCUT12_vec[i] != NULL && var_ZCUT12_vec[i]->size() > jet)
    var_ZCUT12_vec[i]->erase(var_ZCUT12_vec[i]->begin()+jet);
  if (var_ZCUT23_vec[i] != NULL && var_ZCUT23_vec[i]->size() > jet)
    var_ZCUT23_vec[i]->erase(var_ZCUT23_vec[i]->begin()+jet);
  if (var_ZCUT34_vec[i] != NULL && var_ZCUT34_vec[i]->size() > jet)
    var_ZCUT34_vec[i]->erase(var_ZCUT34_vec[i]->begin()+jet);
    
  
} // eraseJet

/*
 * Overlap removal between jets and electrons, removes any (groomed) jet that overlaps with an electron.
 *
 */
void overlapRemoval()
{
  for (int it = 0 ; it < (*var_pt_vec[jetType::GROOMED]).size(); it++)
    {
      for (std::vector<TLorentzVector>::iterator el = var_electrons_vec->begin(); el != var_electrons_vec->end(); el++)
	{
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
  if (!hvtllqq)
    {
      std::cout << "running event selection when analysis is not being done, should not call this method" << std::endl;
      return -1;
    }
  int lepType = leptonType::FAIL;
  // 2 electrons or 2 muons of opposite charge  
  if (var_electrons_vec->size() == 2)
    {
      var_mll = (var_electrons_vec->at(0)+var_electrons_vec->at(1)).M();
      var_ptll = (var_electrons_vec->at(0)+var_electrons_vec->at(1)).Pt();
      lepType = leptonType::ELECTRON;
    }
  else if (var_muons_vec->size() == 2 && (*var_mu_charge_vec)[0]*(*var_mu_charge_vec)[1] == -1)
    {
      var_mll = (var_muons_vec->at(0)+var_muons_vec->at(1)).M();
      var_ptll = (var_muons_vec->at(0)+var_muons_vec->at(1)).Pt();
      lepType = leptonType::MUON;
    }
  else
    {
      // it has failed
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
      // et > 25 GeV
      for (int idx = 0; idx < 2; idx++)
	{
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
	  //std::cout << (*var_el_ptcone20_vec)[idx] << std::endl;
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
      // pt > 25 gev
      for (int idx = 0; idx < 2; idx++)
	{
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
void setJetsBranches(TChain * tree, std::string &groomalgo,  std::string & groomIdx)
{
  std::string samplePrefix = "";
  // get a list of all the branches in the tree
  TObjArray * brancharray = tree->GetListOfBranches();
  bool addLC = false; // we use this variable to indicate whether "LC" should be in the algorithm string
  samplePrefix = algorithms.AlgoPrefix[groomIdx];
  
  if (algorithms.AlgoType[groomIdx].find("recluster") == std::string::npos) // we're doing reclustering
    addLC = true; // just add the LC to the name

  tree->SetBranchAddress("mc_event_weight",&mc_event_weight);
  tree->SetBranchAddress("mc_channel_number", &mc_channel_number);

  // double check to make sure the branch exists in the tree
  if (brancharray->FindObject("RunNumber"))
    tree->SetBranchAddress("RunNumber", &runNumberIn);
  if (brancharray->FindObject("nVertices"))
    tree->SetBranchAddress("nVertices", &nvtxIn);
  if (brancharray->FindObject("averageIntPerXing"))
    tree->SetBranchAddress("averageIntPerXing",&avgIntpXingIn);

  
  //
  // set up the branch addresses for the lepton variables if running hvtllqq analysis
  if (hvtllqq)
    {
      setVector(tree, brancharray, var_electrons_vec, std::string("electrons"));
      setVector(tree, brancharray, var_electronPt_vec, std::string("el_pt"));
      setVector(tree, brancharray, var_electronEta_vec, "el_eta");
      setVector(tree, brancharray, var_electronPhi_vec, "el_phi");
      setVector(tree, brancharray, var_el_ptcone20_vec, "el_ptcone20");
      setVector(tree, brancharray, var_el_etcone20_vec, "el_etcone20");
      
      setVector(tree, brancharray, var_muons_vec, "muons");
      setVector(tree, brancharray, var_muonPt_vec, "mu_pt");
      setVector(tree, brancharray, var_muonEta_vec, "mu_eta");
      setVector(tree, brancharray, var_muonPhi_vec, "mu_phi");
      setVector(tree, brancharray, var_mu_ptcone20_vec, "mu_ptcone20");
      setVector(tree, brancharray, var_mu_etcone20_vec, "mu_etcone20");
      setVector(tree, brancharray, var_mu_charge_vec, "mu_charge");
    }

  // there is only yfilt stored for groomed jets
  setVector(tree, brancharray, var_YFilt_vec, std::string(returnJetType(samplePrefix, groomalgo, addLC,2)+"YFilt") );
  // loop through the truth, toppo and groomed jets and set up the branches for the different variables
  for (int i = 0; i < jetType::MAX; i++) // truth, topo, groomed
    {
      std::string jetString = returnJetType(samplePrefix, groomalgo, addLC,i); //set to truth/ topo/ groomed
      setVector(tree, brancharray, var_E_vec.at(i), std::string(jetString+"E") );
      setVector(tree, brancharray, var_pt_vec.at(i), std::string(jetString+"pt") );
      setVector(tree, brancharray, var_m_vec.at(i), std::string(jetString+"m") );
      setVector(tree, brancharray, var_eta_vec.at(i), std::string(jetString+"eta") );
      setVector(tree, brancharray, var_phi_vec.at(i), std::string(jetString+"phi") );
      setVector(tree, brancharray, var_emfrac_vec.at(i), std::string(jetString+"emfrac") );
      setVector(tree, brancharray, var_Tau1_vec.at(i), std::string(jetString+"Tau1") );
      setVector(tree, brancharray, var_Tau2_vec.at(i), std::string(jetString+"Tau2") );
      setVector(tree, brancharray, var_Tau3_vec.at(i), std::string(jetString+"Tau3") );
      setVector(tree, brancharray, var_WIDTH_vec.at(i), std::string(jetString+"WIDTH") );
      setVector(tree, brancharray, var_SPLIT12_vec.at(i), std::string(jetString+"SPLIT12") );
      setVector(tree, brancharray, var_SPLIT23_vec.at(i), std::string(jetString+"SPLIT23") );
      setVector(tree, brancharray, var_SPLIT34_vec.at(i), std::string(jetString+"SPLIT34") );
      setVector(tree, brancharray, var_Dip12_vec.at(i), std::string(jetString+"Dip12") );
      setVector(tree, brancharray, var_Dip13_vec.at(i), std::string(jetString+"Dip13") );
      setVector(tree, brancharray, var_Dip23_vec.at(i), std::string(jetString+"Dip23") );
      setVector(tree, brancharray, var_DipExcl12_vec.at(i), std::string(jetString+"DipExcl12") );
      setVector(tree, brancharray, var_PlanarFlow_vec.at(i), std::string(jetString+"PlanarFlow") );
      setVector(tree, brancharray, var_Angularity_vec.at(i), std::string(jetString+"Angularity") );
      setVector(tree, brancharray, var_QW_vec.at(i), std::string(jetString+"QW") );
      setVector(tree, brancharray, var_PullMag_vec.at(i), std::string(jetString+"PullMag") );
      setVector(tree, brancharray, var_PullPhi_vec.at(i), std::string(jetString+"PullPhi") );
      setVector(tree, brancharray, var_Pull_C00_vec.at(i), std::string(jetString+"Pull_C00") );
      setVector(tree, brancharray, var_Pull_C01_vec.at(i), std::string(jetString+"Pull_C01") );
      setVector(tree, brancharray, var_Pull_C10_vec.at(i), std::string(jetString+"Pull_C10") );
      setVector(tree, brancharray, var_Pull_C11_vec.at(i), std::string(jetString+"Pull_C11") );

      setVector(tree, brancharray, var_ActiveArea_vec.at(i), std::string(jetString+"ActiveArea") );
      setVector(tree, brancharray, var_Aplanarity_vec.at(i), std::string(jetString+"Aplanarity") );
      setVector(tree, brancharray, var_Sphericity_vec.at(i), std::string(jetString+"Sphericity") );
      setVector(tree, brancharray, var_ThrustMaj_vec.at(i), std::string(jetString+"ThrustMaj") );
      setVector(tree, brancharray, var_ThrustMin_vec.at(i), std::string(jetString+"ThrustMin") );
      setVector(tree, brancharray, var_VoronoiArea_vec.at(i), std::string(jetString+"VoronoiArea") );

      setVector(tree, brancharray, var_TauWTA1_vec.at(i), std::string(jetString+"TauWTA1") );
      setVector(tree, brancharray, var_TauWTA2_vec.at(i), std::string(jetString+"TauWTA2") );
      setVector(tree, brancharray, var_TauWTA3_vec.at(i), std::string(jetString+"TauWTA3") );
      setVector(tree, brancharray, var_ZCUT12_vec.at(i), std::string(jetString+"ZCUT12") );
      setVector(tree, brancharray, var_ZCUT23_vec.at(i), std::string(jetString+"ZCUT23") );
      setVector(tree, brancharray, var_ZCUT34_vec.at(i), std::string(jetString+"ZCUT34") );
      
    } // end for loop over topo/truth/groom

  std::string jetString = returnJetType( samplePrefix, groomalgo, addLC, 2); //set to truth/ topo/ groomed

  if (subjetspre) // if the sample has the pre-calculated subjet variables
    {
      if(brancharray->FindObject(std::string(jetString+"config_massFraction").c_str() ) )
	tree->SetBranchAddress(std::string(jetString+"config_massFraction").c_str(),&var_massFraction_vec);
      else
	std::cout << "branch not found: " << std::string(jetString+"config_massFraction") << std::endl;
      if(brancharray->FindObject(std::string(jetString+"config_ktycut2").c_str()))
	tree->SetBranchAddress(std::string(jetString+"config_ktycut2").c_str(),&var_ktycut2_vec);
      else
	std::cout << "branch not found: " << std::string(jetString+"config_ktycut2") << std::endl;
    }

  if (!subjetscalc) // if the sample has the subjet branches available to create massdrop and yt
    return;

  std::string subjetType = returnSubJetType(samplePrefix, groomIdx, addLC);

  // set the branch status and address for the branch name that has the subjets for the grooming algorithm
  tree->SetBranchStatus(std::string(algorithms.subjetIndex[groomIdx]).c_str(),1);
  tree->SetBranchAddress(std::string(algorithms.subjetIndex[groomIdx]).c_str(),&subjet_index);

  if(!tree->GetListOfBranches()->FindObject(std::string(algorithms.subjetIndex[groomIdx]).c_str())) {
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
  tree->Branch("ptcone20", &var_ptcone20);//, "ptcone20/F");
  tree->Branch("etcone20", &var_etcone20);//, "etcone20/F");
  tree->Branch("charge", &var_charge);//, "charge/F");
  tree->Branch(string(jetString+"mass_llj").c_str(), &var_mllj, "mass_mllj/F");
  tree->Branch("mll", &var_mll, "mll/F");
  tree->Branch("ptll", &var_ptll, "ptll/F");
  tree->Branch("isElectronEvent", &var_isElectronEvent, "isElectronEvent/I");

} // addLeptonBranches

/*
 * Set up some vectors for easier access to the eta, phi and pt of the truth, topo and groomed jets.
 */
void setSelectionVectors()
{  
 
  jet_eta_truth = var_eta_vec[0];
  jet_phi_truth = var_phi_vec[0];
  jet_pt_truth = var_pt_vec[0];
  jet_m_truth = var_m_vec[0];
  
  jet_eta_topo = var_eta_vec[1];
  jet_phi_topo = var_phi_vec[1];
  jet_pt_topo = var_pt_vec[1];
  jet_emfrac_topo = var_emfrac_vec[1];
  
  jet_eta_groomed = var_eta_vec[2];
  jet_phi_groomed = var_phi_vec[2];
  jet_pt_groomed = var_pt_vec[2];
  jet_emfrac_groomed = var_emfrac_vec[2];
 
}// setSelectionVectors()

/*
 * Initialise all of the vectors used for reading in the input TTree.
 *
 */
void initVectors()
{

  // have the vectors for the above histograms so we can do the reading in stuff from the TTree
  jet_eta_truth = 0;
  jet_phi_truth = 0;
  jet_pt_truth = 0;
  jet_m_truth = 0;
  jet_eta_topo = 0;
  jet_phi_topo = 0;
  jet_pt_topo = 0;
  jet_emfrac_topo = 0;
  jet_eta_groomed = 0;
  jet_phi_groomed = 0;
  jet_pt_groomed = 0;
  jet_emfrac_groomed = 0;
  for (int i = 0; i < 3; i++) // for jetType::TRUTH/TOPO/GROOMED
    {

      var_E_vec[i] = 0;
      var_pt_vec[i] = 0;
      var_m_vec[i] = 0;
      var_eta_vec[i] = 0;
      var_phi_vec[i] = 0;
      var_emfrac_vec[i] = 0;
      var_constit_index[i] = 0;
      var_Tau1_vec[i] = 0;
      var_Tau2_vec[i] = 0;
      var_Tau3_vec[i] = 0;
      var_WIDTH_vec[i] = 0;
      var_SPLIT12_vec[i] = 0;
      var_SPLIT23_vec[i] = 0;
      var_SPLIT34_vec[i] = 0;
      var_Dip12_vec[i] = 0;
      var_Dip13_vec[i] = 0;
      var_Dip23_vec[i] = 0;
      var_DipExcl12_vec[i] = 0;
      var_PlanarFlow_vec[i] = 0;
      var_Angularity_vec[i] = 0;
      var_QW_vec[i] = 0;
      var_PullMag_vec[i] = 0;
      var_PullPhi_vec[i] = 0;
      var_Pull_C00_vec[i] = 0;
      var_Pull_C01_vec[i] = 0;
      var_Pull_C10_vec[i] = 0;
      var_Pull_C11_vec[i]= 0;//std::map<i, 0>;

      var_ActiveArea_vec[i] = 0;
      var_Aplanarity_vec[i] = 0;
      var_Sphericity_vec[i] = 0;
      var_ThrustMaj_vec[i] = 0;
      var_ThrustMin_vec[i] = 0;
      var_VoronoiArea_vec[i] = 0;
      
  
      var_TauWTA1_vec[i] = 0;
      var_TauWTA2_vec[i] = 0;
      var_TauWTA3_vec[i] = 0;
      
      
      var_ZCUT12_vec[i] = 0;
      var_ZCUT23_vec[i] = 0;
      var_ZCUT34_vec[i] = 0;
      
    }  
  var_YFilt_vec = 0;
  var_massFraction_vec = 0;
  var_ktycut2_vec = 0;

  subjet_index = 0;    
  var_subjets_E_vec = 0;
  var_subjets_pt_vec = 0;
  var_subjets_m_vec = 0;
  var_subjets_eta_vec = 0;
  var_subjets_phi_vec = 0;
  
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
    elsize = var_electronPt_vec->size();

  for (int x = 0; x < elsize; x++) // loop through all electrons
    {
      TLorentzVector e = TLorentzVector();
      e.SetPtEtaPhiM((*var_electronPt_vec)[x], (*var_electronEta_vec)[x], (*var_electronPhi_vec)[x], ELMASS);
      electrons.push_back(e);
      
    }

  int musize = 0;
  if (var_muonPt_vec != NULL)
    musize = var_muonPt_vec->size();
  TLorentzVector m = TLorentzVector();
  for (int x = 0; x < musize; x++)
    {
      
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
 * @param groomIdx The full name of the algorithm.
 */
void setOutputVariables( int jet_idx_truth, int jet_idx_topo, int jet_idx_groomed, int subjet_idx, std::string & groomalgo, std::string &  groomIdx)
{
  // set to 0 for now
  int jet_idx = 0;
  // set some of the variables that are not collections/ vectors
  mc_event_weight_out = mc_event_weight;
  mc_channel_number_out = mc_channel_number;
  runNumberOut = runNumberIn;
  nvtxOut = nvtxIn;
  avgIntpXingOut = avgIntpXingIn;

  std::string samplePrefix = ""; // AntiKt10/ CamKt12 for example
  bool addLC = false; // some algorithms have "LC" in their name
  samplePrefix = algorithms.AlgoPrefix[groomIdx];

  if (algorithms.AlgoType[groomIdx].find("recluster") == std::string::npos) // we're doing reclustering
    addLC = true; // just add the LC to the name

  // default weights
  var_k_factor = 1.0;  
  var_filter_eff = 1.0;
  var_xs = 1.0;
  // check if the weights read in from the file have entries for this channel number
  if (k_factors.count(long(mc_channel_number)) > 0 )
    {
      var_k_factor = k_factors[long(mc_channel_number)];
    }

  if (filt_eff.count(long(mc_channel_number)) > 0 )
    {
      var_filter_eff = filt_eff[long(mc_channel_number)];
    }

  if (xs.count(long(mc_channel_number)) > 0)
    {
      var_xs = xs[long(mc_channel_number)];
    }

  // yfilt only exists for groomed jets
  if (var_YFilt_vec != NULL && useBranch(string(returnJetType(samplePrefix, groomalgo, addLC,2)+"YFilt")))
    var_YFilt=(*var_YFilt_vec)[jet_idx];
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
      
      string jetString = returnJetType( samplePrefix, groomalgo, addLC, x); //set to truth/ topo/ groomed;
      
      // before setting the output variable we check to see if the variable we are setting up is
      // not null and if we are wanting to write out that information anyway.
      if (var_E_vec[x] != NULL && useBranch(string(jetString+"E")))
	var_E[x]=(*var_E_vec[x])[jet_idx];
      if (var_pt_vec[x] != NULL && useBranch(string(jetString+"pt")))
	var_pt[x]=(*var_pt_vec[x])[jet_idx];
      if (var_m_vec[x] != NULL && useBranch(string(jetString+"m")))
	var_m[x]=(*var_m_vec[x])[jet_idx];
      if (var_eta_vec[x] != NULL && useBranch(string(jetString+"eta")))
	var_eta[x]=(*var_eta_vec[x])[jet_idx];
      if (var_phi_vec[x] != NULL && useBranch(string(jetString+"phi")))
	var_phi[x]=(*var_phi_vec[x])[jet_idx];
      if (x!=0 && xAODemfrac && useBranch(string(jetString +"emfrac")))
	var_emfrac[x]=(*var_emfrac_vec[x])[jet_idx];
      if (var_Tau1_vec[x] != NULL && useBranch(string(jetString+"Tau1")))
	var_Tau1[x]=(*var_Tau1_vec[x])[jet_idx];
      if (var_Tau2_vec[x] != NULL && useBranch(string(jetString+"Tau2")))
	var_Tau2[x]=(*var_Tau2_vec[x])[jet_idx];
      if (var_Tau3_vec[x] != NULL && useBranch(string(jetString+"Tau3")))
	var_Tau3[x]=(*var_Tau3_vec[x])[jet_idx];
      if (var_WIDTH_vec[x] != NULL && useBranch(string(jetString+"WIDTH")))
	var_WIDTH[x]=(*var_WIDTH_vec[x])[jet_idx];
      if (var_SPLIT12_vec[x] != NULL && useBranch(string(jetString+"SPLIT12")))
	var_SPLIT12[x]=(*var_SPLIT12_vec[x])[jet_idx];
      if (var_SPLIT23_vec[x] != NULL && useBranch(string(jetString+"SPLIT23")))
	var_SPLIT23[x]=(*var_SPLIT23_vec[x])[jet_idx];
      if (var_SPLIT34_vec[x] != NULL && useBranch(string(jetString+"SPLIT34")))
	var_SPLIT34[x]=(*var_SPLIT34_vec[x])[jet_idx];
      if (var_Dip12_vec[x] != NULL && useBranch(string(jetString+"Dip12")))
	var_Dip12[x]=(*var_Dip12_vec[x])[jet_idx];
      if (var_Dip13_vec[x] != NULL && useBranch(string(jetString+"Dip13")))
	var_Dip13[x]=(*var_Dip13_vec[x])[jet_idx];
      if (var_Dip23_vec[x] != NULL && useBranch(string(jetString+"Dip23")))
	var_Dip23[x]=(*var_Dip23_vec[x])[jet_idx];
      if (var_DipExcl12_vec[x] != NULL && useBranch(string(jetString+"DipExcl12")))
	var_DipExcl12[x]=(*var_DipExcl12_vec[x])[jet_idx];
      if (var_PlanarFlow_vec[x] != NULL && useBranch(string(jetString+"PlanarFlow")))
	var_PlanarFlow[x]=(*var_PlanarFlow_vec[x])[jet_idx];
      if (var_Angularity_vec[x] != NULL && useBranch(string(jetString+"Angularity")))
	var_Angularity[x]=(*var_Angularity_vec[x])[jet_idx];
      if (var_QW_vec[x] != NULL && useBranch(string(jetString+"QW")))
	var_QW[x]=(*var_QW_vec[x])[jet_idx];
      if (var_PullMag_vec[x] != NULL && useBranch(string(jetString+"PullMag")))
	var_PullMag[x]=(*var_PullMag_vec[x])[jet_idx];
      if (var_PullPhi_vec[x] != NULL && useBranch(string(jetString+"PullPhi")))
	var_PullPhi[x]=(*var_PullPhi_vec[x])[jet_idx];
      if (var_Pull_C00_vec[x] != NULL && useBranch(string(jetString+"Pull_C00")))
	var_Pull_C00[x]=(*var_Pull_C00_vec[x])[jet_idx];
      if (var_Pull_C01_vec[x] != NULL && useBranch(string(jetString+"Pull_C01")))
	var_Pull_C01[x]=(*var_Pull_C01_vec[x])[jet_idx];
      if (var_Pull_C10_vec[x] != NULL && useBranch(string(jetString+"Pull_C10")))
	var_Pull_C10[x]=(*var_Pull_C10_vec[x])[jet_idx];
      if (var_Pull_C11_vec[x] != NULL && useBranch(string(jetString+"Pull_C11")))
	var_Pull_C11[x]=(*var_Pull_C11_vec[x])[jet_idx];

      if (var_ActiveArea_vec[x] != NULL && useBranch(string(jetString+"ActiveArea")))
	var_ActiveArea[x] = (*var_ActiveArea_vec[x])[jet_idx];
      if (var_Aplanarity_vec[x] != NULL && useBranch(string(jetString+"Aplanarity")))
	var_Aplanarity[x] = (*var_Aplanarity_vec[x])[jet_idx];
      if (var_Sphericity_vec[x] != NULL && useBranch(string(jetString+"Sphericity")))
	var_Sphericity[x] = (*var_Sphericity_vec[x])[jet_idx];
      if (var_ThrustMaj_vec[x] != NULL && useBranch(string(jetString+"ThrustMaj")))
	var_ThrustMaj[x] = (*var_ThrustMaj_vec[x])[jet_idx];
      if (var_ThrustMin_vec[x] != NULL && useBranch(string(jetString+"ThrustMin")))
	var_ThrustMin[x] = (*var_ThrustMin_vec[x])[jet_idx];
      if (var_VoronoiArea_vec[x] != NULL && useBranch(string(jetString+"VoronoiArea")))
	var_VoronoiArea[x] = (*var_VoronoiArea_vec[x])[jet_idx];
      
      // tau21 and tauwta21 are set in the main loop, not here, because we have to calculate them
      if (useBranch(string(jetString+"TauWTA1")))
	var_TauWTA1[x]=(*var_TauWTA1_vec[x])[jet_idx];
      if (useBranch(string(jetString+"TauWTA2")))
	var_TauWTA2[x]=(*var_TauWTA2_vec[x])[jet_idx];
      if (useBranch(string(jetString+"TauWTA3")))
	var_TauWTA3[x]=(*var_TauWTA3_vec[x])[jet_idx];
      if (useBranch(string(jetString+"ZCUT12")))
	var_ZCUT12[x]=(*var_ZCUT12_vec[x])[jet_idx];
      if (useBranch(string(jetString+"ZCUT23")))
	var_ZCUT23[x]=(*var_ZCUT23_vec[x])[jet_idx];
      if (useBranch(string(jetString+"ZCUT34")))
	var_ZCUT34[x]=(*var_ZCUT34_vec[x])[jet_idx];
	
    } // end for loop

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
  var_Tau3.clear();
  var_WIDTH.clear();
  var_SPLIT12.clear();
  var_SPLIT23.clear();
  var_SPLIT34.clear();
  var_Dip12.clear();
  var_Dip13.clear();
  var_Dip23.clear();
  var_DipExcl12.clear();
  var_PlanarFlow.clear();
  var_Angularity.clear();
  var_QW.clear();
  var_PullMag.clear();
  var_PullPhi.clear();
  var_Pull_C00.clear();
  var_Pull_C01.clear();
  var_Pull_C10.clear();
  var_Pull_C11.clear();

  var_ActiveArea.clear();
  var_Aplanarity.clear();
  var_Sphericity.clear();
  var_ThrustMaj.clear();
  var_ThrustMin.clear();
  var_VoronoiArea.clear();

  var_Tau21.clear();
  var_TauWTA2TauWTA1.clear();
  var_TauWTA1.clear();
  var_TauWTA2.clear();
  var_TauWTA3.clear();
  var_ZCUT12.clear();
  var_ZCUT23.clear();
  var_ZCUT34.clear();



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
      var_Tau3.push_back(-999);
      var_WIDTH.push_back(-999);
      var_SPLIT12.push_back(-999);
      var_SPLIT23.push_back(-999);
      var_SPLIT34.push_back(-999);
      var_Dip12.push_back(-999);
      var_Dip13.push_back(-999);
      var_Dip23.push_back(-999);
      var_DipExcl12.push_back(-999);
      var_PlanarFlow.push_back(-999);
      var_Angularity.push_back(-999);
      var_QW.push_back(-999);
      var_PullMag.push_back(-999);
      var_PullPhi.push_back(-999);
      var_Pull_C00.push_back(-999);
      var_Pull_C01.push_back(-999);
      var_Pull_C10.push_back(-999);
      var_Pull_C11.push_back(-999);

      var_ActiveArea.push_back(-999);
      var_VoronoiArea.push_back(-999);
      var_Aplanarity.push_back(-999);
      var_Sphericity.push_back(-999);
      var_ThrustMaj.push_back(-999);
      var_ThrustMin.push_back(-999);

      var_Tau21.push_back(-999);
      var_TauWTA2TauWTA1.push_back(-999);
      var_TauWTA1.push_back(-999);
      var_TauWTA2.push_back(-999);
      var_TauWTA3.push_back(-999);
      var_ZCUT12.push_back(-999);
      var_ZCUT23.push_back(-999);
      var_ZCUT34.push_back(-999);

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

  tree->Branch("mc_event_weight",&mc_event_weight_out,"mc_event_weight/F");
  tree->Branch("mc_channel_number", &mc_channel_number_out,"mc_channel_number/I");
  tree->Branch("vxp_n", &nvtxOut, "vxp_n/I");
  tree->Branch("averageIntPerXing",&avgIntpXingOut,"averageIntPerXing/F");

  for (int i = 0; i < jetType::MAX; i++) // truth, topo, groomed
    {
     
      std::string jetString = returnJetType(samplePrefix, groomalgo, addLC,i); //set to truth/ topo/ groomed

      tree->Branch(std::string(jetString+"E").c_str(),&var_E.at(i),std::string(jetString+"E/F").c_str());
      tree->Branch(std::string(jetString+"pt").c_str(),&var_pt.at(i),std::string(jetString+"pt/F").c_str());
      tree->Branch(std::string(jetString+"m").c_str(),&var_m.at(i),std::string(jetString+"m/F").c_str());
      tree->Branch(std::string(jetString+"eta").c_str(),&var_eta.at(i),std::string(jetString+"eta/F").c_str());
      tree->Branch(std::string(jetString+"phi").c_str(),&var_phi.at(i),std::string(jetString+"phi/F").c_str());
      if (i != jetType::TRUTH) // emfrac doesn't exist for truth jets
	tree->Branch(std::string(jetString+"emfrac").c_str(),&var_emfrac.at(i),std::string(jetString+"emfrac"+"/F").c_str());
      tree->Branch(std::string(jetString+"Tau1").c_str(),&var_Tau1.at(i),std::string(jetString+"Tau1/F").c_str());
      tree->Branch(std::string(jetString+"Tau2").c_str(),&var_Tau2.at(i),std::string(jetString+"Tau2/F").c_str());
      tree->Branch(std::string(jetString+"Tau3").c_str(),&var_Tau3.at(i),std::string(jetString+"Tau3/F").c_str());
      tree->Branch(std::string(jetString+"WIDTH").c_str(),&var_WIDTH.at(i),std::string(jetString+"WIDTH/F").c_str());
      tree->Branch(std::string(jetString+"SPLIT12").c_str(),&var_SPLIT12.at(i),std::string(jetString+"SPLIT12/F").c_str());
      tree->Branch(std::string(jetString+"SPLIT23").c_str(),&var_SPLIT23.at(i),std::string(jetString+"SPLIT23/F").c_str());
      tree->Branch(std::string(jetString+"SPLIT34").c_str(),&var_SPLIT34.at(i),std::string(jetString+"SPLIT34/F").c_str());
      tree->Branch(std::string(jetString+"Dip12").c_str(),&var_Dip12.at(i),std::string(jetString+"Dip12/F").c_str());
      tree->Branch(std::string(jetString+"Dip13").c_str(),&var_Dip13.at(i),std::string(jetString+"Dip13/F").c_str());
      tree->Branch(std::string(jetString+"Dip23").c_str(),&var_Dip23.at(i),std::string(jetString+"Dip23/F").c_str());
      tree->Branch(std::string(jetString+"DipExcl12").c_str(),&var_DipExcl12.at(i),std::string(jetString+"DipExcl12/F").c_str());
      tree->Branch(std::string(jetString+"PlanarFlow").c_str(),&var_PlanarFlow.at(i),std::string(jetString+"PlanarFlow/F").c_str());
      tree->Branch(std::string(jetString+"Angularity").c_str(),&var_Angularity.at(i),std::string(jetString+"Angularity/F").c_str());
      tree->Branch(std::string(jetString+"QW").c_str(),&var_QW.at(i),std::string(jetString+"QW/F").c_str());
      tree->Branch(std::string(jetString+"PullMag").c_str(),&var_PullMag.at(i),std::string(jetString+"PullMag/F").c_str());
      tree->Branch(std::string(jetString+"PullPhi").c_str(),&var_PullPhi.at(i),std::string(jetString+"PullPhi/F").c_str());
      tree->Branch(std::string(jetString+"Pull_C00").c_str(),&var_Pull_C00.at(i),std::string(jetString+"Pull_C00/F").c_str());
      tree->Branch(std::string(jetString+"Pull_C01").c_str(),&var_Pull_C01.at(i),std::string(jetString+"Pull_C01/F").c_str());
      tree->Branch(std::string(jetString+"Pull_C10").c_str(),&var_Pull_C10.at(i),std::string(jetString+"Pull_C10/F").c_str());
      tree->Branch(std::string(jetString+"Pull_C11").c_str(),&var_Pull_C11.at(i),std::string(jetString+"Pull_C11/F").c_str());
      tree->Branch(std::string(jetString+"YFilt").c_str(),&var_YFilt,std::string(jetString+"YFilt/F").c_str());

tree->Branch(std::string(jetString+"ActiveArea").c_str(),&var_ActiveArea,std::string(jetString+"ActiveArea/F").c_str());
tree->Branch(std::string(jetString+"VoronoiArea").c_str(),&var_VoronoiArea,std::string(jetString+"VoronoiArea/F").c_str());
tree->Branch(std::string(jetString+"Aplanarity").c_str(),&var_Aplanarity,std::string(jetString+"Aplanarity/F").c_str());
tree->Branch(std::string(jetString+"Sphericity").c_str(),&var_Sphericity,std::string(jetString+"Sphericity/F").c_str());
tree->Branch(std::string(jetString+"ThrustMaj").c_str(),&var_ThrustMaj,std::string(jetString+"ThrustMaj/F").c_str());
tree->Branch(std::string(jetString+"ThrustMin").c_str(),&var_ThrustMin,std::string(jetString+"ThrustMin/F").c_str());
      

 tree->Branch(std::string(jetString+"TauWTA1").c_str(),&var_TauWTA1.at(i),std::string(jetString+"TauWTA1/F").c_str());
 tree->Branch(std::string(jetString+"TauWTA2").c_str(),&var_TauWTA2.at(i),std::string(jetString+"TauWTA2/F").c_str());
 tree->Branch(std::string(jetString+"TauWTA3").c_str(),&var_TauWTA3.at(i),std::string(jetString+"TauWTA3/F").c_str());
 tree->Branch(std::string(jetString+"ZCUT12").c_str(),&var_ZCUT12.at(i),std::string(jetString+"ZCUT12/F").c_str());
 tree->Branch(std::string(jetString+"ZCUT23").c_str(),&var_ZCUT23.at(i),std::string(jetString+"ZCUT23/F").c_str());
 tree->Branch(std::string(jetString+"ZCUT34").c_str(),&var_ZCUT34.at(i),std::string(jetString+"ZCUT34/F").c_str());
	  
 tree->Branch(std::string(returnJetType(samplePrefix, groomalgo, addLC,0)+"TauWTA2/TauWTA1").c_str(),&var_TauWTA2TauWTA1.at(0),std::string(jetString+"TauWTA2TauWTA1/F").c_str());
 tree->Branch(std::string(returnJetType(samplePrefix, groomalgo, addLC,1)+"TauWTA2/TauWTA1").c_str(),&var_TauWTA2TauWTA1.at(1),std::string(jetString+"TauWTA2TauWTA1/F").c_str());
 tree->Branch(std::string(jetString+"TauWTA2/TauWTA1").c_str(),&var_TauWTA2TauWTA1.at(2),std::string(jetString+"TauWTA2TauWTA1/F").c_str());  
	  
      // add a calculated variable Tau2/Tau1
      tree->Branch(std::string(returnJetType(samplePrefix, groomalgo, addLC,0)+"Tau21").c_str(),&var_Tau21.at(0),std::string(jetString+"Tau21/F").c_str());
      tree->Branch(std::string(returnJetType(samplePrefix, groomalgo, addLC,1)+"Tau21").c_str(),&var_Tau21.at(1),std::string(jetString+"Tau21/F").c_str());
      tree->Branch(std::string(jetString+"Tau21").c_str(),&var_Tau21.at(2),std::string(jetString+"Tau21/F").c_str());  

    }
  std::string jetString = returnJetType(samplePrefix, groomalgo, addLC,jetType::GROOMED); 
  // add the lepton branches if running hvtllqq analysis
  if (hvtllqq)
    addLeptonBranches(jetString, tree);
} // setOutputBranches


/*
 * Normalise all of the histograms for Wprime and QCD samples.
 */
void scaleHists()
{

    Wprime_Lead_CA12_scaled_pt->Scale(1.0/Wprime_Lead_CA12_scaled_pt->Integral());

    for (int i = 0; i < 3; i ++)
      {
	qcd_Lead_CA12_pt[i]->Scale(1.0/qcd_Lead_CA12_pt[i]->Integral());  
	Wprime_Lead_CA12_pt[i]->Scale(1.0/Wprime_Lead_CA12_pt[i]->Integral());	
	for (int j=0; j<nPtBins; j++){
	  qcd_Lead_CA12_mass[i][j]->Scale(1.0/qcd_Lead_CA12_mass[i][j]->Integral());
	  Wprime_Lead_CA12_mass[i][j]->Scale(1.0/Wprime_Lead_CA12_mass[i][j]->Integral());
	}
	
	for (int j=0; j<nFineBins; j++){	  
	  qcd_finePtBin_mass[i][j]->Scale(1.0/qcd_finePtBin_mass[i][j]->Integral());
	  Wprime_finePtBin_mass[i][j]->Scale(1.0/Wprime_finePtBin_mass[i][j]->Integral());
	}	
      }

} // scaleHists

/*
 * Calculate the most probable value for the WPrime sample.
 */
void getMPV()
{
    for (int j=0; j<nPtBins; j++){
      myMPV[j]=mpv(Wprime_Lead_CA12_mass[histType::GROOMEDJET][j]);
    }
    
    for (int j=0; j<nFineBins; j++){
      myMPV_finePt[j]=mpv(Wprime_finePtBin_mass[histType::GROOMEDJET][j]);
    }
} //getMPV


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
  ifstream wf ("weightings_input.csv");
  string line;
  std::cout << "reading weights" << std::endl;
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
  int bins = bkg->GetNbinsX();
  ofstream out(fname);
  double weight;

  for (int b = 1; b <= bins; b++) // loop through all bins
    {
      weight = bkg->GetBinContent(b);
      if (sig->GetBinContent(b) != 0)
	weight/=sig->GetBinContent(b);
      else
	{
	  std::cout << "oh no, we've been rumbled!  We have no signal in this bin and now we have to reweight by over 9000! But actually we're just going to make it 0." << std::endl;
	  weight = 0;
	}
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
  BOOST_FOREACH( ptree::value_type & v, pt.get_child("config.algorithms"))
    {
      std::string name =  v.second.get<std::string>("Algorithm","");
      if (name != "")
	{
	  std::string groomalgo = v.second.get<std::string>("GroomingAlgorithm","");
	  if (groomalgo == "")
	    {
	      groomalgo = name.substr(name.find("Topo"),name.length()-1); // shortened version of name, - AntiKt10 for eg.
	    }
	  AlgoNames[name] = groomalgo;
	  std::string jetalgo = v.second.get<std::string>("JetAlgorithm","");
	  if (jetalgo == "")
	    {
	      jetalgo = name.substr(0,name.find("Topo"));
	      if (jetalgo.find("LC") != std::string::npos)
		jetalgo.erase(jetalgo.end()-2,jetalgo.end()); // erase LC from the name
	    }
	  AlgoPrefix[name] = jetalgo;
	  AlgoList[name] = v.second.get<std::string>("Abbreviation","");
	  AlgoListN[name] = v.second.get<std::string>("PlotLabel","");
	  AlgoType[name] = v.second.get<std::string>("Type","NONE");
	  
	  subjetMap[name] = v.second.get<std::string>("SubjetGroomingAlgorithm");
	  subjetIndex[name] = v.second.get<std::string>("SubjetIndexBranch");
	  binLabel[name] = v.second.get<std::string>("BinLabel","");
	}
    };

} // algorithms::load()

