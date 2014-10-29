///
///F. A. Dias (flavia.dias@cern.ch)
///Analysis code for Jet Substructure - Munich workshop
///September 2014
///
///v0.1 - August 14th, 2014
///v0.2 - August 25th, 2014
///v0.3 - August 26th, 2014
///v0.4 - August 30th, 2014
///
/// usage: bin/PlotsMunichWorkshop QCD_SplitFiltered.root Wprime_SplitFiltered.root QCD_Trim.root Wprime_Trim.root QCD_Prune.root Wprime_Prune.root QCD_Reclust.root Wprime_Reclust.root
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
  //subjets = true; // should be set from command line argument TODO
  gROOT->ProcessLine("#include <vector>");
  gInterpreter->GenerateDictionary("std::vector<std::vector<int> >", "vector");
  gInterpreter->GenerateDictionary("vector<vector<double> >", "vector");
  
  AtlasStyle();
  
  bool makePtPlotsFlag = false;
  bool extendedVars = false;
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
      ("extendedvars", po::value<bool>(&extendedVars)->default_value(false),"add extra variables - TauWTA1/2 and ZCUT12")
      ("mass-window", po::value<bool>(&applyMassWindowFlag)->default_value(false),"apply mass window cuts")
      ("make-plots", po::value<bool>(&makePlotsFlag)->default_value(false),"create plots")
      ("make-ptplots", po::value<bool>(&makePtPlotsFlag)->default_value(false),"create pT plots")
      ("scale-hists", po::value<bool>(&scaleHistsFlag)->default_value(false),"scale mass/ pt histograms")
      ("mpv", po::value<bool>(&getMPVFlag)->default_value(false),"Find the MPV")
      ("fileid", po::value<string>(&fileid)->default_value(""),"Add an identifier to the output folder/ file names")
      ("bkg-frac", po::value<bool>(&checkBkgFrac)->default_value(false),"Check the background fraction in the signal")
      ("tree-name", po::value<string>(&treeName)->default_value("physics"),"Name of tree to be read in from input file")
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

  algorithms.load("algorithms.xml");

  bool massHistos = applyMassWindowFlag;

  fileid_global = fileid;

  vector<string> inputBkgFiles = vm["background-file"].as<vector<string> > ();
  vector<string> inputSigFiles = vm["signal-file"].as<vector<string> > ();

  std::string alg_in = vm["algorithm"].as<string>();
  std::string alg_in_orig = vm["algorithm"].as<string>();

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

  std::transform(alg_in.begin(), alg_in.end(), alg_in.begin(), ::tolower);
  vector<int> subset;

  int nArg = 1; // TODO: this will need to be changed when moving from a system where only one algorithm is run at once

  // right now the idea of filteridx, Xidx is not fully complete... The code as it stands runs on one type of algorithm at once, but 
  // I would like it to run on a number of types.  This is why this is here.  It used to do this in fact, before switching stuff around

  std::cout << alg_in_orig << std::endl;
  //for (std::map<std::string,std::string>::iterator it = algorithms.AlgoNames.begin(); it!= algorithms.AlgoNames.end(); it++){
  //std::cout << (*it) << std::endl;
  //}

  if (algorithms.AlgoNames.find(alg_in_orig) == algorithms.AlgoNames.end())
    {
      std::cout << "Algorithm not defined in config file!" << std::endl;
      return 1;
    }
  
  
  
  //}
  readWeights();
  defineStrings(pTbins, finePtBins);
  createHistos(alg_in_orig);
  
  // TODO: option to run all algorithms of one type????

  runAlgorithm(inputTChain[sampleType::BACKGROUND], inputTChain[sampleType::SIGNAL], algorithms.AlgoNames[alg_in_orig], alg_in_orig, massHistos);
  

  //Plot things
  //Normalise all histograms of interest
  // this is not the right way to do the normalisation according to this - https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/BoostedBosonTaggingD3PD#MC_only_studies

  // NEvents = Sum[mcevt_weight[0][0]*PileupWeight
  // Note, for NEvents for SherpaW + jets we need to call getNormSherpaW() because the weighting is slightly different
  // Normalization = CrossSection(at NLO) * Luminosity / NEvents


  if (scaleHistsFlag)
    {
      scaleHists();
    } // end scaleHists if
  


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
    

  //3. Check background fraction in this window  
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

  makeMassWindowFile(applyMassWindowFlag, extendedVars, alg_in_orig);

  return 0;
}



float DeltaR(float eta1,float phi1,float eta2,float phi2){ 
  // we could actually inline this function
  float deltaPhi = TMath::Abs(phi1-phi2);
  float deltaEta = eta1-eta2;
  if(deltaPhi > TMath::Pi())
    deltaPhi = TMath::TwoPi() - deltaPhi;
  return TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
}


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
  //cout<<"MPV: "<<mostprob<<endl;                    
                                  
  return mostprob;
}


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
	if ((*qcd_CA12_topo_emfrac)[n]<0.99 && (*qcd_CA12_topo_pt)[n]/1000>15.0){
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
      for (int i=0; i<(*qcd_CA12_topo_pt).size(); i++){
	if (!hasTopoJet && (*qcd_CA12_topo_emfrac)[i]<0.99 && fabs((*qcd_CA12_topo_eta)[i])<1.2) {
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

      if (chosenLeadTruthJetIndex>=0 && chosenLeadGroomedIndex<0 && DeltaR((*qcd_CA12_truth_eta)[chosenLeadTruthJetIndex],(*qcd_CA12_truth_phi)[chosenLeadTruthJetIndex],(*qcd_CA12_groomed_eta)[i],(*qcd_CA12_groomed_phi)[i])<0.9 && (*qcd_CA12_groomed_emfrac)[i]<0.99 && fabs((*qcd_CA12_groomed_eta)[i])<1.2){
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
	TLorentzVector tempJet;
	if ((*Wp_CA12_topo_emfrac)[n]<0.99){
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
	if (!hasTopoJet && (*Wp_CA12_topo_emfrac)[i]<0.99 && fabs((*Wp_CA12_topo_eta)[i])<1.2) {
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

      if (chosenLeadTruthJetIndex>=0 && chosenLeadGroomedIndex<0 && DeltaR((*Wp_CA12_truth_eta)[chosenLeadTruthJetIndex],(*Wp_CA12_truth_phi)[chosenLeadTruthJetIndex],(*Wp_CA12_groomed_eta)[i],(*Wp_CA12_groomed_phi)[i])<0.9 && (*Wp_CA12_groomed_emfrac)[i]<0.99 && fabs((*Wp_CA12_groomed_eta)[i])<1.2){
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
    inputTree->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LCTopo_emfrac").c_str(), &qcd_CA12_topo_emfrac);
    
    inputTree->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LC"+groomAlgo+"_pt").c_str(), &qcd_CA12_groomed_pt);
    inputTree->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LC"+groomAlgo+"_eta").c_str(), &qcd_CA12_groomed_eta);
    inputTree->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LC"+groomAlgo+"_phi").c_str(), &qcd_CA12_groomed_phi);
    inputTree->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LC"+groomAlgo+"_m").c_str(), &qcd_CA12_groomed_mass);
    inputTree->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LC"+groomAlgo+"_E").c_str(), &qcd_CA12_groomed_E);
    inputTree->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LC"+groomAlgo+"_emfrac").c_str(), &qcd_CA12_groomed_emfrac);
    
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
    inputTree1->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LCTopo_emfrac").c_str(), &Wp_CA12_topo_emfrac);
    
    inputTree1->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LC"+groomAlgo+"_pt").c_str(), &Wp_CA12_groomed_pt);
    inputTree1->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LC"+groomAlgo+"_eta").c_str(), &Wp_CA12_groomed_eta);
    inputTree1->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LC"+groomAlgo+"_phi").c_str(), &Wp_CA12_groomed_phi);
    inputTree1->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LC"+groomAlgo+"_m").c_str(), &Wp_CA12_groomed_mass);
    inputTree1->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LC"+groomAlgo+"_E").c_str(), &Wp_CA12_groomed_E);
    inputTree1->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LC"+groomAlgo+"_emfrac").c_str(), &Wp_CA12_groomed_emfrac);


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
    inputTree->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+""+groomAlgo+"_emfrac").c_str(), &qcd_CA12_groomed_emfrac);
    
    
    
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
    inputTree1->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+"LCTopo_emfrac").c_str(), &Wp_CA12_topo_emfrac);
    
    inputTree1->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+""+groomAlgo+"_pt").c_str(), &Wp_CA12_groomed_pt);
    inputTree1->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+""+groomAlgo+"_eta").c_str(), &Wp_CA12_groomed_eta);
    inputTree1->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+""+groomAlgo+"_phi").c_str(), &Wp_CA12_groomed_phi);
    inputTree1->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+""+groomAlgo+"_m").c_str(), &Wp_CA12_groomed_mass);
    inputTree1->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+""+groomAlgo+"_E").c_str(), &Wp_CA12_groomed_E);
    inputTree1->SetBranchAddress(std::string("jet_"+algorithms.AlgoPrefix[groomAlgoIndex]+""+groomAlgo+"_emfrac").c_str(), &Wp_CA12_groomed_emfrac);
    

  }


 
}


void makeROC(int type, TH1F *&S, TH1F *&B, TGraph &curve, TString name, bool draw){//, TGraph &curve_up, TGraph &curve_do){
  //Usage:
  //TH1D *S = new TH1D();
  //TH1D *B = new TH1D();
  //...Fill S and B with signal (S) and background (B)
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

}



void makePlots(std::string & algo){
  ////////PLOT THINGS -> This should go into it's own class or at the very least it's own method, TODO

  // CANVAS FOR THE MASS PLOTS
  
  for (int i=0; i<3; i++){//-2; i++){
  //for (int ii=0; ii<nAlgos; ii++){//-2; i++){
    //std::cout << "ii: " << ii << std::endl;
    //int i = algoMap[ii];
    //std::cout << "i: " << i << std::endl;
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


void makeMassWindowFile(bool applyMassWindow,bool extendedVars, std::string & algorithm)
{
  std::cout << "starting make mass window output files " << std::endl;
  double mass_max = 0.0;
  double mass_min = 0.0;

  bool signal = false;

  // need branches for specific algorithms because otherwise we end up with a million branches that we don't need
  vector<std::string> branches;

  // loop through different algorithms
  //for (int ii=0; ii<nAlgos; ii++){

  std::string i = algorithm;//algoMap[ii];
    std::cout << "Doing mass window plots for " << algorithms.AlgoNames[i] << std::endl;

    std::string groomAlgoIndex = i;
    std::string prefix = "";

    branches = getListOfJetBranches(algorithms.AlgoNames[i]);
    // loop through the different pt bins. j == 0 is the inclusive one
    for (int j=0; j<1; j++){//nPtBins; j++){

      signal = false;
      std::stringstream ss;
      ss << algorithms.AlgoNames[i] << "_" << i << "_" << pTbins[j];
      // loop through background and signal
      for (int k = 0; k < 2; k++)
	{
	  signal = k == 0 ? false : true;
	  // Initialise all the vectors to.. something 
	  initVectors(extendedVars);
	  
	  int tchainIdx = signal ? sampleType::SIGNAL : sampleType::BACKGROUND;
	  std::stringstream ss_fname; // store the name of the output file and include the i and j indices!
	  std::string bkg = signal ? "sig": "bkg";
	  ss_fname << ss.str() << "_" << bkg << ".root";
	  boost::filesystem::path dir(std::string(algorithms.AlgoNames[i]+fileid_global));
	  boost::filesystem::create_directory(dir);

	  inputTChain[tchainIdx]->GetEntries();
	  //TTree * intree = (TTree*)inputTChain[tchainIdx]->GetTree();


	  // I can't remember if I need intree anymore... it was for debugging something
	  try
	    {
	      inputTChain[tchainIdx]->SetBranchStatus("*",0);
	    }
	  catch (exception &e)
	    {
	      std::cout << "There is a missing algorithm in the tree" << std::endl;
	      return;
	    }
	  //intree->SetBranchStatus("*",0);
	  // turn on the branches we're interested in

	  for (std::vector<string>::iterator it = branches.begin(); it != branches.end(); it++)
	    {
	      //inputTChain[tchainIdx]->SetBranchStatus((*it).c_str(),1);
	      //intree->SetBranchStatus((*it).c_str(),1);
	      if(!inputTChain[tchainIdx]->GetListOfBranches()->FindObject((*it).c_str())) {
		std::cout << "could not find branch: " << (*it) << std::endl;
	      }
	      else
		{
		  inputTChain[tchainIdx]->SetBranchStatus((*it).c_str(),1);
		  //intree->SetBranchStatus((*it).c_str(),1);
		}
	      
	    }
	  setJetsBranches(inputTChain[tchainIdx], algorithms.AlgoNames[i], signal, i, extendedVars); //set all of the branches for the output tree for the jets	  

	  std::cout << "set jetsbranches" << std::endl;
	  //mass_max = TopEdgeMassWindow[i][j];
	  mass_max = TopEdgeMassWindow[j];
	  //mass_min = BottomEdgeMassWindow[i][j];
	  mass_min = BottomEdgeMassWindow[j];

	  TFile * outfile = new TFile(std::string(algorithms.AlgoNames[i]+fileid_global+"/"+ss_fname.str()).c_str(),"RECREATE");	  
	  TTree * outTree = new TTree(treeName.c_str(),treeName.c_str());

	  resetOutputVariables();
	  setOutputBranches(outTree, algorithms.AlgoNames[i], i, extendedVars);
	  //addJetsBranches(outTree, AlgoNames[i], signal, i);

	  if (subjetscalc || subjetspre)
	    addSubJets(outTree, algorithms.AlgoNames[i], signal, i);
	  addJetsBranches(outTree, algorithms.AlgoNames[i], signal);

	  long entries = (long)inputTChain[tchainIdx]->GetEntries();

	  // Setting up this with a high limit of 3.5 TeV so we don't miss anything.  Lots of bins - 200, so we can
	  // do a lot of tuning of the scale factor regions later on!
	  pt_reweight = new TH1F(std::string("pt_reweight"+bkg).c_str(),std::string("pt_reweight_"+bkg).c_str(), 200, 0, 3500);
	  
	  double mass = 0;

	  NEvents = entries;
	  NEvents_weighted.clear();

	  std::cout << "total entries: " << entries << std::endl;
	  for (long n = 0; n < entries; n++)
	    {
	      inputTChain[tchainIdx]->GetEntry(n);

	      setSelectionVectors(signal, algorithms.AlgoNames[i]); // again reading from tree, but just setting jet_pt_x whatever
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


	      if ((*jet_pt_truth)[0] / 1000.0 < 100)// || (*jet_m_truth)[0]/1000.0 < 100) 
		continue;
	      
	      if (algorithms.AlgoType[groomAlgoIndex].find("truthmatch") == std::string::npos) // check the pt is in the correct bin
		{
		  if ((*jet_pt_truth)[0]/1000.0 < ptrange[j].first || (*jet_pt_truth)[0]/1000.0 > ptrange[j].second)
		    continue;
		}
	      
	      if (algorithms.AlgoType[groomAlgoIndex].find("truthmatch") != std::string::npos)
		{
		  //check which is my CA12 ungroomed reco jet with EMfrac<0.99 and |eta|<1.2
		  //loop over all topo jets and get the leading with cuts
		  bool hasTopoJet=false;
		
		  for (int jet_i=0; jet_i<(*jet_pt_topo).size(); jet_i++)
		    {
		      if (!hasTopoJet && (*jet_emfrac_topo)[jet_i]<0.99 && fabs((*jet_eta_topo)[jet_i])<1.2) 
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

	      if (chosenLeadTruthJetIndex < 0)
		continue;
	      //Now I have which events to make my pt reweight with, and to match to, etc
	      int chosenLeadGroomedIndex=-99;

	      for (int jet_i=0; jet_i<(*jet_pt_groomed).size(); jet_i++)
		{
		
		  if (chosenLeadGroomedIndex<0 && DeltaR((*jet_eta_truth)[chosenLeadTruthJetIndex],(*jet_phi_truth)[chosenLeadTruthJetIndex],(*jet_eta_groomed)[jet_i],(*jet_phi_groomed)[jet_i])<0.9 && (*jet_emfrac_groomed)[jet_i]<0.99 && fabs((*jet_eta_groomed)[jet_i])<1.2)
		    {
		      chosenLeadGroomedIndex=jet_i;
		      continue;
		    }     
		} // end loop over jet_pt_groomed
	      
	      if (chosenLeadGroomedIndex == -99) // failed selection
		continue;	      
	      
	      leadGroomedIndex = chosenLeadGroomedIndex;
	      leadTruthIndex = chosenLeadTruthJetIndex;
	      leadTopoIndex = chosenLeadTopoJetIndex;
	      runNumberOut = runNumberIn;

	      mass = (*var_m_vec[2])[leadGroomedIndex]/1000.0 ;

	      if (applyMassWindow && (mass > mass_max && mass < mass_min))
		{
		  continue;
		}

	      int lead_subjet = 0;

	      if (subjetscalc)
		{

		  std::vector<int> subjet_idx = (*subjet_index).at(leadGroomedIndex); // only groomed ones.....
		  std::pair<int,int> subjet_leading = getTwoLeadingSubjets(subjet_idx,var_subjets_pt_vec);

		  // calculate the mass drop
		  if (subjet_idx.size() <= 1)
		    {
		      var_massdrop = -999;
		      var_yt = -999;
		      //break;
		    }
		  else
		    {
		      lead_subjet = subjet_leading.first;
		    
		      double pt1 =  (*var_subjets_pt_vec)[subjet_leading.first];
		      double pt2 =  (*var_subjets_pt_vec)[subjet_leading.second];
		      double eta_1 =  (*var_subjets_eta_vec)[subjet_leading.first];
		      double eta_2 =  (*var_subjets_eta_vec)[subjet_leading.second];
		      double phi_1 = (*var_subjets_phi_vec)[subjet_leading.first];
		      double phi_2 = (*var_subjets_phi_vec)[subjet_leading.second];
		      double e1 = (*var_subjets_E_vec)[subjet_leading.first];
		    
		      double subjet_mass =  (*var_subjets_m_vec)[subjet_leading.first];
		    
		      double mu12 = subjet_mass/(mass*1000);
		      var_massdrop = mu12;
		    
		      // momentum balance
		      double dRsub12 = DeltaR (eta_1, phi_1, eta_2, phi_2);
		      Float_t yt = (pt2*dRsub12)/(mass*1000);		    
		      yt*=yt;
		      var_yt = yt;
		    }
		} // end if(subjets)

	      else if (subjetspre)
		{
		  var_massdrop = var_massFraction_vec;//[leadGroomedIndex];
		  //std::cout << var_massFraction_vec << std::endl;
		  var_yt = var_ktycut2_vec;//[leadGroomedIndex];
		}

	      // tau21
	      var_Tau21[0]=(*var_Tau2_vec[0])[leadTruthIndex]/(*var_Tau1_vec[0])[leadTruthIndex];
	      var_Tau21[1]=(*var_Tau2_vec[1])[leadTopoIndex]/(*var_Tau1_vec[1])[leadTopoIndex];
	      var_Tau21[2]=(*var_Tau2_vec[2])[leadGroomedIndex]/(*var_Tau1_vec[2])[leadGroomedIndex];
	    
	      if (extendedVars)
		{
		  var_TauWTA2TauWTA1[0]=(*var_TauWTA2_vec[0])[leadTruthIndex]/(*var_TauWTA1_vec[0])[leadTruthIndex];
		  var_TauWTA2TauWTA1[2]=(*var_TauWTA2_vec[2])[leadGroomedIndex]/(*var_TauWTA1_vec[2])[leadGroomedIndex];
		  if (leadTopoIndex == -99)
		    var_TauWTA2TauWTA1[1]=-99;
		  else
		    var_TauWTA2TauWTA1[1]=(*var_TauWTA2_vec[1])[leadTopoIndex]/(*var_TauWTA1_vec[1])[leadTopoIndex];
		}
	    
	      setOutputVariables(extendedVars, leadTruthIndex, leadTopoIndex, leadGroomedIndex, lead_subjet);

	      if (k_factors.count(long(mc_channel_number)) > 0 )
		{
		  var_k_factor = k_factors[long(mc_channel_number)];
		}
	      else
		{
		  var_k_factor = 1.0;
		  //std::cout << "Warning: k-factor not known for " << mc_channel_number << endl;
		}
	      if (filt_eff.count(long(mc_channel_number)) > 0 )
		{
		  var_filter_eff = filt_eff[long(mc_channel_number)];
		}
	      else
		{
		  var_filter_eff = 1.0;
		  //std::cout << "Warning: filter eff not known for " << mc_channel_number << endl;
		}
	      if (xs.count(long(mc_channel_number)) > 0)
		{
		  var_xs = xs[long(mc_channel_number)];
		}
	      else
		{
		  var_xs = 1.0;
		  //std::cout << "Warning: xs not known for " << mc_channel_number << endl;
		}

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
	  for (std::map<int,float>::iterator it = NEvents_weighted.begin(); it!= NEvents_weighted.end(); it++)
	    ev_out << it->first << "," << NEvents_weighted[it->first] << std::endl;
	  ev_out.close();
	  delete outfile;
	  //delete intree;
	} // end loop of datatype
      // create the reweighting histogram
      std::string fname = algorithms.AlgoNames[i]+fileid_global+"/" + ss.str() + ".ptweights";
      createPtReweightFile(pt_reweight_arr[sampleType::BACKGROUND], pt_reweight_arr[sampleType::SIGNAL], fname);
    } // end loop over pt bins
    //} // end loop over algorithms


} // makeMassWindowFile()


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


vector<std::string> getListOfJetBranches(std::string &algorithm)
{
  // Ideally we want this stored in an XML file, but for now it'll have to be a standard text file because I'm short on time!
  vector<string> branches;
  ifstream in(algorithm+"_branches.txt");
  string line;
  // what OBto do about branches that have * in them?
  while (getline(in, line))
    {
      branches.push_back(rtrim(line));
    }
  // need to add some essential ones in case they get forgotten in that config file :)
  branches.push_back("RunNumber");
  branches.push_back("mc_channel_number");
  branches.push_back("mc_event_weight");
  in.close();
  return branches;
} // getListOfBranches()




std::string returnJetType(std::string & samplePrefix, std::string & groomalgo, bool addLC, int i)
{
  std::string jetType = "";
  switch (i)
    {
    case 0: // truth
      jetType="jet_CamKt12Truth_";
      break;
    case 1: // topo
      if (addLC)
	jetType = "jet_" + samplePrefix + "LCTopo_";
      else
	jetType = "jet_" + samplePrefix + "Topo_";
      break;
    default: // groomed
      if (addLC)
	jetType = "jet_" + samplePrefix +"LC" + groomalgo + "_";
      else
	jetType = "jet_" + samplePrefix + groomalgo + "_"; 	  
    }
  return jetType;
} // returnJetType


std::string returnSubJetType(std::string & samplePrefix, std::string & groomalgo, bool addLC)
{
  if (addLC)
    return std::string("jet_"+samplePrefix+"LC"+algorithms.subjetMap[groomalgo]+"_");
  return std::string("jet_"+samplePrefix+algorithms.subjetMap[groomalgo]+"_");

}// return subjettype



void addJetsBranches(TTree * tree, std::string & groomalgo, bool signal)
{
  tree->Branch("leadTruthIndex",&leadTruthIndex, "leadTruthIndex/I");
  tree->Branch("leadTopoIndex",&leadTopoIndex, "leadTopoIndex/I");
  tree->Branch("leadGroomedIndex",&leadGroomedIndex, "leadGroomedIndex/I");
  //tree->Branch("pt_reweight",&pt_reweight, "pt_reweight/F");
  tree->Branch("normalisation",&normalisation, "normalisation/F");
  tree->Branch("NEvents",&NEvents,"NEvents/I");
  tree->Branch("RunNumber",&runNumberOut, "RunNumber/I");
  tree->Branch("k_factor", &var_k_factor, "k_factor/F");
  tree->Branch("filter_eff", &var_filter_eff, "filter_eff/F");
  tree->Branch("xs", &var_xs, "xs/F");
  //tree->Branch("NEvents_weighted",&NEvents_weighted,"NEvents_weighted/F");
}// addJetsBranches

void setJetsBranches(TChain * tree, std::string &groomalgo, bool signal, std::string & groomIdx, bool extendedVars)
{
  std::string samplePrefix = "";
  bool addLC = false;
  samplePrefix = algorithms.AlgoPrefix[groomIdx];
  
  if (algorithms.AlgoType[groomIdx].find("recluster") == std::string::npos) // we're doing reclustering
    addLC = true; // just add the LC to the name

  tree->SetBranchAddress("mc_event_weight",&mc_event_weight);
  tree->SetBranchAddress("mc_channel_number", &mc_channel_number);
  tree->SetBranchAddress("RunNumber", &runNumberIn);

  TObjArray * brancharray = tree->GetListOfBranches();

  for (int i = 0; i < jetType::MAX; i++) // truth, topo, groomed
    {

      std::string jetType = returnJetType(samplePrefix, groomalgo, addLC,i); //set to truth/ topo/ groomed
      // I want to call GetListOfBranches once, but the ROOT site is down so I'm not sure what datatype it is....
      // whenever something fails we should set a bool so that we know not to add this to the output tree
      if(brancharray->FindObject(std::string(jetType+"E").c_str()))
	tree->SetBranchAddress(std::string(jetType+"E").c_str(),&var_E_vec.at(i));
      else
	std::cout << "branch not found: " << std::string(jetType+"E") << std::endl;
      if(brancharray->FindObject(std::string(jetType+"pt").c_str()))
	tree->SetBranchAddress(std::string(jetType+"pt").c_str(),&var_pt_vec.at(i));
      else
	std::cout << "branch not found: " << std::string(jetType+"pt") << std::endl;
      if(brancharray->FindObject(std::string(jetType+"m").c_str()))
	tree->SetBranchAddress(std::string(jetType+"m").c_str(),&var_m_vec.at(i));
      else
	std::cout << "branch not found: " << std::string(jetType+"m") << std::endl;
      if(brancharray->FindObject(std::string(jetType+"eta").c_str()))
	tree->SetBranchAddress(std::string(jetType+"eta").c_str(),&var_eta_vec.at(i));
      else
	std::cout << "branch not found: " << std::string(jetType+"eta") << std::endl;
      if(brancharray->FindObject(std::string(jetType+"phi").c_str()))
	tree->SetBranchAddress(std::string(jetType+"phi").c_str(),&var_phi_vec.at(i));
      else
	std::cout << "branch not found: " << std::string(jetType+"phi") << std::endl;
      //tree->SetBranchStatus(std::string(jetType+"constit_index").c_str(),1);
      //tree->SetBranchAddress(std::string(jetType+"constit_index").c_str(),&var_constit_index.at(i));
      if (i != 0 && brancharray->FindObject(std::string(jetType+"emfrac").c_str())) // no truth emfrac
	tree->SetBranchAddress(std::string(jetType+"emfrac").c_str(),&var_emfrac_vec.at(i));
      else
	std::cout << "branch not found: " << std::string(jetType+"emfrac") << std::endl;
      if(brancharray->FindObject(std::string(jetType+"Tau1").c_str()))
	tree->SetBranchAddress(std::string(jetType+"Tau1").c_str(),&var_Tau1_vec.at(i));
      else
	std::cout << "branch not found: " << std::string(jetType+"Tau1") << std::endl;
      if(brancharray->FindObject(std::string(jetType+"Tau2").c_str()))
	tree->SetBranchAddress(std::string(jetType+"Tau2").c_str(),&var_Tau2_vec.at(i));
      else
	std::cout << "branch not found: " << std::string(jetType+"Tau2") << std::endl;
      if(brancharray->FindObject(std::string(jetType+"Tau3").c_str()))
	tree->SetBranchAddress(std::string(jetType+"Tau3").c_str(),&var_Tau3_vec.at(i));
      else
	std::cout << "branch not found: " << std::string(jetType+"Tau3") << std::endl;
      if(brancharray->FindObject(std::string(jetType+"WIDTH").c_str()))
	tree->SetBranchAddress(std::string(jetType+"WIDTH").c_str(),&var_WIDTH_vec.at(i));
      else
	std::cout << "branch not found: " << std::string(jetType+"WIDTH") << std::endl;
      if(brancharray->FindObject(std::string(jetType+"SPLIT12").c_str()))
	tree->SetBranchAddress(std::string(jetType+"SPLIT12").c_str(),&var_SPLIT12_vec.at(i));
      else
	std::cout << "branch not found: " << std::string(jetType+"SPLIT12") << std::endl;
      if(brancharray->FindObject(std::string(jetType+"SPLIT23").c_str()))
	tree->SetBranchAddress(std::string(jetType+"SPLIT23").c_str(),&var_SPLIT23_vec.at(i));
      else
	std::cout << "branch not found: " << std::string(jetType+"SPLIT23") << std::endl;
      if(brancharray->FindObject(std::string(jetType+"SPLIT34").c_str()))
	tree->SetBranchAddress(std::string(jetType+"SPLIT34").c_str(),&var_SPLIT34_vec.at(i));
      else
	std::cout << "branch not found: " << std::string(jetType+"SPLIT34") << std::endl;
      if(brancharray->FindObject(std::string(jetType+"Dip12").c_str()))
	tree->SetBranchAddress(std::string(jetType+"Dip12").c_str(),&var_Dip12_vec.at(i));
      else
	std::cout << "branch not found: " << std::string(jetType+"Dip12") << std::endl;
      if(brancharray->FindObject(std::string(jetType+"Dip13").c_str()))
	tree->SetBranchAddress(std::string(jetType+"Dip13").c_str(),&var_Dip13_vec.at(i));
      else
	std::cout << "branch not found: " << std::string(jetType+"Dip13") << std::endl;
      if(brancharray->FindObject(std::string(jetType+"Dip23").c_str()))
	tree->SetBranchAddress(std::string(jetType+"Dip23").c_str(),&var_Dip23_vec.at(i));
      else
	std::cout << "branch not found: " << std::string(jetType+"Dip23") << std::endl;
      if(brancharray->FindObject(std::string(jetType+"DipExcl12").c_str()))
	tree->SetBranchAddress(std::string(jetType+"DipExcl12").c_str(),&var_DipExcl12_vec.at(i));
      else
	std::cout << "branch not found: " << std::string(jetType+"DipExcl12") << std::endl;
      if(brancharray->FindObject(std::string(jetType+"PlanarFlow").c_str()))
	tree->SetBranchAddress(std::string(jetType+"PlanarFlow").c_str(),&var_PlanarFlow_vec.at(i));
      else
	std::cout << "branch not found: " << std::string(jetType+"PlanarFlow") << std::endl;
      if(brancharray->FindObject(std::string(jetType+"Angularity").c_str()))
	tree->SetBranchAddress(std::string(jetType+"Angularity").c_str(),&var_Angularity_vec.at(i));
      else
	std::cout << "branch not found: " << std::string(jetType+"Angularity") << std::endl;
      if(brancharray->FindObject(std::string(jetType+"QW").c_str()))
	tree->SetBranchAddress(std::string(jetType+"QW").c_str(),&var_QW_vec.at(i));
      else
	std::cout << "branch not found: " << std::string(jetType+"QW") << std::endl;
      if(brancharray->FindObject(std::string(jetType+"PullMag").c_str()))
	tree->SetBranchAddress(std::string(jetType+"PullMag").c_str(),&var_PullMag_vec.at(i));
      else
	std::cout << "branch not found: " << std::string(jetType+"PullMag") << std::endl;
      if(brancharray->FindObject(std::string(jetType+"PullPhi").c_str()))
	tree->SetBranchAddress(std::string(jetType+"PullPhi").c_str(),&var_PullPhi_vec.at(i));
      else
	std::cout << "branch not found: " << std::string(jetType+"PullPhi") << std::endl;
      if(brancharray->FindObject(std::string(jetType+"Pull_C00").c_str()))
	tree->SetBranchAddress(std::string(jetType+"Pull_C00").c_str(),&var_Pull_C00_vec.at(i));
      else
	std::cout << "branch not found: " << std::string(jetType+"Pull_C00") << std::endl;
      if(brancharray->FindObject(std::string(jetType+"Pull_C01").c_str()))
	tree->SetBranchAddress(std::string(jetType+"Pull_C01").c_str(),&var_Pull_C01_vec.at(i));
      else
	std::cout << "branch not found: " << std::string(jetType+"Pull_C01") << std::endl;
      if(brancharray->FindObject(std::string(jetType+"Pull_C10").c_str()))
	tree->SetBranchAddress(std::string(jetType+"Pull_C10").c_str(),&var_Pull_C10_vec.at(i));
      else
	std::cout << "branch not found: " << std::string(jetType+"Pull_C10") << std::endl;
      if(brancharray->FindObject(std::string(jetType+"Pull_C11").c_str()))
	tree->SetBranchAddress(std::string(jetType+"Pull_C11").c_str(),&var_Pull_C11_vec.at(i));
      else
	std::cout << "branch not found: " << std::string(jetType+"Pull_C11") << std::endl;


      if (extendedVars)
	{
	  if(brancharray->FindObject(std::string(jetType+"TauWTA1").c_str()))
	    tree->SetBranchAddress(std::string(jetType+"TauWTA1").c_str(),&var_TauWTA1_vec.at(i));
	  else
	    std::cout << "branch not found: " << std::string(jetType+"TauWTA1") << std::endl;
	  if(brancharray->FindObject(std::string(jetType+"TauWTA2").c_str()))
	    tree->SetBranchAddress(std::string(jetType+"TauWTA2").c_str(),&var_TauWTA2_vec.at(i));
	  else
	    std::cout << "branch not found: " << std::string(jetType+"TauWTA2") << std::endl;
	  if(brancharray->FindObject(std::string(jetType+"ZCUT12").c_str()))
	    tree->SetBranchAddress(std::string(jetType+"ZCUT12").c_str(),&var_ZCUT12_vec.at(i));
	  else
	    std::cout << "branch not found: " << std::string(jetType+"ZCUT12") << std::endl;
	}


      /*setAddress(tree,std::string(jetType+"E"),var_E_vec.at(i));
	setAddress(tree,std::string(jetType+"pt"),var_pt_vec.at(i));
	setAddress(tree,std::string(jetType+"m"),var_m_vec.at(i));
	setAddress(tree,std::string(jetType+"eta"),var_eta_vec.at(i));
	setAddress(tree,std::string(jetType+"phi"),var_phi_vec.at(i));
	if (i != 0)
	setAddress(tree,std::string(jetType+"emfrac"),var_emfrac_vec.at(i));
	setAddress(tree,std::string(jetType+"Tau1"),var_Tau1_vec.at(i));
	setAddress(tree,std::string(jetType+"Tau2"),var_Tau2_vec.at(i));
	setAddress(tree,std::string(jetType+"Tau3"),var_Tau3_vec.at(i));
	setAddress(tree,std::string(jetType+"WIDTH"),var_WIDTH_vec.at(i));
	setAddress(tree,std::string(jetType+"SPLIT12"),var_SPLIT12_vec.at(i));
	setAddress(tree,std::string(jetType+"SPLIT23"),var_SPLIT23_vec.at(i));
	setAddress(tree,std::string(jetType+"SPLIT34"),var_SPLIT34_vec.at(i));
	setAddress(tree,std::string(jetType+"Dip12"),var_Dip12_vec.at(i));
	setAddress(tree,std::string(jetType+"Dip13"),var_Dip13_vec.at(i));
	setAddress(tree,std::string(jetType+"Dip23"),var_Dip23_vec.at(i));
	setAddress(tree,std::string(jetType+"DipExcl12"),var_DipExcl12_vec.at(i));
	setAddress(tree,std::string(jetType+"PlanarFlow"),var_PlanarFlow_vec.at(i));
	setAddress(tree,std::string(jetType+"Angularity"),var_Angularity_vec.at(i));
	setAddress(tree,std::string(jetType+"QW"),var_QW_vec.at(i));
	setAddress(tree,std::string(jetType+"PullMag"),var_PullMag_vec.at(i));
	setAddress(tree,std::string(jetType+"PullPhi"),var_PullPhi_vec.at(i));
	setAddress(tree,std::string(jetType+"Pull_C00"),var_Pull_C00_vec.at(i));
	setAddress(tree,std::string(jetType+"Pull_C01"),var_Pull_C01_vec.at(i));
	setAddress(tree,std::string(jetType+"Pull_C10"),var_Pull_C10_vec.at(i));
	setAddress(tree,std::string(jetType+"Pull_C11"),var_Pull_C11_vec.at(i));
	if (extendedVars)
	{
	setAddress(tree,std::string(jetType+"TauWTA1"),var_TauWTA1_vec.at(i));
	setAddress(tree,std::string(jetType+"TauWTA2"),var_TauWTA2_vec.at(i));
	setAddress(tree,jetType+"ZCUT12",var_ZCUT12_vec.at(i));
	}*/


    } // end for loop over topo/truth/groom

  std::string jetType = returnJetType( samplePrefix, groomalgo, addLC, 2); //set to truth/ topo/ groomed
  if (subjetspre)
    {

      if(brancharray->FindObject(std::string(jetType+"config_massFraction").c_str() ) )
	tree->SetBranchAddress(std::string(jetType+"config_massFraction").c_str(),&var_massFraction_vec);
      else
	std::cout << "branch not found: " << std::string(jetType+"config_massFraction") << std::endl;
      if(brancharray->FindObject(std::string(jetType+"config_ktycut2").c_str()))
	tree->SetBranchAddress(std::string(jetType+"config_ktycut2").c_str(),&var_ktycut2_vec);
      else
	std::cout << "branch not found: " << std::string(jetType+"config_ktycut2") << std::endl;
    }




  if (!subjetscalc)
    return;


  std::string subjetType = returnSubJetType(samplePrefix, groomalgo, addLC);
  
  tree->SetBranchStatus(std::string(algorithms.subjetIndex[groomalgo]).c_str(),1);
  tree->SetBranchAddress(std::string(algorithms.subjetIndex[groomalgo]).c_str(),&subjet_index);

  if(!tree->GetListOfBranches()->FindObject(std::string(algorithms.subjetIndex[groomalgo]).c_str())) {
    std::cout << "subjet branch is not here" << std::endl;
    return;
  }
  
  tree->SetBranchAddress(std::string(subjetType+"E").c_str(),&var_subjets_E_vec);
  tree->SetBranchAddress(std::string(subjetType+"pt").c_str(),&var_subjets_pt_vec);
  tree->SetBranchAddress(std::string(subjetType+"m").c_str(),&var_subjets_m_vec);
  tree->SetBranchAddress(std::string(subjetType+"eta").c_str(),&var_subjets_eta_vec);
  tree->SetBranchAddress(std::string(subjetType+"phi").c_str(),&var_subjets_phi_vec);
  
}//setJetsBranches

void addSubJets(TTree * tree, std::string & groomalgo, bool signal, std::string &  groomIdx)
{

  std::string samplePrefix = "";
  bool addLC = false;
  samplePrefix = algorithms.AlgoPrefix[groomIdx];

  if (algorithms.AlgoType[groomIdx].find("recluster") == std::string::npos) // we're doing reclustering
    addLC = true; // just add the LC to the name

  std::string jetType = returnJetType( samplePrefix, groomalgo, addLC, 2); //set to truth/ topo/ groomed
  std::string subjetType = returnSubJetType(samplePrefix, groomalgo, addLC);
  if (subjetscalc)
    {
      tree->Branch(std::string(subjetType+"E").c_str(),&var_subjets_E,std::string(subjetType+"E/F").c_str());
      tree->Branch(std::string(subjetType+"pt").c_str(),&var_subjets_pt,std::string(subjetType+"pt/F").c_str());
      tree->Branch(std::string(subjetType+"m").c_str(),&var_subjets_m,std::string(subjetType+"m/F").c_str());
      tree->Branch(std::string(subjetType+"eta").c_str(),&var_subjets_eta,std::string(subjetType+"eta/F").c_str());
      tree->Branch(std::string(subjetType+"phi").c_str(),&var_subjets_phi,std::string(subjetType+"phi/F").c_str());
    }

  tree->Branch(std::string(jetType+"massdrop").c_str(),&var_massdrop, std::string(jetType+"massdrop/F").c_str());
  tree->Branch(std::string(jetType+"yt").c_str(),&var_yt,std::string(jetType+"yt/F").c_str());

} //addSubJets()


//void getBranchesSelection(TTree * tree, std::string & algorithm)
void setSelectionVectors(bool signal, std::string & algorithm)
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

void initVectors(bool extendedVars)
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
  
      if (extendedVars)
	{
	  var_TauWTA1_vec[i] = 0;
	  var_TauWTA2_vec[i] = 0;
	  var_ZCUT12_vec[i] = 0;
	}

    }  
  var_massFraction_vec = 0;
  var_ktycut2_vec = 0;

  subjet_index = 0;    
  var_subjets_E_vec = 0;
  var_subjets_pt_vec = 0;
  var_subjets_m_vec = 0;
  var_subjets_eta_vec = 0;
  var_subjets_phi_vec = 0;

} // initVEctors

void setOutputVariables(bool extendedVars, int jet_idx_truth, int jet_idx_topo, int jet_idx_groomed, int subjet_idx)
{
  //clearOutputVariables();
  int jet_idx = 0;
  mc_event_weight_out = mc_event_weight;
  mc_channel_number_out = mc_channel_number;
  
  for (int x = 0; x < jetType::MAX ; x++)
    {
      switch (x)
	{
	case 0:
	  jet_idx = jet_idx_truth;
	  break;
	case 1:
	  jet_idx = jet_idx_topo;
	  break;
	case 2:
	  jet_idx = jet_idx_groomed;
	  break;
	}
      
      if (jet_idx == -99)
	continue;

      var_E[x]=(*var_E_vec[x])[jet_idx];
      var_pt[x]=(*var_pt_vec[x])[jet_idx];
      var_m[x]=(*var_m_vec[x])[jet_idx];
      var_eta[x]=(*var_eta_vec[x])[jet_idx];
      var_phi[x]=(*var_phi_vec[x])[jet_idx];
      if (x!=0)
	var_emfrac[x]=(*var_emfrac_vec[x])[jet_idx];
      var_Tau1[x]=(*var_Tau1_vec[x])[jet_idx];
      var_Tau2[x]=(*var_Tau2_vec[x])[jet_idx];
      var_Tau3[x]=(*var_Tau3_vec[x])[jet_idx];
      var_WIDTH[x]=(*var_WIDTH_vec[x])[jet_idx];
      var_SPLIT12[x]=(*var_SPLIT12_vec[x])[jet_idx];
      var_SPLIT23[x]=(*var_SPLIT23_vec[x])[jet_idx];
      var_SPLIT34[x]=(*var_SPLIT34_vec[x])[jet_idx];
      var_Dip12[x]=(*var_Dip12_vec[x])[jet_idx];
      var_Dip13[x]=(*var_Dip13_vec[x])[jet_idx];
      var_Dip23[x]=(*var_Dip23_vec[x])[jet_idx];
      var_DipExcl12[x]=(*var_DipExcl12_vec[x])[jet_idx];
      var_PlanarFlow[x]=(*var_PlanarFlow_vec[x])[jet_idx];
      var_Angularity[x]=(*var_Angularity_vec[x])[jet_idx];
      var_QW[x]=(*var_QW_vec[x])[jet_idx];
      var_PullMag[x]=(*var_PullMag_vec[x])[jet_idx];
      var_PullPhi[x]=(*var_PullPhi_vec[x])[jet_idx];
      var_Pull_C00[x]=(*var_Pull_C00_vec[x])[jet_idx];
      var_Pull_C01[x]=(*var_Pull_C01_vec[x])[jet_idx];
      var_Pull_C10[x]=(*var_Pull_C10_vec[x])[jet_idx];
      var_Pull_C11[x]=(*var_Pull_C11_vec[x])[jet_idx];
      // tau21 is set in the main loop, not here, because we have to calculate it

      if (extendedVars)
	{
	  var_TauWTA1[x]=(*var_TauWTA1_vec[x])[jet_idx];
	  var_TauWTA2[x]=(*var_TauWTA2_vec[x])[jet_idx];
	  var_ZCUT12[x]=(*var_ZCUT12_vec[x])[jet_idx];
	  //var_massFraction[x]=(*var_massFraction_vec[x])[jet_idx];
	  //var_ytcut2[x]=(*var_ZCUT12_vec[x])[jet_idx];
	}
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

void clearOutputVariables()
{

  var_massdrop = 0;
  var_yt = 0;

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
  var_Tau21.clear();
  var_TauWTA2TauWTA1.clear();
  var_TauWTA1.clear();
  var_TauWTA2.clear();
  var_ZCUT12.clear();


} // clearOutputVariables

void resetOutputVariables()
{
  clearOutputVariables();
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
      var_Tau21.push_back(-999);
      var_TauWTA2TauWTA1.push_back(-999);
      var_TauWTA1.push_back(-999);
      var_TauWTA2.push_back(-999);
      var_ZCUT12.push_back(-999);
    }
} //resetOutputVariables


void setOutputBranches(TTree * tree, std::string & groomalgo, std::string & groomIdx, bool extendedVars)
{

  std::string samplePrefix = "";
  bool addLC = false;
  samplePrefix = algorithms.AlgoPrefix[groomIdx];

  if (algorithms.AlgoType[groomIdx].find("recluster") == std::string::npos) // we're doing reclustering
    addLC = true; // just add the LC to the name

  tree->Branch("mc_event_weight",&mc_event_weight_out,"mc_event_weight/F");
  tree->Branch("mc_channel_number", &mc_channel_number_out,"mc_channel_number/I");

  for (int i = 0; i < jetType::MAX; i++) // truth, topo, groomed
    {
      
      std::string jetType = returnJetType(samplePrefix, groomalgo, addLC,i); //set to truth/ topo/ groomed

      tree->Branch(std::string(jetType+"E").c_str(),&var_E.at(i),std::string(jetType+"E/F").c_str());
      tree->Branch(std::string(jetType+"pt").c_str(),&var_pt.at(i),std::string(jetType+"pt/F").c_str());
      tree->Branch(std::string(jetType+"m").c_str(),&var_m.at(i),std::string(jetType+"m/F").c_str());
      tree->Branch(std::string(jetType+"eta").c_str(),&var_eta.at(i),std::string(jetType+"eta/F").c_str());
      tree->Branch(std::string(jetType+"phi").c_str(),&var_phi.at(i),std::string(jetType+"phi/F").c_str());
      if (i != 0)
	tree->Branch(std::string(jetType+"emfrac").c_str(),&var_emfrac.at(i),std::string(jetType+"emfrac"+"/F").c_str());
      tree->Branch(std::string(jetType+"Tau1").c_str(),&var_Tau1.at(i),std::string(jetType+"Tau1/F").c_str());
      tree->Branch(std::string(jetType+"Tau2").c_str(),&var_Tau2.at(i),std::string(jetType+"Tau2/F").c_str());
      tree->Branch(std::string(jetType+"Tau3").c_str(),&var_Tau3.at(i),std::string(jetType+"Tau3/F").c_str());
      tree->Branch(std::string(jetType+"WIDTH").c_str(),&var_WIDTH.at(i),std::string(jetType+"WIDTH/F").c_str());
      tree->Branch(std::string(jetType+"SPLIT12").c_str(),&var_SPLIT12.at(i),std::string(jetType+"SPLIT12/F").c_str());
      tree->Branch(std::string(jetType+"SPLIT23").c_str(),&var_SPLIT23.at(i),std::string(jetType+"SPLIT23/F").c_str());
      tree->Branch(std::string(jetType+"SPLIT34").c_str(),&var_SPLIT34.at(i),std::string(jetType+"SPLIT34/F").c_str());
      tree->Branch(std::string(jetType+"Dip12").c_str(),&var_Dip12.at(i),std::string(jetType+"Dip12/F").c_str());
      tree->Branch(std::string(jetType+"Dip13").c_str(),&var_Dip13.at(i),std::string(jetType+"Dip13/F").c_str());
      tree->Branch(std::string(jetType+"Dip23").c_str(),&var_Dip23.at(i),std::string(jetType+"Dip23/F").c_str());
      tree->Branch(std::string(jetType+"DipExcl12").c_str(),&var_DipExcl12.at(i),std::string(jetType+"DipExcl12/F").c_str());
      tree->Branch(std::string(jetType+"PlanarFlow").c_str(),&var_PlanarFlow.at(i),std::string(jetType+"PlanarFlow/F").c_str());
      tree->Branch(std::string(jetType+"Angularity").c_str(),&var_Angularity.at(i),std::string(jetType+"Angularity/F").c_str());
      tree->Branch(std::string(jetType+"QW").c_str(),&var_QW.at(i),std::string(jetType+"QW/F").c_str());
      tree->Branch(std::string(jetType+"PullMag").c_str(),&var_PullMag.at(i),std::string(jetType+"PullMag/F").c_str());
      tree->Branch(std::string(jetType+"PullPhi").c_str(),&var_PullPhi.at(i),std::string(jetType+"PullPhi/F").c_str());
      tree->Branch(std::string(jetType+"Pull_C00").c_str(),&var_Pull_C00.at(i),std::string(jetType+"Pull_C00/F").c_str());
      tree->Branch(std::string(jetType+"Pull_C01").c_str(),&var_Pull_C01.at(i),std::string(jetType+"Pull_C01/F").c_str());
      tree->Branch(std::string(jetType+"Pull_C10").c_str(),&var_Pull_C10.at(i),std::string(jetType+"Pull_C10/F").c_str());
      tree->Branch(std::string(jetType+"Pull_C11").c_str(),&var_Pull_C11.at(i),std::string(jetType+"Pull_C11/F").c_str());


      if (extendedVars)
	{
	  tree->Branch(std::string(jetType+"TauWTA1").c_str(),&var_TauWTA1.at(i),std::string(jetType+"TauWTA1/F").c_str());
	  tree->Branch(std::string(jetType+"TauWTA2").c_str(),&var_TauWTA2.at(i),std::string(jetType+"TauWTA2/F").c_str());
	  tree->Branch(std::string(jetType+"ZCUT12").c_str(),&var_ZCUT12.at(i),std::string(jetType+"ZCUT12/F").c_str());
	  
	  tree->Branch(std::string(returnJetType(samplePrefix, groomalgo, addLC,0)+"TauWTA2/TauWTA1").c_str(),&var_TauWTA2TauWTA1.at(0),std::string(jetType+"TauWTA2TauWTA1/F").c_str());
	  tree->Branch(std::string(returnJetType(samplePrefix, groomalgo, addLC,1)+"TauWTA2/TauWTA1").c_str(),&var_TauWTA2TauWTA1.at(1),std::string(jetType+"TauWTA2TauWTA1/F").c_str());
	  tree->Branch(std::string(jetType+"TauWTA2/TauWTA1").c_str(),&var_TauWTA2TauWTA1.at(2),std::string(jetType+"TauWTA2TauWTA1/F").c_str());  
	  
	}

      tree->Branch(std::string(returnJetType(samplePrefix, groomalgo, addLC,0)+"Tau21").c_str(),&var_Tau21.at(0),std::string(jetType+"Tau21/F").c_str());
      tree->Branch(std::string(returnJetType(samplePrefix, groomalgo, addLC,1)+"Tau21").c_str(),&var_Tau21.at(1),std::string(jetType+"Tau21/F").c_str());
      tree->Branch(std::string(jetType+"Tau21").c_str(),&var_Tau21.at(2),std::string(jetType+"Tau21/F").c_str());  

    }
  
} // setOutputBranches



void scaleHists()
{
  //for (int ii=0; ii<nAlgos; ii++){//-1; ii++){
  //int i = algoMap[ii]; // this is here because we don't always run all algorithms, just a subsample...
	

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
    //qcd_PtReweight->Scale(1.0/qcd_PtReweight->Integral());
    //Wp_PtReweight->Scale(1.0/Wp_PtReweight->Integral());

		
    //}

} // scaleHists


void getMPV()
{
  //for (int ii=0; ii<nAlgos; ii++){//-1; ii++){
  //int i = algoMap[ii];
    for (int j=0; j<nPtBins; j++){
      myMPV[j]=mpv(Wprime_Lead_CA12_mass[histType::GROOMEDJET][j]);
    }
    
    for (int j=0; j<nFineBins; j++){
      myMPV_finePt[j]=mpv(Wprime_finePtBin_mass[histType::GROOMEDJET][j]);
    }
    //}
} //getMPV

void setAddress(TChain * tree, std::string name, std::vector<Float_t> * var_vec)
{
  if(tree->GetListOfBranches()->FindObject(name.c_str()) )
    tree->SetBranchAddress(name.c_str(),&var_vec);
  else
    std::cout << "branch not found: " << name << std::endl;
  
} //setAddress


void readWeights()
{
  // read in the weights for the different samples.  Weights are based on RunNumbers, the file should have
  // RunNumber,SampleName,xs,k-factor,filter-efficiency
  // weights come from these two twiki pages: https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/BosonTaggingPaper2012#7th_July_2014_Mixture_with_small and https://twiki.cern.ch/twiki/bin/view/AtlasProtected/BoostedBosonTaggingD3PD#Properly_normalizing_MC_samples
  ifstream wf ("weightings_input.csv");
  string line;
  std::cout << "reading weights" << std::endl;
  while (getline( wf, line))
    {
      trim(line);
      std::vector<std::string> strs;
      boost::split(strs, line, boost::is_any_of(","));
      //std::cout << strs.size() << std::endl;
      int cntr = 0;
      long runNumber = 0;
      for (std::vector<std::string>::iterator it = strs.begin(); it != strs.end(); it ++)
	{
	  switch (cntr)
	    {
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
	    case 1:
	      // This will be the sample name, but not using it right now.
	      break;
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
	  //std::cout << (*it) << std::endl;
	}
    }

  wf.close();
  
} // readWeights

void createPtReweightFile(TH1F * bkg, TH1F * sig, std::string & fname)
{
  // create an output file containing a number of datapoints for each algorithm + pt range.
  // each datapoint contains bkg/signal for each bin in bkg and sig arguments
  // this output file will be read in by the plotting code later and used for reweighting the signal sample
  int bins = bkg->GetNbinsX();
  std::cout << "number of bins: " << bins << std::endl;
  ofstream out(fname);
  double weight;
  for (int b = 1; b <= bins; b++)
    {
      weight = bkg->GetBinContent(b);
      if (sig->GetBinContent(b) != 0)
	weight/=sig->GetBinContent(b);
      else
	{
	  std::cout << "oh no, we've been rumbled!  We have no signal in this bin and now we have to reweight by over 9000! But actually we're just going to make it 0." << std::endl;
	  weight = 0;
	}
      out << bkg->GetXaxis()->GetBinLowEdge(b) << "," << weight << endl;     
    }
  out.close();
  
} // end createPtReweightHistogram


// loads algorithms from config file
void Algorithms::load(const std::string & filename)
{
  // Create an empty property tree object
  using boost::property_tree::ptree;
  ptree pt;

  read_xml(filename,pt);
  BOOST_FOREACH( ptree::value_type & v, pt.get_child("config.algorithms"))
    {
      std::string name =  v.second.get<std::string>("Algorithm","");
      if (name != "")
	{
	  std::string groomalgo = v.second.get<std::string>("GroomingAlgorithm","");
	  std::cout << groomalgo << std::endl;
	  if (groomalgo == "")
	    {
	      std::cout << name.find("Topo") << std::endl;
	      groomalgo = name.substr(name.find("Topo"),name.length()-1);
	    }
	  AlgoNames[name] = groomalgo;
	  std::string jetalgo = v.second.get<std::string>("JetAlgorithm","");
	  if (jetalgo == "")
	    {
	      jetalgo = name.substr(0,name.find("Topo"));
	      if (jetalgo.find("LC") != std::string::npos)
		jetalgo.erase(jetalgo.end()-2,jetalgo.end());
	    }
	  AlgoPrefix[name] = jetalgo;
	  AlgoList[name] = v.second.get<std::string>("Abbreviation","");
	  AlgoListN[name] = v.second.get<std::string>("PlotLabel","");
	  AlgoType[name] = v.second.get<std::string>("Type","NONE");
	  
	  subjetMap[groomalgo] = v.second.get<std::string>("SubjetGroomingAlgorithm");
	  subjetIndex[groomalgo] = v.second.get<std::string>("SubjetIndexBranch");
	  binLabel[groomalgo] = v.second.get<std::string>("BinLabel","");
	}
    };

} // algorithms::load()

