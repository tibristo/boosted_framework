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


using namespace std;

bool ComparePt(TLorentzVector a, TLorentzVector b) { return a.Pt() > b.Pt(); }


int main( int argc, char * argv[] ) {
  
  gROOT->ProcessLine("#include <vector>");
  gInterpreter->GenerateDictionary("vector<vector<int> >", "vector");
  gInterpreter->GenerateDictionary("vector<vector<double> >", "vector");
  
  AtlasStyle();
  bool filtering = false, trimmed = false, pruned = false, recluster = false, allsamples = false; // we may not want to run this for all at once, so I added these variables to indicate what we're running.  Currently this just looks at the file name, but it might be better to add an arg to the command line indicating which ones we're running, something like -f (split/filter), -t (trim), -ft (both). TODO
  int filteridx = 0, trimidx = 0, pruneidx = 0, reclusteridx = 0; // continuing from the idea above, we set the InputFile[idx] based on the input files that we give the program
  
  //Make an array of TFiles with all the relevant inputs and add them to my array of trees
  for (int nArg=1; nArg < argc; nArg++) {
    //cout << nArg << " " << argv[nArg] << endl;
    inputFile[nArg-1] = new TFile(argv[nArg], "READ");
    inputTree[nArg-1] = ( TTree* ) inputFile[nArg-1]->Get( "physics" );    
    std::string arg  = std::string(argv[nArg]);
    std::transform(arg.begin(), arg.end(), arg.begin(), ::tolower);
    vector<int> subset;
    if (arg.find("all")!=string::npos)
      {
	allsamples = true;
	filteridx = trimidx = pruneidx = reclusteridx = 0;
	cout << "all samples is true" << endl;
	subset = {0,1,2,3,4,5,6,7,8,9};
	algoMap.insert(algoMap.end(),subset.begin(),subset.end());
      }
    else if (false)//arg.find("filter")!=string::npos || arg.find("split")!=string::npos)
      {
	if (!filtering)
	  {
	    filteridx = nArg-1;
	    subset = {1,2};//{0,1,2};
	    algoMap.insert(algoMap.end(),subset.begin(),subset.end());
	  }
	filtering = true;
      }
    else if (true)//arg.find("trim")!=string::npos)
      {
	if(!trimmed)
	  {
	    trimidx = nArg-1;
	    subset = {3,4};
	    algoMap.insert(algoMap.end(),subset.begin(),subset.end());
	  }
	trimmed = true;
      }
    else if (arg.find("prune")!=string::npos)
      {
	if (!pruned)
	  {
	    pruneidx = nArg-1;
	    subset = {0,5,6};
	    algoMap.insert(algoMap.end(),subset.begin(),subset.end());
	  }
	pruned = true;
      }
    else if (arg.find("reclust")!=string::npos)
      {
	if (!recluster)
	  {
	    reclusteridx = nArg -1;	
	    subset = {7,8,9}; // removed 0 for now
	    algoMap.insert(algoMap.end(),subset.begin(),subset.end());
	  }
	recluster = true;
      }
    if (!subset.empty())
      {
	for (std::vector<int>::iterator it = subset.begin(); it != subset.end(); it++)
	  fileMap[(*it)] = nArg-1; // what about index 0?  That might cause problems because we set it for ALL files!
      }
    
  }

  defineStrings(AlgoList, binLabel, pTbins, finePtBins);
  createHistos();
  

  if (filtering || allsamples)
    {
      //runAlgorithm(inputTree[filteridx], inputTree[filteridx+1], "TopoSplitFilteredMu67SmallR0YCut9", 0);
      runAlgorithm(inputTree[filteridx], inputTree[filteridx+1], "TopoSplitFilteredMu67SmallR0YCut9", 1);
      runAlgorithm(inputTree[filteridx], inputTree[filteridx+1], "TopoSplitFilteredMu100SmallR30YCut4", 2);
      //runAlgorithm(inputTree[filteridx], inputTree[filteridx+1], "TrimmedPtFrac5SmallR30", 2);
    }
  
  if (trimmed || allsamples)
    {
      runAlgorithm(inputTree[trimidx], inputTree[trimidx+1], "TopoTrimmedPtFrac5SmallR30", 3);
      runAlgorithm(inputTree[trimidx], inputTree[trimidx+1], "TopoTrimmedPtFrac5SmallR20", 4);
    }
  
  if (pruned || allsamples)
    {
      runAlgorithm(inputTree[pruneidx], inputTree[pruneidx+1], "TopoPrunedCaRcutFactor50Zcut10", 5);
      runAlgorithm(inputTree[pruneidx], inputTree[pruneidx+1], "TopoPrunedCaRcutFactor50Zcut20", 6);
    }
  
  //reclustering

  if (false)//recluster || allsamples)
    {
      runAlgorithm(inputTree[reclusteridx], inputTree[reclusteridx+1], "AntiKt2LCTopo", 7);
      runAlgorithm(inputTree[reclusteridx], inputTree[reclusteridx+1], "AntiKt3LCTopo", 8);      
      runAlgorithm(inputTree[reclusteridx], inputTree[reclusteridx+1], "AntiKt4LCTopo", 9);
    }
  

  //Plot things
  //Normalise all histograms of interest
  // this is not the right way to do the normalisation according to this - https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/BoostedBosonTaggingD3PD#MC_only_studies

  // NEvents = Sum[mcevt_weight[0][0]*PileupWeight
  // Note, for NEvents for SherpaW + jets we need to call getNormSherpaW() because the weighting is slightly different
  // Normalization = CrossSection(at NLO) * Luminosity / NEvents


  for (int ii=0; ii<nAlgos; ii++){//-1; ii++){
    int i = algoMap[ii]; // this is here because we don't always run all algorithms, just a subsample...
    
    qcd_Lead_CA12_pt[i]->Scale(1.0/qcd_Lead_CA12_pt[i]->Integral());  
    Wprime_Lead_CA12_pt[i]->Scale(1.0/Wprime_Lead_CA12_pt[i]->Integral());
    Wprime_Lead_CA12_scaled_pt[i]->Scale(1.0/Wprime_Lead_CA12_scaled_pt[i]->Integral());
    
    for (int j=0; j<nPtBins; j++){
      
      qcd_Lead_CA12_mass[i][j]->Scale(1.0/qcd_Lead_CA12_mass[i][j]->Integral());
      Wprime_Lead_CA12_mass[i][j]->Scale(1.0/Wprime_Lead_CA12_mass[i][j]->Integral());
    }

    for (int j=0; j<nFineBins; j++){

      qcd_finePtBin_mass[i][j]->Scale(1.0/qcd_finePtBin_mass[i][j]->Integral());
      Wprime_finePtBin_mass[i][j]->Scale(1.0/Wprime_finePtBin_mass[i][j]->Integral());
    }

    //qcd_PtReweight[i]->Scale(1.0/qcd_PtReweight[i]->Integral());
    //Wp_PtReweight[i]->Scale(1.0/Wp_PtReweight[i]->Integral());
    Wprime_Lead_CA12_scaled_pt[i]->Scale(1.0/Wprime_Lead_CA12_scaled_pt[i]->Integral());


  }
  
  makePlots();

  //NOW THAT I HAVE ALL MY MASS PLOTS NORMALISED, I WANT TO GET THE MPV AND MASS WINDOW ANALYSIS RIGHT

  // QCD VECTOR: qcd_Lead_CA12_mass[i][j]
  // WPRIME VECTOR: Wprime_Lead_CA12_mass[i][j]

  //1. Get the MPV

  for (int ii=0; ii<nAlgos; ii++){//-1; ii++){
    int i = algoMap[ii];
    for (int j=0; j<nPtBins; j++){
      myMPV[i][j]=mpv(Wprime_Lead_CA12_mass[i][j]);
    }

    for (int j=0; j<nFineBins; j++){
      myMPV_finePt[i][j]=mpv(Wprime_finePtBin_mass[i][j]);
    }
  }
  
  //2. Get the mass window which gives 68% W mass efficiency 

  for (int ii=0; ii<nAlgos; ii++){//-1; ii++){
    int i = algoMap[ii];
    for (int j=0; j<nPtBins; j++){
      Qw(WidthMassWindow[i][j],TopEdgeMassWindow[i][j],Wprime_Lead_CA12_mass[i][j], 0.68);
      //std::cout << "TopEdge: " << TopEdgeMassWindow[i][j] << std::endl;
      //one and two are outputs, three and four are inputs
      BottomEdgeMassWindow[i][j]=TopEdgeMassWindow[i][j]-WidthMassWindow[i][j];    
      //cout << AlgoList[i] << ": top edge " << TopEdgeMassWindow[i][j] << " bottom edge " << BottomEdgeMassWindow[i][j] << " mass window " << WidthMassWindow[i][j] << endl;
    }

    for (int j=0; j<nFineBins; j++){
      Qw(WidthMassWindow_finePt[i][j],TopEdgeMassWindow_finePt[i][j],Wprime_finePtBin_mass[i][j], 0.68);
      //one and two are outputs, three and four are inputs
      BottomEdgeMassWindow_finePt[i][j]=TopEdgeMassWindow_finePt[i][j]-WidthMassWindow_finePt[i][j];    
      //cout << AlgoList[i] << ": top edge " << TopEdgeMassWindow[i][j] << " bottom edge " << BottomEdgeMassWindow[i][j] << " mass window " << WidthMassWindow[i][j] << endl;
    }
   
  }

  //3. Check background fraction in this window

  for (int ii=0; ii<nAlgos; ii++){//-1; ii++){
    int i = algoMap[ii];
    for (int j=0; j<nPtBins; j++){
      QCDfrac[i][j]=qcd_Lead_CA12_mass[i][j]->Integral(qcd_Lead_CA12_mass[i][j]->FindBin(BottomEdgeMassWindow[i][j]),qcd_Lead_CA12_mass[i][j]->FindBin(TopEdgeMassWindow[i][j]));
      //cout << "QCD fraction " << QCDfrac[i][j] << endl;
    }

    for (int j=0; j<nFineBins; j++){
      QCDfrac_finePt[i][j]=qcd_finePtBin_mass[i][j]->Integral(qcd_finePtBin_mass[i][j]->FindBin(BottomEdgeMassWindow_finePt[i][j]),qcd_finePtBin_mass[i][j]->FindBin(TopEdgeMassWindow_finePt[i][j]));
      //cout << "QCD fraction " << QCDfrac[i][j] << endl;

    }
    
  }
  
  
  makePtPlots();


  bool applyMassWindow = false; // should be read in from a command line option!
  //for (int ii = 0; ii < nAlgos-1; ii++)
  //{
  //int groomIdx = algoMap[ii];
  //int fileIdx = fileMap[ii]; // basically points to filteridx, reclusteridx, etc.
  makeMassWindowFile(applyMassWindow);//inputTree[fileIdx], inputTree[fileIdx]); // pass input file indices too?
  //  }
  

  return 0;
}



float DeltaR(float eta1,float phi1,float eta2,float phi2){
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


void runAlgorithm(TTree *inputTree, TTree *inputTree1, TString groomAlgo, int groomAlgoIndex)//, int fileidx)
{
  initializeVariables();
  //algoMap.push_back(groomAlgoIndex);
  //fileMap.push_back(fileidx);
  getMassHistograms(inputTree, inputTree1, groomAlgo, groomAlgoIndex);
  nAlgos++;
  deleteVectors();
}


void getMassHistograms(TTree *inputTree, TTree *inputTree1, TString groomAlgo, int groomAlgoIndex){
  
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
    
    //add weight to the QCD event
    if (RunNumber == 147913) weight*=1.6664E+03*1.9139E-03; 
    if (RunNumber == 147914) weight*=2.7646E+01*1.4296E-03; 
    if (RunNumber == 147915) weight*=3.0317E-01*5.5040E-03; 
    if (RunNumber == 147916) weight*=7.5078E-03*1.5252E-02;          
    if (RunNumber == 147917) weight*=1.3760E-03*7.6369E-02;  
    //cout << "weight after re-weight " << weight << endl;
    

    //NOW I WILL RECLUSTER MY JETS AND FEED THE APPROPRIATE VARIABLES 
    //FOR THE RECLUSTERING OPTIONS

    if (groomAlgoIndex>6){

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
    
    if (groomAlgoIndex == 0){
      
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
	qcd_Lead_CA12_pt[0]->Fill((*qcd_CA12_truth_pt)[chosenLeadTruthJetIndex]/1000, weight);
	qcd_Lead_CA12_mass[0][0]->Fill((*qcd_CA12_truth_mass)[chosenLeadTruthJetIndex]/1000, weight);
	
	if (((*qcd_CA12_truth_pt)[0]/1000)>1000.0){
	  qcd_Lead_CA12_mass[0][1]->Fill((*qcd_CA12_truth_mass)[chosenLeadTruthJetIndex]/1000, weight);
	}
	
	if (((*qcd_CA12_truth_pt)[0]/1000)>500.0 && ((*qcd_CA12_truth_pt)[0]/1000)<1000.0){
	  qcd_Lead_CA12_mass[0][2]->Fill((*qcd_CA12_truth_mass)[chosenLeadTruthJetIndex]/1000, weight);
	}
	if (((*qcd_CA12_truth_pt)[0]/1000)>1000.0 && ((*qcd_CA12_truth_pt)[0]/1000)<1500.0){
	  qcd_Lead_CA12_mass[0][3]->Fill((*qcd_CA12_truth_mass)[chosenLeadTruthJetIndex]/1000, weight);
	}
	if (((*qcd_CA12_truth_pt)[0]/1000)>1500.0 && ((*qcd_CA12_truth_pt)[0]/1000)<2000.0){
	  qcd_Lead_CA12_mass[0][4]->Fill((*qcd_CA12_truth_mass)[chosenLeadTruthJetIndex]/1000, weight);
	}
	if (((*qcd_CA12_truth_pt)[0]/1000)>2000.0){
	  qcd_Lead_CA12_mass[0][5]->Fill((*qcd_CA12_truth_mass)[chosenLeadTruthJetIndex]/1000, weight);
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

      qcd_Lead_CA12_pt[groomAlgoIndex]->Fill((*qcd_CA12_groomed_pt)[chosenLeadGroomedIndex]/1000, weight);
      qcd_Lead_CA12_mass[groomAlgoIndex][0]->Fill((*qcd_CA12_groomed_mass)[chosenLeadGroomedIndex]/1000, weight);
      
      if (((*qcd_CA12_truth_pt)[0]/1000)>1000.0){
	qcd_Lead_CA12_mass[groomAlgoIndex][1]->Fill((*qcd_CA12_groomed_mass)[chosenLeadGroomedIndex]/1000, weight);
      }
      
      if (((*qcd_CA12_truth_pt)[0]/1000)>500.0 && ((*qcd_CA12_truth_pt)[0]/1000)<1000.0){
	qcd_Lead_CA12_mass[groomAlgoIndex][2]->Fill((*qcd_CA12_groomed_mass)[chosenLeadGroomedIndex]/1000, weight);
      }
      if (((*qcd_CA12_truth_pt)[0]/1000)>1000.0 && ((*qcd_CA12_truth_pt)[0]/1000)<1500.0){
	qcd_Lead_CA12_mass[groomAlgoIndex][3]->Fill((*qcd_CA12_groomed_mass)[chosenLeadGroomedIndex]/1000, weight);
      }
      if (((*qcd_CA12_truth_pt)[0]/1000)>1500.0 && ((*qcd_CA12_truth_pt)[0]/1000)<2000.0){
	qcd_Lead_CA12_mass[groomAlgoIndex][4]->Fill((*qcd_CA12_groomed_mass)[chosenLeadGroomedIndex]/1000, weight);
      }
      if (((*qcd_CA12_truth_pt)[0]/1000)>2000.0){
	qcd_Lead_CA12_mass[groomAlgoIndex][5]->Fill((*qcd_CA12_groomed_mass)[chosenLeadGroomedIndex]/1000, weight);
      }
      
      //now for fine pT bins
      if (((*qcd_CA12_truth_pt)[0]/1000)>0.0 && ((*qcd_CA12_truth_pt)[0]/1000)<250.0 ){
	qcd_finePtBin_mass[groomAlgoIndex][0]->Fill((*qcd_CA12_groomed_mass)[chosenLeadGroomedIndex]/1000, weight);
      }

      for (int k=0; k<11; k++){
	if (k==10){
	  if (((*qcd_CA12_truth_pt)[0]/1000)>250.0*k){
	    qcd_finePtBin_mass[groomAlgoIndex][k]->Fill((*qcd_CA12_groomed_mass)[chosenLeadGroomedIndex]/1000, weight);
	  }
	}
	else	if (((*qcd_CA12_truth_pt)[0]/1000)>250.0*k && ((*qcd_CA12_truth_pt)[0]/1000)<250.0*(k+1)){
	  qcd_finePtBin_mass[groomAlgoIndex][k]->Fill((*qcd_CA12_groomed_mass)[chosenLeadGroomedIndex]/1000, weight);
	}
      }


      
    }
    
    //for pt reweighing purposes
    qcd_PtReweight[groomAlgoIndex]->Fill((*qcd_CA12_truth_pt)[0]/1000, weight);
    
  } // end loop over QCD events
  //end loop over QCD events
  
  //cout << "QCD efficiencies for " << groomAlgo << ": " << nEvt_0 << " " << nEvt_1 << " " << nEvt_2  << " " << nEvt_3  << " " << nEvt_4 << endl;
  
  
  //NOW I NEED A LOOP OVER Wp SAMPLES TO GET THEIR WEIGHT
  // I want to put this into another method

  int nentriesW=(int)inputTree1->GetEntries();
  cout<<"start Wp weight loop"<<endl;
  
  for (int ientry=0;ientry<nentriesW;ientry++) {
    inputTree1->GetEntry(ientry);

    Wp_PtReweight[groomAlgoIndex]->Fill((*Wp_CA12_truth_pt)[0]/1000);

  }


  //NOW I HAVE THE HISTOGRAMS TO DIVIDE AND GET THE WEIGHT HISTOGRAM

  for (int j=1; j<21; j++){
    //loop bins 1-20
    float weight=0.0;
    if (Wp_PtReweight[groomAlgoIndex]->GetBinContent(j)>0){
      weight = qcd_PtReweight[groomAlgoIndex]->GetBinContent(j)/Wp_PtReweight[groomAlgoIndex]->GetBinContent(j);
      pTweights[groomAlgoIndex]->SetBinContent(j,weight);
    }
    else pTweights[groomAlgoIndex]->SetBinContent(j,0.0);
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
	weight = pTweights[groomAlgoIndex]->GetBinContent(j+1); // nothing is being scaled here?!  This doesn't do anything
      }
    }
    
    

    //NOW I WILL RECLUSTER MY JETS AND FEED THE APPROPRIATE VARIABLES 
    //FOR THE RECLUSTERING OPTIONS

    if (groomAlgoIndex>6){

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
      
    if (groomAlgoIndex == 0){
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
	Wprime_Lead_CA12_pt[0]->Fill((*Wp_CA12_truth_pt)[chosenLeadTruthJetIndex]/1000, weight);
	Wprime_Lead_CA12_mass[0][0]->Fill((*Wp_CA12_truth_mass)[chosenLeadTruthJetIndex]/1000, weight);
	
	if (((*Wp_CA12_truth_pt)[0]/1000)>1000.0){
	  Wprime_Lead_CA12_mass[0][1]->Fill((*Wp_CA12_truth_mass)[chosenLeadTruthJetIndex]/1000, weight);
	}
	
	if (((*Wp_CA12_truth_pt)[0]/1000)>500.0 && ((*Wp_CA12_truth_pt)[0]/1000)<1000.0){
	  Wprime_Lead_CA12_mass[0][2]->Fill((*Wp_CA12_truth_mass)[chosenLeadTruthJetIndex]/1000, weight);
	}
	if (((*Wp_CA12_truth_pt)[0]/1000)>1000.0 && ((*Wp_CA12_truth_pt)[0]/1000)<1500.0){
	  Wprime_Lead_CA12_mass[0][3]->Fill((*Wp_CA12_truth_mass)[chosenLeadTruthJetIndex]/1000, weight);
	}
	if (((*Wp_CA12_truth_pt)[0]/1000)>1500.0 && ((*Wp_CA12_truth_pt)[0]/1000)<2000.0){
	  Wprime_Lead_CA12_mass[0][4]->Fill((*Wp_CA12_truth_mass)[chosenLeadTruthJetIndex]/1000, weight);
	}
	if (((*Wp_CA12_truth_pt)[0]/1000)>2000.0){
	  Wprime_Lead_CA12_mass[0][5]->Fill((*Wp_CA12_truth_mass)[chosenLeadTruthJetIndex]/1000, weight);
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

      Wprime_Lead_CA12_pt[groomAlgoIndex]->Fill((*Wp_CA12_groomed_pt)[chosenLeadGroomedIndex]/1000, weight);
      Wprime_Lead_CA12_mass[groomAlgoIndex][0]->Fill((*Wp_CA12_groomed_mass)[chosenLeadGroomedIndex]/1000, weight);
      
      if (((*Wp_CA12_truth_pt)[0]/1000)>1000.0){
	Wprime_Lead_CA12_mass[groomAlgoIndex][1]->Fill((*Wp_CA12_groomed_mass)[chosenLeadGroomedIndex]/1000, weight);
      }
      
      if (((*Wp_CA12_truth_pt)[0]/1000)>500.0 && ((*Wp_CA12_truth_pt)[0]/1000)<1000.0){
	Wprime_Lead_CA12_mass[groomAlgoIndex][2]->Fill((*Wp_CA12_groomed_mass)[chosenLeadGroomedIndex]/1000, weight);
      }
      if (((*Wp_CA12_truth_pt)[0]/1000)>1000.0 && ((*Wp_CA12_truth_pt)[0]/1000)<1500.0){
	Wprime_Lead_CA12_mass[groomAlgoIndex][3]->Fill((*Wp_CA12_groomed_mass)[chosenLeadGroomedIndex]/1000, weight);
      }
      if (((*Wp_CA12_truth_pt)[0]/1000)>1500.0 && ((*Wp_CA12_truth_pt)[0]/1000)<2000.0){
	Wprime_Lead_CA12_mass[groomAlgoIndex][4]->Fill((*Wp_CA12_groomed_mass)[chosenLeadGroomedIndex]/1000, weight);
      }
      if (((*Wp_CA12_truth_pt)[0]/1000)>2000.0){
	//if ((*Wp_CA12_groomed_pt)[chosenLeadGroomedIndex]/1000<500.00)
	//cout << (*Wp_CA12_groomed_pt)[chosenLeadGroomedIndex]/1000 << endl;
	Wprime_Lead_CA12_mass[groomAlgoIndex][5]->Fill((*Wp_CA12_groomed_mass)[chosenLeadGroomedIndex]/1000, weight);
      }

      //fine pt bins      
      for (int k=0; k<12; k++){
	if (((*Wp_CA12_truth_pt)[0]/1000)>250.0*k && ((*Wp_CA12_truth_pt)[0]/1000)<250.0*(k+1)){
	  Wprime_finePtBin_mass[groomAlgoIndex][k]->Fill((*Wp_CA12_groomed_mass)[chosenLeadGroomedIndex]/1000, weight);
	}
      }
      
    }
    //for rescaling check purposes
    
    Wprime_Lead_CA12_scaled_pt[groomAlgoIndex]->Fill((*Wp_CA12_truth_pt)[0]/1000, weight);
    /*
      Wprime_Lead_CA12_scaled_pt[8]->Fill((*Wp_CA12_truth_pt)[0]/1000, weight);
      Wprime_Lead_CA12_mass[8][0]->Fill((*Wp_CA12_truth_mass)[0]/1000, weight);
    */
  }
  //end loop over Wprime events
  
  cout << "Wprime efficiencies for " << groomAlgo << ": " << nEvt1_0 << " " << nEvt1_1 << " " << nEvt1_2  << " " << nEvt1_3  << " " << nEvt1_4 << endl;
  

}



void createHistos(){

  //create my histograms
  // make an array with histograms for all sorts of algorithms used
  int M = 100;
  for (int i=0; i<nAlgosMax-1; i++){
    //std::cout << nAlgos << std::endl;
    //for (int ii=0; ii<nAlgos; ii++){
    //int i = algoMap[ii];
    //std::cout << i << std::endl;
    for (int j=0; j<nFineBins; j++){
      qcd_finePtBin_mass[i][j] = new TH1F("qcd_finePtBin_mass_"+AlgoList[i]+finePtBins[j], "qcd_finePtBin_mass_"+AlgoList[i]+finePtBins[j], 150, 0.0, 300.0);
      //qcd_finePtBin_mass_curve[i][j] = new TH1D("qcd_finePtBin_mass_"+AlgoList[i]+finePtBins[j]+"_ROC","qcd_finePtBin_mass_"+AlgoList[i]+finePtBins[j]+"_ROC",M,0,1);
      Wprime_finePtBin_mass[i][j] = new TH1F("Wprime_finePtBin_mass_"+AlgoList[i]+finePtBins[j], "Wprime_finePtBin_mass_"+AlgoList[i]+finePtBins[j], 150, 0.0, 300.0);
      //finePtBin_mass_curve[i][j] = new TH1F("finePtBin_mass_"+AlgoList[i]+finePtBins[j]+"_ROC", "finePtBin_mass_"+AlgoList[i]+finePtBins[j]+"_ROC", M, 0, 1);
      //finePtBin_mass_curve[i][j] = new TGraph(M, 0, 1);
      qcd_finePtBin_mass[i][j]->Sumw2();
      Wprime_finePtBin_mass[i][j]->Sumw2();
    }


    for (int j=0; j<nPtBins; j++){
      qcd_Lead_CA12_mass[i][j] = new TH1F("qcd_Lead_CA12_mass_"+AlgoList[i]+pTbins[j], "qcd_Lead_CA12_mass_"+AlgoList[i]+pTbins[j], 300, 0.0, 300.0);
      //qcd_Lead_CA12_mass_curve[i][j] = new TH1F("qcd_Lead_CA12_mass_"+AlgoList[i]+pTbins[j]+"_ROC", "qcd_Lead_CA12_mass_"+AlgoList[i]+pTbins[j]+"_ROC", M, 0, 1);
      qcd_Lead_CA12_mass[i][j]->Sumw2();
      Wprime_Lead_CA12_mass[i][j] = new TH1F("Wprime_Lead_CA12_mass_"+AlgoList[i]+pTbins[j], "Wprime_Lead_CA12_mass_"+AlgoList[i]+pTbins[j], 300, 0.0, 300.0);
      //Lead_CA12_mass_curve[i][j] = new TH1F("Lead_CA12_mass_"+AlgoList[i]+pTbins[j]+"_ROC", "Lead_CA12_mass_"+AlgoList[i]+pTbins[j]+"_ROC",M, 0, 1);
      //Lead_CA12_mass_curve[i][j] = new TGraph(M, 0, 1);
      Wprime_Lead_CA12_mass[i][j]->Sumw2();
    }

    qcd_Lead_CA12_pt[i] = new TH1F("qcd_Lead_CA12_pt_"+AlgoList[i], "qcd_Lead_CA12_pt_"+AlgoList[i], 20, 0.0, 3500.0);
    Wprime_Lead_CA12_pt[i] = new TH1F("Wprime_Lead_CA12_pt_"+AlgoList[i], "Wprime_Lead_CA12_pt_"+AlgoList[i], 20, 0.0, 3500.0);
    Wprime_Lead_CA12_scaled_pt[i] = new TH1F("Wprime_Lead_CA12_scaled_pt_"+AlgoList[i], "Wprime_Lead_CA12_scaled_pt_"+AlgoList[i], 20, 0.0, 3500.0);
    qcd_PtReweight[i] = new TH1F(std::string("qcd_ptreweight_"+AlgoNames[i]).c_str(), std::string("qcd_ptreweight_"+AlgoNames[i]).c_str(), 20, 0.0, 3500.0);
    Wp_PtReweight[i] = new TH1F(std::string("Wp_ptreweight_"+AlgoNames[i]).c_str(), std::string("Wp_ptreweight_"+AlgoNames[i]).c_str(), 20, 0.0, 3500.0);
    pTweights[i] = new TH1F("pT_weights_"+AlgoList[i], "pt weights_"+AlgoList[i], 20, 0.0, 3500.0);

    // ROC curves
    //qcd_Lead_CA12_pt_curve[i] = new TH1F("qcd_Lead_CA12_pt_"+AlgoList[i]+"_ROC", "qcd_Lead_CA12_pt_"+AlgoList[i]+"_ROC", M, 0, 1);
    //Lead_CA12_pt_curve[i] = new TGraph(M,0,1);//TH1F("Lead_CA12_pt_"+AlgoList[i]+"_ROC", "Lead_CA12_pt_"+AlgoList[i]+"_ROC", M, 0, 1);
    //Lead_CA12_scaled_pt_curve[i] = new TH1F("Lead_CA12_scaled_pt_"+AlgoList[i]+"_ROC", "Lead_CA12_scaled_pt_"+AlgoList[i]+"_ROC", M, 0, 1);
    //Lead_CA12_scaled_pt_curve[i] = new TGraph(M, 0, 1);
    //qcd_PtReweight_curve[i] = new TH1F("qcd_ptreweight_"+AlgoList[i]+"_ROC", "qcd_ptreweight_"+AlgoList[i]+"_ROC", M, 0, 1);
    //PtReweight_curve[i] = new TH1F("ptreweight_"+AlgoList[i]+"_ROC", "ptreweight_"+AlgoList[i]+"_ROC", M, 0, 1);
    //PtReweight_curve[i] = new TGraph(M, 0, 1);
    //pTweights_curve[i] = new TH1F("pT_weights_"+AlgoList[i]+"_ROC", "pt weights_"+AlgoList[i]+"_ROC", M, 0, 1);
    //pTweights_curve[i] = new TGraph( M, 0, 1);
      
  }
  
  qcd_Lead_CA12_mass[nAlgosMax-1][0] = new TH1F("qcd_Lead_CA12_mass_"+AlgoList[nAlgosMax-1], "qcd_Lead_CA12_mass_"+AlgoList[nAlgosMax-1], 100, 0.0, 300.0);
  qcd_Lead_CA12_mass[nAlgosMax-1][0]->Sumw2();
  Wprime_Lead_CA12_mass[nAlgosMax-1][0] = new TH1F("Wprime_Lead_CA12_mass_"+AlgoList[nAlgosMax-1], "Wprime_Lead_CA12_mass_"+AlgoList[nAlgosMax-1], 100, 0.0, 300.0);
  Wprime_Lead_CA12_mass[nAlgosMax-1][0]->Sumw2();

  //qcd_Lead_CA12_mass_curve[nAlgosMax-1][0] = new TH1F("qcd_Lead_CA12_mass_"+AlgoList[nAlgosMax-1]+"_ROC", "qcd_Lead_CA12_mass_"+AlgoList[nAlgosMax-1]+"_ROC", M, 0, 1);
  //Lead_CA12_mass_curve[nAlgosMax-1][0] = new TH1F("Lead_CA12_mass_"+AlgoList[nAlgosMax-1]+"_ROC", "Lead_CA12_mass_"+AlgoList[nAlgosMax-1]+"_ROC", M, 0, 1);
  //Lead_CA12_mass_curve[nAlgosMax-1][0] = new TGraph( M, 0, 1);
  
  qcd_Lead_CA12_pt[nAlgosMax-1] = new TH1F("qcd_Lead_CA12_pt_"+AlgoList[nAlgosMax-1], "qcd_Lead_CA12_pt_"+AlgoList[nAlgosMax-1], 350, 0.0, 3500.0);
  Wprime_Lead_CA12_pt[nAlgosMax-1] = new TH1F("Wprime_Lead_CA12_pt_"+AlgoList[nAlgosMax-1], "Wprime_Lead_CA12_pt_"+AlgoList[nAlgosMax-1], 350, 0.0, 3500.0);
  Wprime_Lead_CA12_scaled_pt[nAlgosMax-1] = new TH1F("Wprime_Lead_CA12_scaled_pt_"+AlgoList[nAlgosMax-1], "Wprime_Lead_CA12_scaled_pt_"+AlgoList[nAlgosMax-1], 350, 0.0, 3500.0);

  //qcd_Lead_CA12_pt_curve[nAlgosMax-1] = new TH1F("qcd_Lead_CA12_pt_"+AlgoList[nAlgosMax-1]+"_ROC", "qcd_Lead_CA12_pt_"+AlgoList[nAlgosMax-1]+"_ROC", M, 0, 1);
  //Lead_CA12_pt_curve[nAlgosMax-1] = new TH1F("Lead_CA12_pt_"+AlgoList[nAlgosMax-1]+"_ROC", "Lead_CA12_pt_"+AlgoList[nAlgosMax-1]+"_ROC", M, 0, 1);
  //Lead_CA12_pt_curve[nAlgosMax-1] = new TGraph( M, 0, 1);
  //Lead_CA12_scaled_pt_curve[nAlgosMax-1] = new TH1F("Lead_CA12_scaled_pt_"+AlgoList[nAlgosMax-1]+"_ROC", "Lead_CA12_scaled_pt_"+AlgoList[nAlgosMax-1]+"_ROC", M, 0, 1);
  //Lead_CA12_scaled_pt_curve[nAlgosMax-1] = new TGraph( M, 0, 1);

}


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
/*
  void getNormSherpaW(TString InputFile, unsigned long & evnum, double & weight)
  {
  //Setup TChain that will access the proper metadata location in the TFile
  TChain *chCutFlow = new TChain("physicsMeta/CutFlowTree");

  //Load the desired TFile into the TChain
  //TString InputFile="path/to/file.root"
  chCutFlow->Add(InputFile);
    
  //Using the loaded TChain from the input file, you access the sum of the events and the sum of weighted events before the filter
  if(chCutFlow!=0 && chCutFlow->GetEntries() != 0){
  std::vector<unsigned long> *nAcceptedEvents = new std::vector<unsigned long>;
  std::vector<double> *nWeightedAcceptedEvents = new std::vector<double>;
  std::vector<std::string> *name = new std::vector<std::string>;
  chCutFlow->SetBranchAddress("nAcceptedEvents", &nAcceptedEvents);
  chCutFlow->SetBranchAddress("nWeightedAcceptedEvents", &nWeightedAcceptedEvents);
  chCutFlow->SetBranchAddress("name", &name);
  for(int i=0;i<chCutFlow->GetEntries();i++){
  chCutFlow->GetEntry(i);
  for(unsigned int j=0; j<name->size(); j++){
  if((*name)[j] == "AllExecutedEvents"){
  evnum+=(*nAcceptedEvents)[j];
  weight+=(*nWeightedAcceptedEvents)[j];
  break;
  }
  }
  }
  }
  else{
  evnum=-1;
  weight=-1;
  }
  cout<<"NTuple Count: "<<evnum<<"  "<<weight<<endl;

  }
*/
void getBranches(TTree *inputTree, TTree *inputTree1, TString groomAlgo, int groomAlgoIndex){

  //if not reclustering
  if (groomAlgoIndex<7){
    
    inputTree->SetBranchAddress("mc_channel_number", &qcd_mc_channel_number); // segfaults for groomalgo = <incomplete type>, groomalgoindex = 3, TopoSplitFilteredMu100SmallR30YCut4
    inputTree->SetBranchAddress("mc_event_weight", &qcd_mc_event_weight);
    inputTree->SetBranchAddress(std::string("jet_CamKt12Truth_pt").c_str(), &qcd_CA12_truth_pt);
    inputTree->SetBranchAddress(std::string("jet_CamKt12Truth_eta").c_str(), &qcd_CA12_truth_eta);
    inputTree->SetBranchAddress(std::string("jet_CamKt12Truth_phi").c_str(), &qcd_CA12_truth_phi);
    inputTree->SetBranchAddress(std::string("jet_CamKt12Truth_m").c_str(), &qcd_CA12_truth_mass);
    inputTree->SetBranchAddress(std::string("jet_CamKt12Truth_E").c_str(), &qcd_CA12_truth_E);
    
    inputTree->SetBranchAddress(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+"LCTopo_pt").c_str(), &qcd_CA12_topo_pt);
    inputTree->SetBranchAddress(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+"LCTopo_eta").c_str(), &qcd_CA12_topo_eta);
    inputTree->SetBranchAddress(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+"LCTopo_phi").c_str(), &qcd_CA12_topo_phi);
    inputTree->SetBranchAddress(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+"LCTopo_m").c_str(), &qcd_CA12_topo_mass);
    inputTree->SetBranchAddress(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+"LCTopo_E").c_str(), &qcd_CA12_topo_E);
    inputTree->SetBranchAddress(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+"LCTopo_emfrac").c_str(), &qcd_CA12_topo_emfrac);
    
    inputTree->SetBranchAddress(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+"LC"+groomAlgo+"_pt").c_str(), &qcd_CA12_groomed_pt);
    inputTree->SetBranchAddress(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+"LC"+groomAlgo+"_eta").c_str(), &qcd_CA12_groomed_eta);
    inputTree->SetBranchAddress(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+"LC"+groomAlgo+"_phi").c_str(), &qcd_CA12_groomed_phi);
    inputTree->SetBranchAddress(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+"LC"+groomAlgo+"_m").c_str(), &qcd_CA12_groomed_mass);
    inputTree->SetBranchAddress(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+"LC"+groomAlgo+"_E").c_str(), &qcd_CA12_groomed_E);
    inputTree->SetBranchAddress(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+"LC"+groomAlgo+"_emfrac").c_str(), &qcd_CA12_groomed_emfrac);
    
    inputTree1->SetBranchAddress("mc_channel_number", &Wp_mc_channel_number);
    inputTree1->SetBranchAddress("mc_event_weight", &Wp_mc_event_weight);
    inputTree1->SetBranchAddress(std::string("jet_CamKt12Truth_pt").c_str(), &Wp_CA12_truth_pt);
    inputTree1->SetBranchAddress(std::string("jet_CamKt12Truth_eta").c_str(), &Wp_CA12_truth_eta);
    inputTree1->SetBranchAddress(std::string("jet_CamKt12Truth_phi").c_str(), &Wp_CA12_truth_phi);
    inputTree1->SetBranchAddress(std::string("jet_CamKt12Truth_m").c_str(), &Wp_CA12_truth_mass);
    inputTree1->SetBranchAddress(std::string("jet_CamKt12Truth_E").c_str(), &Wp_CA12_truth_E);
    
    inputTree1->SetBranchAddress(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+"LCTopo_pt").c_str(), &Wp_CA12_topo_pt);
    inputTree1->SetBranchAddress(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+"LCTopo_eta").c_str(), &Wp_CA12_topo_eta);
    inputTree1->SetBranchAddress(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+"LCTopo_phi").c_str(), &Wp_CA12_topo_phi);
    inputTree1->SetBranchAddress(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+"LCTopo_m").c_str(), &Wp_CA12_topo_mass);
    inputTree1->SetBranchAddress(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+"LCTopo_E").c_str(), &Wp_CA12_topo_E);
    inputTree1->SetBranchAddress(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+"LCTopo_emfrac").c_str(), &Wp_CA12_topo_emfrac);
    
    inputTree1->SetBranchAddress(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+"LC"+groomAlgo+"_pt").c_str(), &Wp_CA12_groomed_pt);
    inputTree1->SetBranchAddress(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+"LC"+groomAlgo+"_eta").c_str(), &Wp_CA12_groomed_eta);
    inputTree1->SetBranchAddress(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+"LC"+groomAlgo+"_phi").c_str(), &Wp_CA12_groomed_phi);
    inputTree1->SetBranchAddress(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+"LC"+groomAlgo+"_m").c_str(), &Wp_CA12_groomed_mass);
    inputTree1->SetBranchAddress(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+"LC"+groomAlgo+"_E").c_str(), &Wp_CA12_groomed_E);
    inputTree1->SetBranchAddress(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+"LC"+groomAlgo+"_emfrac").c_str(), &Wp_CA12_groomed_emfrac);

    ///////////////////////////////////////////////
    /*
    inputTree->SetBranchAddress("mc_channel_number", &qcd_mc_channel_number);
    inputTree->SetBranchAddress("mc_event_weight", &qcd_mc_event_weight);
    inputTree->SetBranchAddress("jet_AntiKt10Truth_pt", &qcd_CA12_truth_pt);
    inputTree->SetBranchAddress("jet_AntiKt10Truth_eta", &qcd_CA12_truth_eta);
    inputTree->SetBranchAddress("jet_AntiKt10Truth_phi", &qcd_CA12_truth_phi);
    inputTree->SetBranchAddress("jet_AntiKt10Truth_m", &qcd_CA12_truth_mass);
    inputTree->SetBranchAddress("jet_AntiKt10Truth_E", &qcd_CA12_truth_E);

    inputTree->SetBranchAddress("jet_AntiKt10"+groomAlgo+"_pt", &qcd_CA12_groomed_pt);
    inputTree->SetBranchAddress("jet_AntiKt10"+groomAlgo+"_eta", &qcd_CA12_groomed_eta);
    inputTree->SetBranchAddress("jet_AntiKt10"+groomAlgo+"_phi", &qcd_CA12_groomed_phi);
    inputTree->SetBranchAddress("jet_AntiKt10"+groomAlgo+"_m", &qcd_CA12_groomed_mass);
    inputTree->SetBranchAddress("jet_AntiKt10"+groomAlgo+"_E", &qcd_CA12_groomed_E);
    inputTree->SetBranchAddress("jet_AntiKt10"+groomAlgo+"_emfrac", &qcd_CA12_groomed_emfrac);

    inputTree1->SetBranchAddress("mc_channel_number", &Wp_mc_channel_number);
    inputTree1->SetBranchAddress("mc_event_weight", &Wp_mc_event_weight);
    inputTree1->SetBranchAddress("jet_AntiKt10Truth_pt", &Wp_CA12_truth_pt);
    inputTree1->SetBranchAddress("jet_AntiKt10Truth_eta", &Wp_CA12_truth_eta);
    inputTree1->SetBranchAddress("jet_AntiKt10Truth_phi", &Wp_CA12_truth_phi);
    inputTree1->SetBranchAddress("jet_AntiKt10Truth_m", &Wp_CA12_truth_mass);
    inputTree1->SetBranchAddress("jet_AntiKt10Truth_E", &Wp_CA12_truth_E);

    inputTree1->SetBranchAddress("jet_AntiKt10LCTopo_pt", &Wp_CA12_topo_pt);
    inputTree1->SetBranchAddress("jet_AntiKt10LCTopo_eta", &Wp_CA12_topo_eta);
    inputTree1->SetBranchAddress("jet_AntiKt10LCTopo_phi", &Wp_CA12_topo_phi);
    inputTree1->SetBranchAddress("jet_AntiKt10LCTopo_m", &Wp_CA12_topo_mass);
    inputTree1->SetBranchAddress("jet_AntiKt10LCTopo_E", &Wp_CA12_topo_E);
    inputTree1->SetBranchAddress("jet_AntiKt10LCTopo_emfrac", &Wp_CA12_topo_emfrac);
    */

  }
 



  //if reclustering algorithm
  else if (groomAlgoIndex>6){

    /*
    inputTree->SetBranchAddress("mc_channel_number", &qcd_mc_channel_number);
    inputTree->SetBranchAddress("mc_event_weight", &qcd_mc_event_weight);
    inputTree->SetBranchAddress("jet_CamKt12Truth_pt", &qcd_CA12_truth_pt);
    inputTree->SetBranchAddress("jet_CamKt12Truth_eta", &qcd_CA12_truth_eta);
    inputTree->SetBranchAddress("jet_CamKt12Truth_phi", &qcd_CA12_truth_phi);
    inputTree->SetBranchAddress("jet_CamKt12Truth_m", &qcd_CA12_truth_mass);
    inputTree->SetBranchAddress("jet_CamKt12Truth_E", &qcd_CA12_truth_E);
    */
    inputTree->SetBranchAddress("mc_channel_number", &qcd_mc_channel_number);
    inputTree->SetBranchAddress("mc_event_weight", &qcd_mc_event_weight);
    inputTree->SetBranchAddress(std::string("jet_CamKt12Truth_pt").c_str(), &qcd_CA12_truth_pt);
    inputTree->SetBranchAddress(std::string("jet_CamKt12Truth_eta").c_str(), &qcd_CA12_truth_eta);
    inputTree->SetBranchAddress(std::string("jet_CamKt12Truth_phi").c_str(), &qcd_CA12_truth_phi);
    inputTree->SetBranchAddress(std::string("jet_CamKt12Truth_m").c_str(), &qcd_CA12_truth_mass);
    inputTree->SetBranchAddress(std::string("jet_CamKt12Truth_E").c_str(), &qcd_CA12_truth_E);
    
    //AntiKt2LCTopo

    /*inputTree->SetBranchAddress("jet_"+groomAlgo+"_pt", &qcd_CA12_topo_pt);
    inputTree->SetBranchAddress("jet_"+groomAlgo+"_eta", &qcd_CA12_topo_eta);
    inputTree->SetBranchAddress("jet_"+groomAlgo+"_phi", &qcd_CA12_topo_phi);
    inputTree->SetBranchAddress("jet_"+groomAlgo+"_m", &qcd_CA12_topo_mass);
    inputTree->SetBranchAddress("jet_"+groomAlgo+"_E", &qcd_CA12_topo_E);
    inputTree->SetBranchAddress("jet_"+groomAlgo+"_emfrac", &qcd_CA12_topo_emfrac);
    */
    inputTree->SetBranchAddress(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+""+groomAlgo+"_pt").c_str(), &qcd_CA12_groomed_pt);
      inputTree->SetBranchAddress(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+""+groomAlgo+"_eta").c_str(), &qcd_CA12_groomed_eta);
      inputTree->SetBranchAddress(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+""+groomAlgo+"_phi").c_str(), &qcd_CA12_groomed_phi);
      inputTree->SetBranchAddress(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+""+groomAlgo+"_m").c_str(), &qcd_CA12_groomed_mass);
      inputTree->SetBranchAddress(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+""+groomAlgo+"_E").c_str(), &qcd_CA12_groomed_E);
      inputTree->SetBranchAddress(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+""+groomAlgo+"_emfrac").c_str(), &qcd_CA12_groomed_emfrac);
    
 
      /*inputTree1->SetBranchAddress("mc_channel_number", &Wp_mc_channel_number);
    inputTree1->SetBranchAddress("mc_event_weight", &Wp_mc_event_weight);
    inputTree1->SetBranchAddress("jet_CamKt12Truth_pt", &Wp_CA12_truth_pt);
    inputTree1->SetBranchAddress("jet_CamKt12Truth_eta", &Wp_CA12_truth_eta);
    inputTree1->SetBranchAddress("jet_CamKt12Truth_phi", &Wp_CA12_truth_phi);
    inputTree1->SetBranchAddress("jet_CamKt12Truth_m", &Wp_CA12_truth_mass);
    inputTree1->SetBranchAddress("jet_CamKt12Truth_E", &Wp_CA12_truth_E);*/
    
    
      inputTree1->SetBranchAddress("mc_channel_number", &Wp_mc_channel_number);
      //inputTree1->SetBranchAddress("mc_event_weight", &Wp_mc_event_weight);
      inputTree1->SetBranchAddress(std::string("jet_CamKt12Truth_pt").c_str(), &Wp_CA12_truth_pt);
      inputTree1->SetBranchAddress(std::string("jet_CamKt12Truth_eta").c_str(), &Wp_CA12_truth_eta);
      inputTree1->SetBranchAddress(std::string("jet_CamKt12Truth_phi").c_str(), &Wp_CA12_truth_phi);
      inputTree1->SetBranchAddress(std::string("jet_CamKt12Truth_m").c_str(), &Wp_CA12_truth_mass);
      inputTree1->SetBranchAddress(std::string("jet_CamKt12Truth_E").c_str(), &Wp_CA12_truth_E);
      
    inputTree1->SetBranchAddress(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+"LCTopo_pt").c_str(), &Wp_CA12_topo_pt);
    inputTree1->SetBranchAddress(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+"LCTopo_eta").c_str(), &Wp_CA12_topo_eta);
    inputTree1->SetBranchAddress(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+"LCTopo_phi").c_str(), &Wp_CA12_topo_phi);
    inputTree1->SetBranchAddress(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+"LCTopo_m").c_str(), &Wp_CA12_topo_mass);
    inputTree1->SetBranchAddress(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+"LCTopo_E").c_str(), &Wp_CA12_topo_E);
    inputTree1->SetBranchAddress(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+"LCTopo_emfrac").c_str(), &Wp_CA12_topo_emfrac);
    
      inputTree1->SetBranchAddress(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+""+groomAlgo+"_pt").c_str(), &Wp_CA12_groomed_pt);
      inputTree1->SetBranchAddress(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+""+groomAlgo+"_eta").c_str(), &Wp_CA12_groomed_eta);
      inputTree1->SetBranchAddress(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+""+groomAlgo+"_phi").c_str(), &Wp_CA12_groomed_phi);
      inputTree1->SetBranchAddress(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+""+groomAlgo+"_m").c_str(), &Wp_CA12_groomed_mass);
      inputTree1->SetBranchAddress(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+""+groomAlgo+"_E").c_str(), &Wp_CA12_groomed_E);
      inputTree1->SetBranchAddress(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+""+groomAlgo+"_emfrac").c_str(), &Wp_CA12_groomed_emfrac);
    

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



void makePlots(){
  ////////PLOT THINGS -> This should go into it's own class or at the very least it's own method, TODO

  // CANVAS FOR THE MASS PLOTS
  
  //for (int i=0; i<nAlgos-1; i++){//-2; i++){
  for (int ii=0; ii<nAlgos; ii++){//-2; i++){
    //std::cout << "ii: " << ii << std::endl;
    int i = algoMap[ii];
    //std::cout << "i: " << i << std::endl;
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
      //Wprime_Lead_CA12_mass[i][j]->DrawCopy("E");
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
      p2.DrawLatex(0.20,0.82,AlgoListN[i]);

    
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

      c1[i][j]->SaveAs("jetMass_"+AlgoList[i]+pTbins[j]+".eps");
      
    }
  }
  
  /////////END OF PLOTTING THINGS
  
  ///// Fine pT bins plots
  
  //for (int i=0; i<nAlgos-1; i++){//-2; i++){
  for (int ii=0; ii<nAlgos-1; ii++){//-2; i++){
    int i = algoMap[ii];
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
      //Wprime_Lead_CA12_mass[i][j]->DrawCopy("E");
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
      p2.DrawLatex(0.20,0.82,AlgoListN[i]);

    
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

      c3[i][j]->SaveAs("jetMass_fineBins_"+AlgoList[i]+finePtBins[j]+".eps");

    }
  }
} // end makePlots()


void makePtPlots(){

  // Now I want to make plots which enclose all those information, for each pT bin

  for (int i=0; i<nPtBins; i++){
    hMassLow[i] = new TH1F("MassLow_"+pTbins[i],"MassLow_"+pTbins[i],nAlgos-2, 0, nAlgos-2);
    hMassHigh[i] = new TH1F("MassHigh_"+pTbins[i],"MassHigh_"+pTbins[i],nAlgos-2, 0, nAlgos-2);
    hMPV[i] = new TH1F("MPV_"+pTbins[i],"MPV_"+pTbins[i],nAlgos-2, 0, nAlgos-2);
    hWmass[i] = new TH1F("Wmass_"+pTbins[i],"Wmass_"+pTbins[i],nAlgos-2, 0, nAlgos-2);
    hQCDeff[i] = new TH1F("QCDeff_"+pTbins[i],"QCDeff_"+pTbins[i], nAlgos-2, 0, nAlgos-2);

  }

  // for (int i=0; i<nAlgos-2; i++){
  //for now only the algos I already have: 0,1,2,3,4,5,6,7,8,9 
  for (int ii=0; ii<nAlgos; ii++){  //-2; ii++){  
    int i = algoMap[ii];
    for (int j=0; j<nPtBins; j++){
      hMassLow[j]->SetBinContent(i+1,BottomEdgeMassWindow[i][j]);
      hMassHigh[j]->SetBinContent(i+1,TopEdgeMassWindow[i][j]);
      hMPV[j]->SetBinContent(i+1,myMPV[i][j]);
      hQCDeff[j]->SetBinContent(i+1,QCDfrac[i][j]);
      hWmass[j]->SetBinContent(i+1,80.0);
    }
  }
  
  //NOW PLOT WINDOW WIDTH VS PT
  //Christo's plot
  for (int ii=0; ii<nAlgos; ii++){//-2; ii++){
    int i = algoMap[ii];
    windowsVsPt[i] = new TH1F("windowVsPt_"+AlgoList[i], "windowVsPt_"+AlgoList[i], 10, 0.0, 2500.);
    windowsVsPt[i]->Sumw2();
    windowsVsPt[i]->SetBinContent(1,WidthMassWindow_finePt[i][0]);
    windowsVsPt[i]->SetBinContent(2,WidthMassWindow_finePt[i][1]);
    windowsVsPt[i]->SetBinContent(3,WidthMassWindow_finePt[i][2]);
    windowsVsPt[i]->SetBinContent(4,WidthMassWindow_finePt[i][3]);
    windowsVsPt[i]->SetBinContent(5,WidthMassWindow_finePt[i][4]);
    windowsVsPt[i]->SetBinContent(6,WidthMassWindow_finePt[i][5]);
    windowsVsPt[i]->SetBinContent(7,WidthMassWindow_finePt[i][6]);
    windowsVsPt[i]->SetBinContent(8,WidthMassWindow_finePt[i][7]);
    windowsVsPt[i]->SetBinContent(9,WidthMassWindow_finePt[i][8]);
    windowsVsPt[i]->SetBinContent(10,WidthMassWindow_finePt[i][9]);
    //windowsVsPt[i]->SetBinContent(11,WidthMassWindow_finePt[i][10]);
    //windowsVsPt[i]->SetBinContent(12,WidthMassWindow_finePt[i][11]);
    
    hQCDeff_finePt[i] = new TH1F("QCD_fineBin_"+AlgoList[i],"QCD_fineBin_"+AlgoList[i],10,0.0, 2500.);
    
    for (int k=1; k<nFineBins-1; k++){
      
      hQCDeff_finePt[i]->SetBinContent(k, QCDfrac_finePt[i][k-1]);
      
    }

    for (int k=0; k<nFineBins; k++){
      windowsVsPt[i]->SetBinError(k+1,2.0);

    }
  }
  

  TCanvas *wVsPt[nAlgos]; // was -2
  TPad *pad11[nAlgos];
  TPad *pad22[nAlgos];
	     
  for (int ii=0; ii<nAlgos; ii++){//ii = 1; ...... -2; ii++){
    int i = algoMap[ii];
    wVsPt[i] = new TCanvas("wVsPt", "wVsPt", 700, 550);
    
    pad11[i] = new TPad("pad1","pad1",0,0.4,1,1);
    pad11[i]->SetBottomMargin(0.009);
    pad11[i]->Draw();
    pad22[i] = new TPad("pad2","pad2",0,0,1,0.4);
    pad22[i]->SetTopMargin(0.009);
    pad22[i]->SetBottomMargin(0.25);
    pad22[i]->Draw();

    pad11[i]->cd();

    windowsVsPt[i]->GetYaxis()->SetTitle("Mass width (68% W eff) [GeV]");  
    windowsVsPt[i]->GetYaxis()->SetTitleSize(0.06);  
    windowsVsPt[i]->GetYaxis()->SetTitleOffset(0.9);

    windowsVsPt[i]->GetXaxis()->SetTitle("Jet CA12 Truth p_{T} [GeV]");  
    windowsVsPt[i]->SetMaximum(windowsVsPt[i]->GetMaximum()*1.3);
    windowsVsPt[i]->Draw("P");

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
    p2.DrawLatex(0.20,0.80,AlgoListN[i]);
    

    

    pad22[i]->cd();
    hQCDeff_finePt[i]->SetLineColor(2);
    //hQCDeff_finePt[i]->SetMaximum(1.1*hQCDeff_finePt[i]->GetMaximum());
    hQCDeff_finePt[i]->SetMinimum(0.05);
    hQCDeff_finePt[i]->SetMaximum(0.5);
    hQCDeff_finePt[i]->Draw();


    hQCDeff_finePt[i]->GetYaxis()->SetTitle("QCD fraction");
    hQCDeff_finePt[i]->GetXaxis()->SetTitle("Jet CA12 Truth p_{T} [GeV]");  
    hQCDeff_finePt[i]->GetXaxis()->SetTitleSize(0.09);
    hQCDeff_finePt[i]->GetXaxis()->SetLabelSize(0.08);
    hQCDeff_finePt[i]->GetYaxis()->SetLabelSize(0.08);
    hQCDeff_finePt[i]->GetYaxis()->SetTitleSize(0.09);
    hQCDeff_finePt[i]->GetYaxis()->SetTitleOffset(0.55);
    

    wVsPt[i]->SaveAs("windowsVsPt_"+AlgoList[i]+".eps");
    
  }

  
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
      hQCDeff[i]->GetXaxis()->SetBinLabel(k,binLabel[k-1]);
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
  
  for (int ii=0; ii<nAlgos; ii++){//-2; ii++){
    int i = algoMap[ii];
    windowsVsPt[i]->Write();
  }

  for (int ii=0; ii<nAlgos; ii++){//-1; ii++){
    int i = algoMap[ii];
    
    qcd_Lead_CA12_pt[i]->Write();
    Wprime_Lead_CA12_pt[i]->Write();
    makeROC(type, Wprime_Lead_CA12_pt[i],qcd_Lead_CA12_pt[i],(*Lead_CA12_pt_curve)[i], Wprime_Lead_CA12_pt[i]->GetName(), true);
    //Lead_CA12_pt_curve[i]->Write();

    pTweights[i]->Write();
    qcd_PtReweight[i]->Write();
    Wp_PtReweight[i]->Write();
    makeROC(type, Wp_PtReweight[i], qcd_PtReweight[i], (*PtReweight_curve)[i], "PTReweight", true);
    //PtReweight_curve[i]->Write();
    Wprime_Lead_CA12_scaled_pt[i]->Write();
    makeROC(type, Wprime_Lead_CA12_scaled_pt[i], qcd_Lead_CA12_pt[i], (*Lead_CA12_scaled_pt_curve)[i], Wprime_Lead_CA12_scaled_pt[i]->GetName(), true);// Is his right?!
    //Lead_CA12_scaled_pt_curve[i]->Write();


    for (int j=0; j<nPtBins; j++){
      qcd_Lead_CA12_mass[i][j]->Write();
      Wprime_Lead_CA12_mass[i][j]->Write();
      makeROC(type, Wprime_Lead_CA12_mass[i][j], qcd_Lead_CA12_mass[i][j], (*Lead_CA12_mass_curve)[i][j], Wprime_Lead_CA12_mass[i][j]->GetName(), true);
      //Lead_CA12_mass_curve[i][j]->Write();
    }

    for (int j=0; j<nFineBins; j++){
      qcd_finePtBin_mass[i][j]->Write();
      Wprime_finePtBin_mass[i][j]->Write();
      makeROC(type, Wprime_finePtBin_mass[i][j], qcd_finePtBin_mass[i][j], (*finePtBin_mass_curve)[i][j], Wprime_finePtBin_mass[i][j]->GetName(), true);
      //finePtBin_mass_curve[i][j]->Write();
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


void makeMassWindowFile(bool applyMassWindow)
{
  std::cout << "starting make mass window output files " << std::endl;
  double mass_max = 0.0;
  double mass_min = 0.0;

  bool signal = false;
  // need branches for specific algorithms because otherwise we end up with a million branches that we don't need
  vector<std::string> branches;
  std::cout << "number of algorithms" << nAlgos<< std::endl;
  for (int ii=0; ii<nAlgos; ii++){

    int i = algoMap[ii];
    std::cout << " doing mass eff plots for " << AlgoNames[i] << std::endl;
    int groomAlgoIndex = i;
    std::string prefix = "";
    branches = getListOfJetBranches(AlgoNames[i]);
    for (int j=0; j<nPtBins; j++){
      signal = false;
      for (int k = 0; k < 2; k++)
	{
	  signal = k == 0 ? false : true;
	  initVectors();
	  int fileIdx = signal ? fileMap[ii]+1 : fileMap[ii];
	  //std::cout << "fileIdx " << fileIdx << endl;
	  TTree * intree = (TTree*) inputFile[fileIdx]->Get("physics");
	  intree->SetBranchStatus("*",0);

	  // turn on the branches we're interested it
	  for (std::vector<string>::iterator it = branches.begin(); it != branches.end(); it++)
	    {
	      intree->SetBranchStatus((*it).c_str(),1);
	      if(!intree->GetListOfBranches()->FindObject((*it).c_str())) {
		std::cout << "could not find branch: " << (*it) << std::endl;
	      }
	    }
	  setMassBranch(intree, AlgoNames[i], i);

	  std::stringstream ss; // store the name of the output file and include the i and j indices!
	  std::string bkg = signal ? "sig": "bkg";
	  ss << AlgoNames[i] << "_" << i << "_" << pTbins[j] << "_" << bkg << ".root";
	  boost::filesystem::path dir(AlgoNames[i]);
	  boost::filesystem::create_directory(dir);
	  TFile * outfile = new TFile(std::string(AlgoNames[i]+"/"+ss.str()).c_str(),"RECREATE");	  
	  TTree * outTree;


	  addJets(intree, AlgoNames[i], signal, i);//, addJetIndices); //set all of the branches for the output tree for the jets
	  outTree = intree->CloneTree(0);
	  mass_max = TopEdgeMassWindow[i][j];
	  mass_min = BottomEdgeMassWindow[i][j];

	  long entries = (long)intree->GetEntries();

	  // these numbers are chosen somewhat arb - they come from the settings I use
	  // in the config files for the plotter() code... such bad coding :(
	  pt_reweight = new TH1F(std::string("pt_reweight").c_str(),std::string("pt_reweight_").c_str(), 20, 0, 1200);
	  
	  double mass = 0;
	  NEvents = entries;
	  NEvents_weighted.clear();
	  //NEvents_weighted = 0;
	  for (long n = 0; n < entries; n++)
	    {
	      intree->GetEntry(n);
	      if (NEvents_weighted.find(mc_channel_number) != NEvents_weighted.end())
		NEvents_weighted[mc_channel_number] += mc_event_weight;
	      else
		NEvents_weighted[mc_channel_number] = mc_event_weight;
	      // what about reclustered jets?! argggg
	      int chosenLeadTruthJetIndex=-99;
	      int chosenLeadTopoJetIndex=-99;
	      setSelectionVectors(signal, AlgoNames[i]);
	      if ((*jet_pt_truth)[0] / 1000.0 < 100) 
		continue;
	      if (groomAlgoIndex != 0) // check the pt is in the correct bin
		{
		  if ((*jet_pt_truth)[0]/1000.0 < ptrange[j].first || (*jet_pt_truth)[0]/1000.0 > ptrange[j].second)
		    continue;
		}
	      if (groomAlgoIndex == 0){
		//check which is my CA12 ungroomed reco jet with EMfrac<0.99 and |eta|<1.2
		//loop over all topo jets and get the leading with cuts
		bool hasTopoJet=false;
		
		for (int jet_i=0; jet_i<(*jet_pt_topo).size(); jet_i++){
		  if (!hasTopoJet && (*jet_emfrac_topo)[jet_i]<0.99 && fabs((*jet_eta_topo)[jet_i])<1.2) {
		    chosenLeadTopoJetIndex=jet_i;
		    hasTopoJet=true;
		  } 
		} // end loop over jet_pt_topo
		
		//now match my truth jet with the chosen topo jet
		
		if (hasTopoJet){
		  for (int jet_i=0; jet_i<(*jet_pt_truth).size(); jet_i++){
		    if (chosenLeadTruthJetIndex<0 && DeltaR((*jet_eta_topo)[chosenLeadTopoJetIndex],(*jet_phi_topo)[chosenLeadTopoJetIndex],(*jet_eta_truth)[jet_i],(*jet_phi_truth)[jet_i])<0.9){ 
		      chosenLeadTruthJetIndex=jet_i;
		    }	  
		  }	// end loop over jet_pt_truth
		} // end if(hasTopoJet)
		//if ((*jet_pt_truth)[chosenLeadTruthJetIndex]/1000.0 > ptrange[j].first || (*jet_pt_truth)[chosenLeadTruthJetIndex]/1000.0 < ptrange[j].second)
		//continue;
	      } // end if(groomAlgoIndex==0)
	      
	      chosenLeadTruthJetIndex = groomAlgoIndex == 0 ? chosenLeadTruthJetIndex : 0;
	      
	      //Now I have which events to make my pt reweight with, and to match to, etc
	      int chosenLeadGroomedIndex=-99;
	      for (int jet_i=0; jet_i<(*jet_pt_groomed).size(); jet_i++){
		
		if (chosenLeadTruthJetIndex>=0 && chosenLeadGroomedIndex<0 && DeltaR((*jet_eta_truth)[chosenLeadTruthJetIndex],(*jet_phi_truth)[chosenLeadTruthJetIndex],(*jet_eta_groomed)[jet_i],(*jet_phi_groomed)[jet_i])<0.9 && (*jet_emfrac_groomed)[jet_i]<0.99 && fabs((*jet_eta_groomed)[jet_i])<1.2){
		  chosenLeadGroomedIndex=jet_i;
		}     
	      } // end loop over jet_pt_groomed
	      
	      if (chosenLeadTruthJetIndex < 0 || chosenLeadGroomedIndex == -99) // failed selection
		  continue;
	      

	      
	    leadGroomedIndex = chosenLeadGroomedIndex;
	    leadTruthIndex = chosenLeadTruthJetIndex;
	    leadTopoIndex = chosenLeadTopoJetIndex;
	    mass = signal ? (*signal_m_vec[2])[leadGroomedIndex]/1000.0 : (*bkg_m_vec[2])[leadGroomedIndex]/1000.0;

	    if (applyMassWindow && (mass > mass_max && mass < mass_min))
	      {
		continue;
	      }
	    outTree->Fill();
	    pt_reweight->Fill((*jet_pt_truth)[chosenLeadTruthJetIndex]/1000.0);

	    }
	  // write the rweight th1f things to the outfile...
	  // calculate the new reweight with a new histogram!

	  outTree->GetCurrentFile()->Write();
	  pt_reweight->Write();
	  if (!signal)
	    {
	      qcd_PtReweight[i]->Write();
	    }
	  else
	    {
	      Wp_PtReweight[i]->Write();
	    }
	  outTree->GetCurrentFile()->Print();
	  outTree->GetCurrentFile()->Close();

	  std::stringstream ss2; // store the name of the output file and include the i and j indices!
	  std::string bkg2 = signal ? "sig": "bkg";
	  ss2 << AlgoNames[i] << "_" << i << "_" << pTbins[j] << "_" << bkg << ".nevents";
	  ofstream ev_out(ss2.str());
	  for (std::map<int,float>::iterator it = NEvents_weighted.begin(); it!= NEvents_weighted.end(); it++)
	    ev_out << it->first << "," << NEvents_weighted[it->first] << std::endl;
	  ev_out.close();
	  delete outfile;
	}
    }
    // not doing this yet.... TODO
    /*for (int j=0; j<nFineBins; j++){
      mass_max = TopEdgeMassWindow_finePt[i][j];
      mass_min = BottomEdgeMassWindow_finePt[i][j];
      }*/
    
  }


} // makeMassWindowFile()


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

/*void plotVariables(TTree * tree, vector<std:string> & branches)
  {
  
  
  
  } // plotVariables()*/


void setMassBranch(TTree * tree, std::string &algorithm, int groomAlgoIndex)
{
  if (algorithm.find("AntiKt") == std::string::npos)
    {
      tree->SetBranchAddress(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+"LC"+algorithm+"_m").c_str(), &currTree_mass);
      tree->SetBranchStatus(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+"LC"+algorithm+"_m").c_str(), 1);
    }
  else
    {
      tree->SetBranchAddress(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+algorithm+"_m").c_str(), &currTree_mass);
      tree->SetBranchStatus(std::string("jet_"+AlgoPrefix[groomAlgoIndex]+""+algorithm+"_m").c_str(), 1);
    }

} // setMassBranch()



void addJets(TTree * tree, std::string &groomalgo, bool signal, int groomIdx)
{
  //if (addJetIndex)
  //{
  tree->Branch("leadTruthIndex",&leadTruthIndex, "leadTruthIndex/I");
  tree->Branch("leadTopoIndex",&leadTopoIndex, "leadTopoIndex/I");
  tree->Branch("leadGroomedIndex",&leadGroomedIndex, "leadGroomedIndex/I");
  //tree->Branch("pt_reweight",&pt_reweight, "pt_reweight/F");
  tree->Branch("normalisation",&normalisation, "normalisation/F");
  tree->Branch("NEvents",&NEvents,"NEvents/I");
  //tree->Branch("NEvents_weighted",&NEvents_weighted,"NEvents_weighted/F");
  //tree->Branch("NEvents_weighted",&NEvents_weighted,"NEvents_weighted/F");
  tree->SetBranchAddress("mc_event_weight",&mc_event_weight);
  tree->SetBranchAddress("mc_channel_number", &mc_channel_number);
  //}
  std::string samplePrefix = "";
  bool addLC = false;
  samplePrefix = AlgoPrefix[groomIdx];
  /*if (groomalgo.find("AntiKt") != string::npos)
    {
      samplePrefix = "AntiKt10";
    }
  
  else
    {
      samplePrefix = "CamKt12"; // this is not a good solution
      addLC = true;
    }
  */

  if (groomIdx > 6) // we're doing reclustering
    addLC = true; // just add the LC to the name

  for (int i = 0; i < 3; i++) // truth, topo, groomed
    {
      
      std::string jetType = ""; //set to truth/ topo/ groomed
      switch (i)
	{
	case 0: // truth
	  jetType="jet_CamKt12Truth_";
	  //jetType = "jet_" + samplePrefix + "Truth_";
	  //jetType = "jet_" + samplePrefix + "Truth"+groomalgo.substr(4,groomalgo.length())+"_"; // added groomalgo+"_"
	  break;
	case 1: // topo
	  //jetType="jet_CamKt12LCTopo_";
	  if (addLC)
	    jetType = "jet_" + samplePrefix + "LCTopo_";
	  else
	    jetType = "jet_" + samplePrefix + "Topo_";
	  break;
	default: // groomed
	  //jetType="jet_CamKt12LC"+groomalgo+"_";
	  if (addLC)
	    jetType = "jet_" + samplePrefix +"LC" + groomalgo + "_";
	  else
	    jetType = "jet_" + samplePrefix + groomalgo + "_"; 
	  
	}
      if (signal)
	{

      tree->SetBranchAddress(std::string(jetType+"E").c_str(),&signal_E_vec[i]);
      tree->SetBranchAddress(std::string(jetType+"pt").c_str(),&signal_pt_vec.at(i));
      tree->SetBranchAddress(std::string(jetType+"m").c_str(),&signal_m_vec.at(i));
      tree->SetBranchAddress(std::string(jetType+"eta").c_str(),&signal_eta_vec.at(i));
      tree->SetBranchAddress(std::string(jetType+"phi").c_str(),&signal_phi_vec.at(i));
      if (i!= 0)
	tree->SetBranchAddress(std::string(jetType+"emfrac").c_str(),&signal_emfrac_vec.at(i));
      tree->SetBranchAddress(std::string(jetType+"Tau1").c_str(),&signal_Tau1_vec.at(i));
      tree->SetBranchAddress(std::string(jetType+"Tau2").c_str(),&signal_Tau2_vec.at(i));
      tree->SetBranchAddress(std::string(jetType+"Tau3").c_str(),&signal_Tau3_vec.at(i));
      tree->SetBranchAddress(std::string(jetType+"WIDTH").c_str(),&signal_WIDTH_vec.at(i));
      tree->SetBranchAddress(std::string(jetType+"SPLIT12").c_str(),&signal_SPLIT12_vec.at(i));
      tree->SetBranchAddress(std::string(jetType+"SPLIT23").c_str(),&signal_SPLIT23_vec.at(i));
      tree->SetBranchAddress(std::string(jetType+"SPLIT34").c_str(),&signal_SPLIT34_vec.at(i));
      tree->SetBranchAddress(std::string(jetType+"Dip12").c_str(),&signal_Dip12_vec.at(i));
      tree->SetBranchAddress(std::string(jetType+"Dip13").c_str(),&signal_Dip13_vec.at(i));
      tree->SetBranchAddress(std::string(jetType+"Dip23").c_str(),&signal_Dip23_vec.at(i));
      tree->SetBranchAddress(std::string(jetType+"DipExcl12").c_str(),&signal_DipExcl12_vec.at(i));
      tree->SetBranchAddress(std::string(jetType+"PlanarFlow").c_str(),&signal_PlanarFlow_vec.at(i));
      tree->SetBranchAddress(std::string(jetType+"Angularity").c_str(),&signal_Angularity_vec.at(i));
      tree->SetBranchAddress(std::string(jetType+"QW").c_str(),&signal_QW_vec.at(i));
      tree->SetBranchAddress(std::string(jetType+"PullMag").c_str(),&signal_PullMag_vec.at(i));
      tree->SetBranchAddress(std::string(jetType+"PullPhi").c_str(),&signal_PullPhi_vec.at(i));
      tree->SetBranchAddress(std::string(jetType+"Pull_C00").c_str(),&signal_Pull_C00_vec.at(i));
      tree->SetBranchAddress(std::string(jetType+"Pull_C01").c_str(),&signal_Pull_C01_vec.at(i));
      tree->SetBranchAddress(std::string(jetType+"Pull_C10").c_str(),&signal_Pull_C10_vec.at(i));
      tree->SetBranchAddress(std::string(jetType+"Pull_C11").c_str(),&signal_Pull_C11_vec.at(i));
	}
      else
	{
	  tree->SetBranchAddress(std::string(jetType+"E").c_str(),&bkg_E_vec.at(i));
	  tree->SetBranchAddress(std::string(jetType+"pt").c_str(),&bkg_pt_vec.at(i));
	  tree->SetBranchAddress(std::string(jetType+"m").c_str(),&bkg_m_vec.at(i));
	  tree->SetBranchAddress(std::string(jetType+"eta").c_str(),&bkg_eta_vec.at(i));
	  tree->SetBranchAddress(std::string(jetType+"phi").c_str(),&bkg_phi_vec.at(i));
	  if (i != 0)
	    tree->SetBranchAddress(std::string(jetType+"emfrac").c_str(),&bkg_emfrac_vec.at(i));
	  tree->SetBranchAddress(std::string(jetType+"Tau1").c_str(),&bkg_Tau1_vec.at(i));
	  tree->SetBranchAddress(std::string(jetType+"Tau2").c_str(),&bkg_Tau2_vec.at(i));
	  tree->SetBranchAddress(std::string(jetType+"Tau3").c_str(),&bkg_Tau3_vec.at(i));
	  tree->SetBranchAddress(std::string(jetType+"WIDTH").c_str(),&bkg_WIDTH_vec.at(i));
	  tree->SetBranchAddress(std::string(jetType+"SPLIT12").c_str(),&bkg_SPLIT12_vec.at(i));
	  tree->SetBranchAddress(std::string(jetType+"SPLIT23").c_str(),&bkg_SPLIT23_vec.at(i));
	  tree->SetBranchAddress(std::string(jetType+"SPLIT34").c_str(),&bkg_SPLIT34_vec.at(i));
	  tree->SetBranchAddress(std::string(jetType+"Dip12").c_str(),&bkg_Dip12_vec.at(i));
	  tree->SetBranchAddress(std::string(jetType+"Dip13").c_str(),&bkg_Dip13_vec.at(i));
	  tree->SetBranchAddress(std::string(jetType+"Dip23").c_str(),&bkg_Dip23_vec.at(i));
	  tree->SetBranchAddress(std::string(jetType+"DipExcl12").c_str(),&bkg_DipExcl12_vec.at(i));
	  tree->SetBranchAddress(std::string(jetType+"PlanarFlow").c_str(),&bkg_PlanarFlow_vec.at(i));
	  tree->SetBranchAddress(std::string(jetType+"Angularity").c_str(),&bkg_Angularity_vec.at(i));
	  tree->SetBranchAddress(std::string(jetType+"QW").c_str(),&bkg_QW_vec.at(i));
	  tree->SetBranchAddress(std::string(jetType+"PullMag").c_str(),&bkg_PullMag_vec.at(i));
	  tree->SetBranchAddress(std::string(jetType+"PullPhi").c_str(),&bkg_PullPhi_vec.at(i));
	  tree->SetBranchAddress(std::string(jetType+"Pull_C00").c_str(),&bkg_Pull_C00_vec.at(i));
	  tree->SetBranchAddress(std::string(jetType+"Pull_C01").c_str(),&bkg_Pull_C01_vec.at(i));
	  tree->SetBranchAddress(std::string(jetType+"Pull_C10").c_str(),&bkg_Pull_C10_vec.at(i));
	  tree->SetBranchAddress(std::string(jetType+"Pull_C11").c_str(),&bkg_Pull_C11_vec.at(i));
	}
    } // end for loop over topo/truth/groom
} //addJets()


//void getBranchesSelection(TTree * tree, std::string & algorithm)
void setSelectionVectors(bool signal, std::string & algorithm)
{  
  
  if (signal)
    {
      jet_eta_truth = signal_eta_vec[0];
      jet_phi_truth = signal_phi_vec[0];
      jet_pt_truth = signal_pt_vec[0];

      jet_eta_topo = signal_eta_vec[1];
      jet_phi_topo = signal_phi_vec[1];
      jet_pt_topo = signal_pt_vec[1];
      jet_emfrac_topo = signal_emfrac_vec[1];
      
      jet_eta_groomed = signal_eta_vec[2];
      jet_phi_groomed = signal_phi_vec[2];
      jet_pt_groomed = signal_pt_vec[2];
      jet_emfrac_groomed = signal_emfrac_vec[2];
    }
  else
    {
      jet_eta_truth = bkg_eta_vec[0];
      jet_phi_truth = bkg_phi_vec[0];
      jet_pt_truth = bkg_pt_vec[0];

      jet_eta_topo = bkg_eta_vec[1];
      jet_phi_topo = bkg_phi_vec[1];
      jet_pt_topo = bkg_pt_vec[1];
      jet_emfrac_topo = bkg_emfrac_vec[1];
      
      jet_eta_groomed = bkg_eta_vec[2];
      jet_phi_groomed = bkg_phi_vec[2];
      jet_pt_groomed = bkg_pt_vec[2];
      jet_emfrac_groomed = bkg_emfrac_vec[2];
    }
  
}// setSelectionVectors()

  void initVectors()
  {

    // have the vectors for the above histograms so we can do the reading in stuff from the TTree
    jet_eta_truth = 0;
    jet_phi_truth = 0;
    jet_pt_truth = 0;
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
	bkg_emfrac_vec[i] = 0;
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
  }
