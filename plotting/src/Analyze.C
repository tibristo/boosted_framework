#include "Helper.h"
#include "TRandom3.h"
#include "TH3F.h"

void analyze(map<int,TFile*>);

void makeplots(){

	//First setup the files
	TString W_file = TString("../data/BB3_lily.WReclusteredCOMMON.jetmet2012.root");
	TFile *W = new TFile(W_file);

	TString QCD_file = TString("../data/BB3_lily.PythJXmc12_14TeV_COMMON.jetmet2012.root");
	TFile *QCD = new TFile(QCD_file);

	map<int,TFile*>files;
	files[20]=QCD; //don't ask why I call QCD 20 and W's 40...
	files[40]=W;

	analyze(files);

}	

void analyze(map<int,TFile*> myfiles){

	gStyle->SetOptStat(0);
	gROOT->Reset();
	gROOT->SetStyle("ATLAS");
	gROOT->ForceStyle();
	gStyle->SetPadLeftMargin(0.16);

	map<int,TTree*> mytrees;	

	Int_t           RunNumber;
	Float_t         Weight;

	//I preprocess the D3PD's so that all I save is the M=mass and P=pT of the relevant algorithms.

	Float_t         rJP4; //RT, fcut=0.05, r=0.4, pTcut=15;
	Float_t         rJM4; 
	Float_t         rJP3; //RT, fcut=0.05, r=0.3, pTcut=15;
	Float_t         rJM3;
	Float_t         rJP2; //RT, fcut=0.05, r=0.2, pTcut=15;
	Float_t         rJM2; 
	Float_t         fJP; //CamKt8LCTopoPrunedCaRcutFactor50Zcut10
	Float_t         fJM; 
	Float_t         fJP2; //AntiKt10LCTopoTrimmedPtFrac5SmallR20
	Float_t         fJM2;
	Float_t         fJP3; //AntiKt10LCTopoTrimmedPtFrac5SmallR30
	Float_t         fJM3; 
	Float_t         fJP4; //AntiKt10LCTopoTrimmedPtFrac5SmallR40
	Float_t         fJM4; 

	//This is the C/A R=1.2 Truth Jet which we match to.  For this, I also save eta, since we bin in that.
	Float_t         fJP12;
	Float_t         fJM12;
	Float_t         fJEta12;

	//Setup the trees with all the branches.
	int processes[2]={20,40};
	for (int i=0; i<2; i++){
		mytrees[processes[i]]=((TTree*)(myfiles[processes[i]]->Get("lilytree")));	

		mytrees[processes[i]]->SetBranchAddress("rJM2", &rJM2);
		mytrees[processes[i]]->SetBranchAddress("rJP2", &rJP2);
		mytrees[processes[i]]->SetBranchAddress("rJM3", &rJM3);
		mytrees[processes[i]]->SetBranchAddress("rJP3", &rJP3);
		mytrees[processes[i]]->SetBranchAddress("rJM4", &rJM4);
		mytrees[processes[i]]->SetBranchAddress("rJP4", &rJP4);
		mytrees[processes[i]]->SetBranchAddress("fJM", &fJM);
		mytrees[processes[i]]->SetBranchAddress("fJP", &fJP);
		mytrees[processes[i]]->SetBranchAddress("fJM2", &fJM2);
		mytrees[processes[i]]->SetBranchAddress("fJP2", &fJP2);
		mytrees[processes[i]]->SetBranchAddress("fJM3", &fJM3);
		mytrees[processes[i]]->SetBranchAddress("fJP3", &fJP3);
		mytrees[processes[i]]->SetBranchAddress("fJM4", &fJM4);
		mytrees[processes[i]]->SetBranchAddress("fJP4", &fJP4);

		mytrees[processes[i]]->SetBranchAddress("fJP12", &fJP12);
		mytrees[processes[i]]->SetBranchAddress("fJM12", &fJM12);
		mytrees[processes[i]]->SetBranchAddress("fJEta12", &fJEta12);		

		mytrees[processes[i]]->SetBranchAddress("RunNumber", &RunNumber);
		mytrees[processes[i]]->SetBranchAddress("Weight", &Weight);		
	}

	map<TString,TH1F*> histos20;
	map<TString,TH1F*> histos40;			
	map<TString,TString> histolabels;
	vector<TString>histonames;
	histonames.push_back("rJM2"); histolabels["rJM2"] = "rJM2";
	histonames.push_back("rJM3"); histolabels["rJM3"] = "rJM3";
	histonames.push_back("rJM4"); histolabels["rJM4"] = "rJM4";
	histonames.push_back("fJM"); histolabels["fJM"] = "fJM";

	histonames.push_back("rJP2"); histolabels["rJP2"] = "rJP2";
	histonames.push_back("rJP3"); histolabels["rJP3"] = "rJP3";
	histonames.push_back("rJP4"); histolabels["rJP4"] = "rJP4";
	histonames.push_back("fJP"); histolabels["fJP"] = "fJP";

	histonames.push_back("fJP12"); histolabels["fJP12"] = "fJP12";
	histonames.push_back("fJM12"); histolabels["fJM12"] = "fJM12";
	histonames.push_back("fJEta12"); histolabels["fJEta12"] = "fJEta12";

	//Set up some histogram settings
	double pTend = 2000; //range of the pT histos
	double massend = 300; //range of the mass histos
	int binnumb=100; //number of pT/mass bins

	for (unsigned int i=0; i< histonames.size(); i++){

		double range_end = pTend;
		if (histonames[i].Contains("M")) range_end = massend;
		TH1F* hold = new TH1F(histonames[i],histonames[i],binnumb,0,range_end);
		hold->Sumw2();
		histos20[histonames[i]]=hold;
	}			
	for (unsigned int i=0; i< histonames.size(); i++){

		double range_end = pTend;
		if (histonames[i].Contains("M")) range_end = massend;		
		TH1F* hold = new TH1F(histonames[i],histonames[i],binnumb,0,range_end);
		hold->Sumw2();		
		histos40[histonames[i]]=hold;
	}

	//Setup some selections for the range.
	double pTcut = 500.*1000;
	double pTcuthigh = 1000*1000;
	double etalow=0.;
	double etahigh=1.2;

	//Now, let's fill the QCD histos.
	for (int i=0; i<mytrees[20]->GetEntries(); i++){
        mytrees[20]->GetEntry(i);
        //These are the cross section and filter efficiency weights for the QCD samples.
        if (RunNumber == 147913) Weight*=1.6664E+03*1.9139E-03; 
        if (RunNumber == 147914) Weight*=2.7646E+01*1.4296E-03; 
        if (RunNumber == 147915) Weight*=3.0317E-01*5.5040E-03; 
        if (RunNumber == 147916) Weight*=7.5078E-03*1.5252E-02;          
        if (RunNumber == 147917) Weight*=1.3760E-03*7.6369E-02;                          
        if (fJP12 > pTcut && fJP12 < pTcuthigh && fabs(fJEta12) > etalow && fabs(fJEta12) < etahigh){
			histos20["rJM4"]->Fill(rJM4/1000,Weight);
        	histos20["rJM3"]->Fill(rJM3/1000,Weight);        
			histos20["rJM2"]->Fill(rJM2/1000,Weight);          						     
			histos20["fJM"]->Fill(fJM/1000,Weight); 
        	histos20["rJP4"]->Fill(rJP4/1000,Weight);
        	histos20["rJP3"]->Fill(rJP3/1000,Weight);        
			histos20["rJP2"]->Fill(rJP2/1000,Weight); 
			histos20["fJP"]->Fill(fJP/1000,Weight); 
			histos20["fJP12"]->Fill(fJP12/1000,Weight); 
		}
    }

    //Now, we fill the signal histos.
	for (int i=0; i<mytrees[40]->GetEntries(); i++){
        mytrees[40]->GetEntry(i);
	    if (fJP12 > pTcut && fJP12 < pTcuthigh && fabs(fJEta12) > etalow && fabs(fJEta12) < etahigh){
	        histos40["rJM4"]->Fill(rJM4/1000,Weight);
	        histos40["rJM3"]->Fill(rJM3/1000,Weight);        
			histos40["rJM2"]->Fill(rJM2/1000,Weight);          	 				  
			histos40["fJM"]->Fill(fJM/1000,Weight); 
	        histos40["rJP4"]->Fill(rJP4/1000,Weight);
	        histos40["rJP3"]->Fill(rJP3/1000,Weight);        
			histos40["rJP2"]->Fill(rJP2/1000,Weight); 
			histos40["fJP"]->Fill(fJP/1000,Weight); 
			histos40["fJP12"]->Fill(fJP12/1000,Weight);
		}
    }       

    //Now, we are going to re-weight to the truth C/A R=1.2 truth pT spectrum.
	std::cout << " Doing the re-weighting " << std::endl;	
	TH1F* weights = (TH1F*)histos20["fJP"]->Clone("");	
	for (int i=1; i<=weights->GetNbinsX(); i++){
		if (histos40["fJP12"]->GetBinContent(i) > 0){
			weights->SetBinContent(i,histos20["fJP12"]->GetBinContent(i)/histos40["fJP12"]->GetBinContent(i));
		}
		else{
			weights->SetBinContent(i,1);
		}
	}
	std::cout << " ====================== " << std::endl;	

	for (unsigned int i=0; i< histonames.size(); i++){
		for (int j=1; j<=histos40[histonames[i]]->GetNbinsX(); j++){
		histos40[histonames[i]]->SetBinContent(j,0.);				
		}
	}

	for (int i=0; i<mytrees[40]->GetEntries(); i++){
        mytrees[40]->GetEntry(i);
        if (fJP12 > pTcut && fJP12 < pTcuthigh && fabs(fJEta12) > etalow && fabs(fJEta12) < etahigh){
        	Weight*=weights->GetBinContent(weights->FindBin(fJP12/1000.));
	        histos40["rJM4"]->Fill(rJM4/1000,Weight);
	        histos40["rJM3"]->Fill(rJM3/1000,Weight);        
			histos40["rJM2"]->Fill(rJM2/1000,Weight);          	 				  
			histos40["fJM"]->Fill(fJM/1000,Weight); 
			histos40["rJP4"]->Fill(rJP4/1000,Weight);
	        histos40["rJP3"]->Fill(rJP3/1000,Weight);        
			histos40["rJP2"]->Fill(rJP2/1000,Weight); 
			histos40["fJP"]->Fill(fJP/1000,Weight); 
			histos40["fJP12"]->Fill(fJP12/1000,Weight);
		}
    } 

    //Normalize the histograms.
	for (unsigned int i=0; i< histonames.size(); i++){
		histos40[histonames[i]]->Scale(1./histos40[histonames[i]]->Integral());		
		histos20[histonames[i]]->Scale(1./histos20[histonames[i]]->Integral());						
	}

	vector<TString> mynames;
	mynames.push_back("rJM4");	
	mynames.push_back("rJM3");		
	mynames.push_back("rJM2");	
	mynames.push_back("fJM");	

	vector<TString> mynamespt;
	mynamespt.push_back("rJP4");	
	mynamespt.push_back("rJP3");		
	mynamespt.push_back("rJP2");	
	mynamespt.push_back("fJP12");

	vector<TString> mynames2;
	mynames2.push_back("W' #rightarrow WZ #rightarrow llqq'");
	mynames2.push_back("Dijets");			
	mynames2.push_back("Jet Mass [GeV]");

	//Plot some distributions
	figure(histos20[mynames[0]],histos40[mynames[0]],"RTMassr04",mynames2,"RT r = 0.4 anti-k_{t}","f_{cut} = 0.05");	
	figure(histos20[mynames[1]],histos40[mynames[1]],"RTMassr03",mynames2,"RT r = 0.3 anti-k_{t}","f_{cut} = 0.05");		
	figure(histos20[mynames[2]],histos40[mynames[2]],"RTMassr02",mynames2,"RT r = 0.2 anti-k_{t}","f_{cut} = 0.05");
	figure(histos20[mynames[3]],histos40[mynames[3]],"CATruthMass",mynames2,"C/A R=1.2","Truth Jet");	

	mynames2[2]="Jet p_{T} [GeV]";
	figure(histos20[mynamespt[0]],histos40[mynamespt[0]],"RTPTr04",mynames2,"RT r = 0.4 anti-k_{t}","f_{cut} = 0.05");	
	figure(histos20[mynamespt[1]],histos40[mynamespt[1]],"RTPTr03",mynames2,"RT r = 0.3 anti-k_{t}","f_{cut} = 0.05");		
	figure(histos20[mynamespt[2]],histos40[mynamespt[2]],"RTPTr02",mynames2,"RT r = 0.2 anti-k_{t}","f_{cut} = 0.05");
	figure(histos20[mynamespt[3]],histos40[mynamespt[3]],"CATruthPT",mynames2,"C/A R=1.2","Truth Jet");	

	map<int,TH1F*> Ws;
	map<int,TH1F*> QCDs;
	map<int,TString> algorithm_labels;
	algorithm_labels[0] = "RT, r=0.4";
	algorithm_labels[1] = "RT, r=0.3";
	algorithm_labels[2] = "RT, r=0.2";
	algorithm_labels[3] = "C/A R=1.2 Truth";			
	Ws[0]=histos20[mynames[0]];
	Ws[1]=histos20[mynames[1]];
	Ws[2]=histos20[mynames[2]];
	Ws[3]=histos20[mynames[3]];	
	QCDs[0]=histos40[mynames[0]];
	QCDs[1]=histos40[mynames[1]];
	QCDs[2]=histos40[mynames[2]];	
	QCDs[3]=histos40[mynames[3]];	
	window(Ws,QCDs,algorithm_labels);

}
