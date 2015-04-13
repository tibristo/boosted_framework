#include "TH1F.h"
#include "TH2F.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TF1.h>
#include <TAxis.h>
#include <iostream>
#include <iomanip>
#include "TGraphErrors.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TVector2.h"
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "AtlasUtils.C"
#include "AtlasStyle.C"
#include "TProfile.h"
#include "TGraphAsymmErrors.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "TExec.h"
#include "TPad.h"
#include "TMatrix.h"
#include "TRandom3.h"

void figure(TH1F* histos2,TH1F* histos1, TString outname,vector<TString> legends,TString aa="", TString bb=""){

	TCanvas *mycaa = new TCanvas("","",500,500);
	histos1->GetXaxis()->SetTitle(legends[2]);
	histos1->GetYaxis()->SetTitleOffset(1.5);		
	histos1->GetYaxis()->SetTitle("Normalized to Unity");

	histos1->GetXaxis()->SetNdivisions(505);
	histos1->SetLineWidth(2);
	histos1->SetLineColor(kViolet+9);

	histos2->SetLineWidth(2);
	histos2->SetLineColor(kMagenta+3);

	histos1->GetYaxis()->SetRangeUser(0.0000000001,histos1->GetMaximum()*2.5);
	histos1->GetXaxis()->SetNdivisions(605);

	histos1->Draw();	
	histos2->Draw("same");			

	TLegend* leggaa = new TLegend(.7,.5,0.9,.7);
	leggaa->SetTextFont(42);
	leggaa->AddEntry(histos1,legends[0],"lep");	
	leggaa->AddEntry(histos2,legends[1],"lep");	
			
	myText(0.2,0.88,kBlack,"#scale[1.1]{#bf{#it{ATLAS}} Simulation Internal}");
	myText(0.2,0.82,kBlack,"#scale[0.8]{#sqrt{s} = 14 TeV}");
	myText(0.2,0.75,kBlack,"#scale[0.8]{500 GeV < p_{T}^{C/A R = 1.2 truth jet} < 1000 GeV}");	
	myText(0.2,0.68,kBlack,"#scale[0.8]{|#eta^{C/A R = 1.2 truth jet}| < 1.2}");	

	myText(0.5,0.4,kBlack,aa);
	myText(0.5,0.35,kBlack,bb);

	leggaa->SetFillStyle(0);
	leggaa->SetFillColor(0);
	leggaa->SetBorderSize(0);
	leggaa->Draw();

	mycaa->Print("../plots/"+outname+".pdf");

}

void window(map<int,TH1F*> histos2,map<int,TH1F*> histos1, map<int,TString> algorithm_labels){

	TCanvas *mycaa = new TCanvas("","",500,500);

	double Ws[4];
	double QCDs[4];
	double windowsize[4];

	for (int k=0; k<4; k++){
		int n = histos1[k]->GetNbinsX();
		int left=-1;
		int right=-1;
		double close=99999999;;
		for (int i=1; i<=n; i++){
			for (int j=i+1; j<=n; j++){
				double frac = histos1[k]->Integral(i,j)/histos1[k]->Integral();
				if (j-i<close && frac > 0.68){
					close = j-i;
					left=i;
					right=j;
				}
			}
		}

		double W1 = histos1[k]->Integral(left,right)/histos1[k]->Integral();
		double QCD1 = histos2[k]->Integral(left,right)/histos2[k]->Integral();

		Ws[k]=W1;
		QCDs[k]=QCD1;
		windowsize[k] = 3*(right-left);
	}

	double Rs[4]={0,1,2,3};

	TGraph* effics = new TGraph(4,Rs,Ws);
	TGraph* effics2 = new TGraph(4,Rs,QCDs);
	TGraph* wsize = new TGraph(4,Rs,windowsize);

	TH1F* dl = new TH1F("","",7,-0.5,6.5);
	for (int i=1; i<=4; i++){
		dl->GetXaxis()->SetBinLabel(i,algorithm_labels[i-1]);
	}	
	dl->GetXaxis()->SetTitleOffset(1.7);	
	dl->GetXaxis()->SetTitle("Jet Algorithm");
	dl->GetYaxis()->SetTitle("Efficiency");	
	dl->GetYaxis()->SetTitleOffset(1.4);		
	dl->GetYaxis()->SetRangeUser(0,1);
	dl->Draw();
	effics->SetMarkerStyle(20);
	effics->SetMarkerSize(1);
	effics2->SetMarkerStyle(20);
	effics2->SetMarkerSize(1);
	effics->Draw("samep");			
	effics->SetMarkerColor(2);
	effics2->Draw("samep");

	TLegend* leggaa = new TLegend(.7,.5,0.9,.7);
	leggaa->SetTextFont(42);
	leggaa->AddEntry(effics,"W","p");		
	leggaa->AddEntry(effics2,"QCD","p");	
			
	myText(0.2,0.88,kBlack,"#scale[1.1]{#bf{#it{ATLAS}} Simulation Internal}");
	myText(0.2,0.82,kBlack,"#scale[0.8]{#sqrt{s} = 14 TeV, m_{W'} = 3 TeV}");

	leggaa->SetFillStyle(0);
	leggaa->SetFillColor(0);
	leggaa->SetBorderSize(0);
	leggaa->Draw();

	mycaa->Print("../plots/effics.pdf");

	/////-----------

	TH1F* dl2 = new TH1F("","",6,-0.5,5.5);
	for (int i=1; i<=4; i++){
		dl2->GetXaxis()->SetBinLabel(i,algorithm_labels[i-1]);
	}	
	dl2->GetXaxis()->SetTitleOffset(1.7);	
	dl2->GetXaxis()->SetTitle("Jet Algorithm");
	dl2->GetYaxis()->SetTitle("Window Size [GeV]");	
	dl2->GetYaxis()->SetTitleOffset(1.4);		
	dl2->GetYaxis()->SetRangeUser(0,100);
	dl2->Draw();
	wsize->SetMarkerStyle(20);
	wsize->SetMarkerSize(1);
	wsize->Draw("samep");			
			
	myText(0.2,0.88,kBlack,"#scale[1.1]{#bf{#it{ATLAS}} Simulation Internal}");
	myText(0.2,0.82,kBlack,"#scale[0.8]{#sqrt{s} = 14 TeV, m_{W'} = 3 TeV}");

	mycaa->Print("../plots/window.pdf");	

}
