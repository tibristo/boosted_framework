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

using namespace std;

bool ComparePt(TLorentzVector a, TLorentzVector b) { return a.Pt() > b.Pt(); }


int main( int argc, char * argv[] ) {
  
  gROOT->ProcessLine("#include <vector>");
  gInterpreter->GenerateDictionary("vector<vector<int> >", "vector");
  gInterpreter->GenerateDictionary("vector<vector<double> >", "vector");
  
  AtlasStyle();
  
  //Make an array of TFiles with all the relevant inputs and add them to my array of trees
  for (int nArg=1; nArg < argc; nArg++) {
    //cout << nArg << " " << argv[nArg] << endl;
    inputFile[nArg-1] = new TFile(argv[nArg], "READ");
    inputTree[nArg-1] = ( TTree* ) inputFile[nArg-1]->Get( "physics" );    
  }

  defineStrings(AlgoList, binLabel, pTbins, finePtBins);
  createHistos();

 
  initializeVariables();
  getMassHistograms(inputTree[0], inputTree[1], "TopoSplitFilteredMu67SmallR0YCut9", 0);
  initializeVariables();
  getMassHistograms(inputTree[0], inputTree[1], "TopoSplitFilteredMu67SmallR0YCut9", 1);
  initializeVariables();
  getMassHistograms(inputTree[0], inputTree[1], "TopoSplitFilteredMu100SmallR30YCut4", 2);
  initializeVariables();
  getMassHistograms(inputTree[2], inputTree[3], "TopoTrimmedPtFrac5SmallR30", 3);
  initializeVariables();
  getMassHistograms(inputTree[2], inputTree[3], "TopoTrimmedPtFrac5SmallR20", 4);
  initializeVariables();
  getMassHistograms(inputTree[4], inputTree[5], "TopoPrunedCaRcutFactor50Zcut10", 5);
  initializeVariables();
  getMassHistograms(inputTree[4], inputTree[5], "TopoPrunedCaRcutFactor50Zcut20", 6);
  //reclustering
  initializeVariables();
  getMassHistograms(inputTree[6], inputTree[7], "AntiKt2LCTopo", 7);
  deleteVectors();
  initializeVariables();
  getMassHistograms(inputTree[6], inputTree[7], "AntiKt3LCTopo", 8);
  deleteVectors();
  initializeVariables();
  getMassHistograms(inputTree[6], inputTree[7], "AntiKt4LCTopo", 9);
  deleteVectors();
  

  //Plot things
  //Normalise all histograms of interest
  // this is not the right way to do the normalisation according to this - https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/BoostedBosonTaggingD3PD#MC_only_studies

  // NEvents = Sum[mcevt_weight[0][0]*PileupWeight
  // Note, for NEvents for SherpaW + jets we need to call getNormSherpaW() because the weighting is slightly different
  // Normalization = CrossSection(at NLO) * Luminosity / NEvents


  for (int i=0; i<nAlgos-1; i++){
    
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

    qcd_PtReweight[i]->Scale(1.0/qcd_PtReweight[i]->Integral());
    Wp_PtReweight[i]->Scale(1.0/Wp_PtReweight[i]->Integral());
    Wprime_Lead_CA12_scaled_pt[i]->Scale(1.0/Wprime_Lead_CA12_scaled_pt[i]->Integral());


  }

  ////////PLOT THINGS -> This should go into it's own class or at the very least it's own method, TODO

  // CANVAS FOR THE MASS PLOTS
  
  for (int i=0; i<nAlgos-2; i++){

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
  
  for (int i=0; i<nAlgos-2; i++){
    
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
  


  //NOW THAT I HAVE ALL MY MASS PLOTS NORMALISED, I WANT TO GET THE MPV AND MASS WINDOW ANALYSIS RIGHT

  // QCD VECTOR: qcd_Lead_CA12_mass[i][j]
  // WPRIME VECTOR: Wprime_Lead_CA12_mass[i][j]

  //1. Get the MPV

  for (int i=0; i<nAlgos-1; i++){
    for (int j=0; j<nPtBins; j++){
      myMPV[i][j]=mpv(Wprime_Lead_CA12_mass[i][j]);
    }

    for (int j=0; j<nFineBins; j++){
      myMPV_finePt[i][j]=mpv(Wprime_finePtBin_mass[i][j]);
    }
  }
  
  //2. Get the mass window which gives 68% W mass efficiency 

  for (int i=0; i<nAlgos-1; i++){
    for (int j=0; j<nPtBins; j++){
      Qw(WidthMassWindow[i][j],TopEdgeMassWindow[i][j],Wprime_Lead_CA12_mass[i][j], 0.68);
      //one and two are outputs, three and four are inputs
      BottomEdgeMassWindow[i][j]=TopEdgeMassWindow[i][j]-WidthMassWindow[i][j];    
      cout << AlgoList[i] << ": top edge " << TopEdgeMassWindow[i][j] << " bottom edge " << BottomEdgeMassWindow[i][j] << " mass window " << WidthMassWindow[i][j] << endl;
    }

    for (int j=0; j<nFineBins; j++){
      Qw(WidthMassWindow_finePt[i][j],TopEdgeMassWindow_finePt[i][j],Wprime_finePtBin_mass[i][j], 0.68);
      //one and two are outputs, three and four are inputs
      BottomEdgeMassWindow_finePt[i][j]=TopEdgeMassWindow_finePt[i][j]-WidthMassWindow_finePt[i][j];    
      //cout << AlgoList[i] << ": top edge " << TopEdgeMassWindow[i][j] << " bottom edge " << BottomEdgeMassWindow[i][j] << " mass window " << WidthMassWindow[i][j] << endl;
    }
   
  }

  //3. Check background fraction in this window

  for (int i=0; i<nAlgos-1; i++){
    for (int j=0; j<nPtBins; j++){
      QCDfrac[i][j]=qcd_Lead_CA12_mass[i][j]->Integral(qcd_Lead_CA12_mass[i][j]->FindBin(BottomEdgeMassWindow[i][j]),qcd_Lead_CA12_mass[i][j]->FindBin(TopEdgeMassWindow[i][j]));
      cout << "QCD fraction " << QCDfrac[i][j] << endl;
    }

    for (int j=0; j<nFineBins; j++){
      QCDfrac_finePt[i][j]=qcd_finePtBin_mass[i][j]->Integral(qcd_finePtBin_mass[i][j]->FindBin(BottomEdgeMassWindow_finePt[i][j]),qcd_finePtBin_mass[i][j]->FindBin(TopEdgeMassWindow_finePt[i][j]));
      //cout << "QCD fraction " << QCDfrac[i][j] << endl;

    }
    
  }
  

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
  for (int i=0; i<nAlgos-2; i++){  
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
  for (int i=0; i<nAlgos-2; i++){
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
  

  TCanvas *wVsPt[nAlgos-2];
  TPad *pad11[nAlgos-2];
  TPad *pad22[nAlgos-2];
	     
  for (int i=1; i<nAlgos-2; i++){
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
  
  /*
  //christos #2 sort of plots - QCD eff superposition
  
  TCanvas * c4 = new TCanvas("canvas", "canvas", 700,550);
  hQCDeff[1]->SetLineColor(1);
  hQCDeff[4]->SetLineColor(2);
  hQCDeff[6]->SetLineColor(3);
  hQCDeff[7]->SetLineColor(4);
  
  hQCDeff[1]->Draw();
  hQCDeff[4]->Draw("same");
  hQCDeff[6]->Draw("same");
  hQCDeff[7]->Draw("same");


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
  
  TLegend *leg = new TLegend(0.65,0.88,0.93,0.97);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.040);
  leg->SetNColumns(2);
  leg->AddEntry( hQCDeff[1] , binLabel[1] , "l" );
  leg->AddEntry( hQCDeff[4] , binLabel[4] , "l" );
  leg->AddEntry( hQCDeff[6] , binLabel[6] , "l" );
  leg->AddEntry( hQCDeff[7] , binLabel[7] , "l" );
  leg->Draw();
  
  
  c4->SaveAs("QCDeff_comparisons.eps");

  //end of christos plots #2
  
  */

  ///// SAVE THINGS TO HISTOGRAM FILE
  
  TFile f("MunichPlots.root","recreate");
  
  for (int i=0; i<nAlgos-2; i++){
    windowsVsPt[i]->Write();
  }

  for (int i=0; i<nAlgos-1; i++){
    qcd_Lead_CA12_pt[i]->Write();
    Wprime_Lead_CA12_pt[i]->Write();
    
    pTweights[i]->Write();
    qcd_PtReweight[i]->Write();
    Wp_PtReweight[i]->Write();
    Wprime_Lead_CA12_scaled_pt[i]->Write();


    for (int j=0; j<nPtBins; j++){
      qcd_Lead_CA12_mass[i][j]->Write();
      Wprime_Lead_CA12_mass[i][j]->Write();
    }

    for (int j=0; j<nFineBins; j++){
      qcd_finePtBin_mass[i][j]->Write();
      Wprime_finePtBin_mass[i][j]->Write();
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
   cout<<"MPV: "<<mostprob<<endl;                    
                                  
   return mostprob;
}


void Qw(double &minWidth, double &topEdge, TH1F* histo, double frac){
  
  minWidth = 100000.0;
  topEdge = 0.0;
  
  int Nbins = histo->GetNbinsX();
  
  for(int i=0; i!=Nbins; i++){
    double tempFrac = 0;
    int imax = i;
    
    while(tempFrac<frac && imax != Nbins){
      //fraction in bin imax=0,1,2,...                                                                                       
      tempFrac+=histo->GetBinContent(imax)/histo->Integral();
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



void getMassHistograms(TTree *inputTree, TTree *inputTree1, TString groomAlgo, int groomAlgoIndex){
  
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
  
  cout << "QCD efficiencies for " << groomAlgo << ": " << nEvt_0 << " " << nEvt_1 << " " << nEvt_2  << " " << nEvt_3  << " " << nEvt_4 << endl;
  
  
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
    cout << weight << endl;
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
	if ((*Wp_CA12_groomed_pt)[chosenLeadGroomedIndex]/1000<500.00)
	  cout << (*Wp_CA12_groomed_pt)[chosenLeadGroomedIndex]/1000 << endl;
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
  
  for (int i=0; i<nAlgos-1; i++){
    
    for (int j=0; j<nFineBins; j++){
      qcd_finePtBin_mass[i][j] = new TH1F("qcd_finePtBin_mass_"+AlgoList[i]+finePtBins[j], "qcd_finePtBin_mass_"+AlgoList[i]+finePtBins[j], 150, 0.0, 300.0);
      Wprime_finePtBin_mass[i][j] = new TH1F("Wprime_finePtBin_mass_"+AlgoList[i]+finePtBins[j], "Wprime_finePtBin_mass_"+AlgoList[i]+finePtBins[j], 150, 0.0, 300.0);
      qcd_finePtBin_mass[i][j]->Sumw2();
      Wprime_finePtBin_mass[i][j]->Sumw2();
    }


    for (int j=0; j<nPtBins; j++){
      qcd_Lead_CA12_mass[i][j] = new TH1F("qcd_Lead_CA12_mass_"+AlgoList[i]+pTbins[j], "qcd_Lead_CA12_mass_"+AlgoList[i]+pTbins[j], 300, 0.0, 300.0);
      qcd_Lead_CA12_mass[i][j]->Sumw2();
      Wprime_Lead_CA12_mass[i][j] = new TH1F("Wprime_Lead_CA12_mass_"+AlgoList[i]+pTbins[j], "Wprime_Lead_CA12_mass_"+AlgoList[i]+pTbins[j], 300, 0.0, 300.0);
      Wprime_Lead_CA12_mass[i][j]->Sumw2();
    }

    qcd_Lead_CA12_pt[i] = new TH1F("qcd_Lead_CA12_pt_"+AlgoList[i], "qcd_Lead_CA12_pt_"+AlgoList[i], 20, 0.0, 3500.0);
    Wprime_Lead_CA12_pt[i] = new TH1F("Wprime_Lead_CA12_pt_"+AlgoList[i], "Wprime_Lead_CA12_pt_"+AlgoList[i], 20, 0.0, 3500.0);
    Wprime_Lead_CA12_scaled_pt[i] = new TH1F("Wprime_Lead_CA12_scaled_pt_"+AlgoList[i], "Wprime_Lead_CA12_scaled_pt_"+AlgoList[i], 20, 0.0, 3500.0);
    qcd_PtReweight[i] = new TH1F("qcd_ptreweight_"+AlgoList[i], "qcd_ptreweight_"+AlgoList[i], 20, 0.0, 3500.0);
    Wp_PtReweight[i] = new TH1F("Wp_ptreweight_"+AlgoList[i], "Wp_ptreweight_"+AlgoList[i], 20, 0.0, 3500.0);
    pTweights[i] = new TH1F("pT_weights_"+AlgoList[i], "pt weights_"+AlgoList[i], 20, 0.0, 3500.0);
      
  }
  
  qcd_Lead_CA12_mass[nAlgos-1][0] = new TH1F("qcd_Lead_CA12_mass_"+AlgoList[nAlgos-1], "qcd_Lead_CA12_mass_"+AlgoList[nAlgos-1], 100, 0.0, 300.0);
  qcd_Lead_CA12_mass[nAlgos-1][0]->Sumw2();
  Wprime_Lead_CA12_mass[nAlgos-1][0] = new TH1F("Wprime_Lead_CA12_mass_"+AlgoList[nAlgos-1], "Wprime_Lead_CA12_mass_"+AlgoList[nAlgos-1], 100, 0.0, 300.0);
  Wprime_Lead_CA12_mass[nAlgos-1][0]->Sumw2();
  
  qcd_Lead_CA12_pt[nAlgos-1] = new TH1F("qcd_Lead_CA12_pt_"+AlgoList[nAlgos-1], "qcd_Lead_CA12_pt_"+AlgoList[nAlgos-1], 350, 0.0, 3500.0);
  Wprime_Lead_CA12_pt[nAlgos-1] = new TH1F("Wprime_Lead_CA12_pt_"+AlgoList[nAlgos-1], "Wprime_Lead_CA12_pt_"+AlgoList[nAlgos-1], 350, 0.0, 3500.0);
  Wprime_Lead_CA12_scaled_pt[nAlgos-1] = new TH1F("Wprime_Lead_CA12_scaled_pt_"+AlgoList[nAlgos-1], "Wprime_Lead_CA12_scaled_pt_"+AlgoList[nAlgos-1], 350, 0.0, 3500.0);
  

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
    inputTree->SetBranchAddress("mc_channel_number", &qcd_mc_channel_number);
    inputTree->SetBranchAddress("mc_event_weight", &qcd_mc_event_weight);
    inputTree->SetBranchAddress("jet_CamKt12Truth_pt", &qcd_CA12_truth_pt);
    inputTree->SetBranchAddress("jet_CamKt12Truth_eta", &qcd_CA12_truth_eta);
    inputTree->SetBranchAddress("jet_CamKt12Truth_phi", &qcd_CA12_truth_phi);
    inputTree->SetBranchAddress("jet_CamKt12Truth_m", &qcd_CA12_truth_mass);
    inputTree->SetBranchAddress("jet_CamKt12Truth_E", &qcd_CA12_truth_E);
    
    inputTree->SetBranchAddress("jet_CamKt12LCTopo_pt", &qcd_CA12_topo_pt);
    inputTree->SetBranchAddress("jet_CamKt12LCTopo_eta", &qcd_CA12_topo_eta);
    inputTree->SetBranchAddress("jet_CamKt12LCTopo_phi", &qcd_CA12_topo_phi);
    inputTree->SetBranchAddress("jet_CamKt12LCTopo_m", &qcd_CA12_topo_mass);
    inputTree->SetBranchAddress("jet_CamKt12LCTopo_E", &qcd_CA12_topo_E);
    inputTree->SetBranchAddress("jet_CamKt12LCTopo_emfrac", &qcd_CA12_topo_emfrac);
    
    inputTree->SetBranchAddress("jet_CamKt12LC"+groomAlgo+"_pt", &qcd_CA12_groomed_pt);
    inputTree->SetBranchAddress("jet_CamKt12LC"+groomAlgo+"_eta", &qcd_CA12_groomed_eta);
    inputTree->SetBranchAddress("jet_CamKt12LC"+groomAlgo+"_phi", &qcd_CA12_groomed_phi);
    inputTree->SetBranchAddress("jet_CamKt12LC"+groomAlgo+"_m", &qcd_CA12_groomed_mass);
    inputTree->SetBranchAddress("jet_CamKt12LC"+groomAlgo+"_E", &qcd_CA12_groomed_E);
    inputTree->SetBranchAddress("jet_CamKt12LC"+groomAlgo+"_emfrac", &qcd_CA12_groomed_emfrac);
    
    inputTree1->SetBranchAddress("mc_channel_number", &Wp_mc_channel_number);
    inputTree1->SetBranchAddress("mc_event_weight", &Wp_mc_event_weight);
    inputTree1->SetBranchAddress("jet_CamKt12Truth_pt", &Wp_CA12_truth_pt);
    inputTree1->SetBranchAddress("jet_CamKt12Truth_eta", &Wp_CA12_truth_eta);
    inputTree1->SetBranchAddress("jet_CamKt12Truth_phi", &Wp_CA12_truth_phi);
    inputTree1->SetBranchAddress("jet_CamKt12Truth_m", &Wp_CA12_truth_mass);
    inputTree1->SetBranchAddress("jet_CamKt12Truth_E", &Wp_CA12_truth_E);
    
    inputTree1->SetBranchAddress("jet_CamKt12LCTopo_pt", &Wp_CA12_topo_pt);
    inputTree1->SetBranchAddress("jet_CamKt12LCTopo_eta", &Wp_CA12_topo_eta);
    inputTree1->SetBranchAddress("jet_CamKt12LCTopo_phi", &Wp_CA12_topo_phi);
    inputTree1->SetBranchAddress("jet_CamKt12LCTopo_m", &Wp_CA12_topo_mass);
    inputTree1->SetBranchAddress("jet_CamKt12LCTopo_E", &Wp_CA12_topo_E);
    inputTree1->SetBranchAddress("jet_CamKt12LCTopo_emfrac", &Wp_CA12_topo_emfrac);
    
    inputTree1->SetBranchAddress("jet_CamKt12LC"+groomAlgo+"_pt", &Wp_CA12_groomed_pt);
    inputTree1->SetBranchAddress("jet_CamKt12LC"+groomAlgo+"_eta", &Wp_CA12_groomed_eta);
    inputTree1->SetBranchAddress("jet_CamKt12LC"+groomAlgo+"_phi", &Wp_CA12_groomed_phi);
    inputTree1->SetBranchAddress("jet_CamKt12LC"+groomAlgo+"_m", &Wp_CA12_groomed_mass);
    inputTree1->SetBranchAddress("jet_CamKt12LC"+groomAlgo+"_E", &Wp_CA12_groomed_E);
    inputTree1->SetBranchAddress("jet_CamKt12LC"+groomAlgo+"_emfrac", &Wp_CA12_groomed_emfrac);
  }
 



  //if reclustering algorithm
  if (groomAlgoIndex>6){


    inputTree->SetBranchAddress("mc_channel_number", &qcd_mc_channel_number);
    inputTree->SetBranchAddress("mc_event_weight", &qcd_mc_event_weight);
    inputTree->SetBranchAddress("jet_CamKt12Truth_pt", &qcd_CA12_truth_pt);
    inputTree->SetBranchAddress("jet_CamKt12Truth_eta", &qcd_CA12_truth_eta);
    inputTree->SetBranchAddress("jet_CamKt12Truth_phi", &qcd_CA12_truth_phi);
    inputTree->SetBranchAddress("jet_CamKt12Truth_m", &qcd_CA12_truth_mass);
    inputTree->SetBranchAddress("jet_CamKt12Truth_E", &qcd_CA12_truth_E);

    /*  inputTree->SetBranchAddress("mc_channel_number", &qcd_mc_channel_number);
    inputTree->SetBranchAddress("mc_event_weight", &qcd_mc_event_weight);
    inputTree->SetBranchAddress("jet_AntiKt10Truth_pt", &qcd_CA12_truth_pt);
    inputTree->SetBranchAddress("jet_AntiKt10Truth_eta", &qcd_CA12_truth_eta);
    inputTree->SetBranchAddress("jet_AntiKt10Truth_phi", &qcd_CA12_truth_phi);
    inputTree->SetBranchAddress("jet_AntiKt10Truth_m", &qcd_CA12_truth_mass);
    inputTree->SetBranchAddress("jet_AntiKt10Truth_E", &qcd_CA12_truth_E);
    */
    //AntiKt2LCTopo

    inputTree->SetBranchAddress("jet_"+groomAlgo+"_pt", &qcd_CA12_topo_pt);
    inputTree->SetBranchAddress("jet_"+groomAlgo+"_eta", &qcd_CA12_topo_eta);
    inputTree->SetBranchAddress("jet_"+groomAlgo+"_phi", &qcd_CA12_topo_phi);
    inputTree->SetBranchAddress("jet_"+groomAlgo+"_m", &qcd_CA12_topo_mass);
    inputTree->SetBranchAddress("jet_"+groomAlgo+"_E", &qcd_CA12_topo_E);
    inputTree->SetBranchAddress("jet_"+groomAlgo+"_emfrac", &qcd_CA12_topo_emfrac);
    
    /*inputTree->SetBranchAddress("jet_CamKt12LC"+groomAlgo+"_pt", &qcd_CA12_groomed_pt);
    inputTree->SetBranchAddress("jet_CamKt12LC"+groomAlgo+"_eta", &qcd_CA12_groomed_eta);
    inputTree->SetBranchAddress("jet_CamKt12LC"+groomAlgo+"_phi", &qcd_CA12_groomed_phi);
    inputTree->SetBranchAddress("jet_CamKt12LC"+groomAlgo+"_m", &qcd_CA12_groomed_mass);
    inputTree->SetBranchAddress("jet_CamKt12LC"+groomAlgo+"_E", &qcd_CA12_groomed_E);
    inputTree->SetBranchAddress("jet_CamKt12LC"+groomAlgo+"_emfrac", &qcd_CA12_groomed_emfrac);
    */
 
    inputTree1->SetBranchAddress("mc_channel_number", &Wp_mc_channel_number);
    inputTree1->SetBranchAddress("mc_event_weight", &Wp_mc_event_weight);
    inputTree1->SetBranchAddress("jet_CamKt12Truth_pt", &Wp_CA12_truth_pt);
    inputTree1->SetBranchAddress("jet_CamKt12Truth_eta", &Wp_CA12_truth_eta);
    inputTree1->SetBranchAddress("jet_CamKt12Truth_phi", &Wp_CA12_truth_phi);
    inputTree1->SetBranchAddress("jet_CamKt12Truth_m", &Wp_CA12_truth_mass);
    inputTree1->SetBranchAddress("jet_CamKt12Truth_E", &Wp_CA12_truth_E);
    
   /*
    inputTree1->SetBranchAddress("mc_channel_number", &Wp_mc_channel_number);
    //inputTree1->SetBranchAddress("mc_event_weight", &Wp_mc_event_weight);
    inputTree1->SetBranchAddress("jet_AntiKt10Truth_pt", &Wp_CA12_truth_pt);
    inputTree1->SetBranchAddress("jet_AntiKt10Truth_eta", &Wp_CA12_truth_eta);
    inputTree1->SetBranchAddress("jet_AntiKt10Truth_phi", &Wp_CA12_truth_phi);
    inputTree1->SetBranchAddress("jet_AntiKt10Truth_m", &Wp_CA12_truth_mass);
    inputTree1->SetBranchAddress("jet_AntiKt10Truth_E", &Wp_CA12_truth_E);
    */
    inputTree1->SetBranchAddress("jet_"+groomAlgo+"_pt", &Wp_CA12_topo_pt);
    inputTree1->SetBranchAddress("jet_"+groomAlgo+"_eta", &Wp_CA12_topo_eta);
    inputTree1->SetBranchAddress("jet_"+groomAlgo+"_phi", &Wp_CA12_topo_phi);
    inputTree1->SetBranchAddress("jet_"+groomAlgo+"_m", &Wp_CA12_topo_mass);
    inputTree1->SetBranchAddress("jet_"+groomAlgo+"_E", &Wp_CA12_topo_E);
    inputTree1->SetBranchAddress("jet_"+groomAlgo+"_emfrac", &Wp_CA12_topo_emfrac);
    /*
    inputTree1->SetBranchAddress("jet_CamKt12LC"+groomAlgo+"_pt", &Wp_CA12_groomed_pt);
    inputTree1->SetBranchAddress("jet_CamKt12LC"+groomAlgo+"_eta", &Wp_CA12_groomed_eta);
    inputTree1->SetBranchAddress("jet_CamKt12LC"+groomAlgo+"_phi", &Wp_CA12_groomed_phi);
    inputTree1->SetBranchAddress("jet_CamKt12LC"+groomAlgo+"_m", &Wp_CA12_groomed_mass);
    inputTree1->SetBranchAddress("jet_CamKt12LC"+groomAlgo+"_E", &Wp_CA12_groomed_E);
    inputTree1->SetBranchAddress("jet_CamKt12LC"+groomAlgo+"_emfrac", &Wp_CA12_groomed_emfrac);
    */

  }


 
}




void makeROC(int type, TH1D *&S, TH1D *&B, TGraph &curve, TGraph &curve_up, TGraph &curve_do){
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

  curve=gr;
  curve_up=gr_up;
  curve_do=gr_do;

  return;

}
