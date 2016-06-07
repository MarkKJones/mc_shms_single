#include <iostream>
#include <TString.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TMath.h>
#include <TProfile.h>
#include <TObjArray.h>
void plot_carbon_elastic(TString inputroot) {
 gROOT->Reset();
gStyle->SetPalette(1,0);
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.12);
 //   inputroot="rootfiles/"+basename+".root";
TFile *fsimc = new TFile(inputroot);
TTree *tsimc = (TTree*)fsimc->Get("h1411");
 Float_t         hsxfp; // position at focal plane ,+X is pointing down
 Float_t         hsyfp; // X x Y = Z so +Y pointing central ray left
 Float_t         hsxpfp; // dx/dz at focal plane
 Float_t         hsypfp; //  dy/dz at focal plane
 Float_t         hsztari; // thrown position along the beam direction
 Float_t         hsytari;  //thrown  horizontal position X x Y = Z so +Y pointing central ray left at plane perpendicular to SHMS at z=0
 Float_t         hsdeltai; // thrown  100*(p - pc)/pc with pc = central SHMS momentum 
 Float_t         hsyptari; // thrown target dy/dz horizontal slope
 Float_t         hsxptari; // thrown target dx/dz vertical slope
 Float_t         hsztar; // reconstructed position along the beam directio
 Float_t         hsytar; //reconstructed horizontal position
 Float_t         hsdelta;//reconstructed
   Float_t         hsyptar;//reconstructed
   Float_t         hsxptar;//reconstructed
   Float_t         hsxtari;// from thrown kinematics , the calculated vertical position at plane perpendicular to SHMS at z=0
   Float_t         yrast;// vertical raster position
   Float_t         xsnum;
   Float_t         ysnum;
   Float_t         xsieve;
   Float_t         ysieve;
   Float_t         evtype;
   Float_t         normfac;
   Float_t         theta;
   Float_t         eprime;
   Float_t         ecalc;
   Float_t         Q2;
   Float_t         W;
   Float_t         xsn;
   // Set branch addresses.
   tsimc->SetBranchAddress("hsxfp",&hsxfp);
   tsimc->SetBranchAddress("hsyfp",&hsyfp);
   tsimc->SetBranchAddress("hsxpfp",&hsxpfp);
   tsimc->SetBranchAddress("hsypfp",&hsypfp);
   tsimc->SetBranchAddress("hsztari",&hsztari);
   tsimc->SetBranchAddress("hsytari",&hsytari);
   tsimc->SetBranchAddress("hsdeltai",&hsdeltai);
   tsimc->SetBranchAddress("hsyptari",&hsyptari);
   tsimc->SetBranchAddress("hsxptari",&hsxptari);
   tsimc->SetBranchAddress("hsztar",&hsztar);
   tsimc->SetBranchAddress("hsytar",&hsytar);
   tsimc->SetBranchAddress("hsdelta",&hsdelta);
   tsimc->SetBranchAddress("hsyptar",&hsyptar);
   tsimc->SetBranchAddress("hsxptar",&hsxptar);
   tsimc->SetBranchAddress("hsxtari",&hsxtari);
   tsimc->SetBranchAddress("yrast",&yrast);
   tsimc->SetBranchAddress("xsnum",&xsnum);
   tsimc->SetBranchAddress("ysnum",&ysnum);
   tsimc->SetBranchAddress("xsieve",&xsieve);
   tsimc->SetBranchAddress("ysieve",&ysieve);
   tsimc->SetBranchAddress("evtype",&evtype);
   tsimc->SetBranchAddress("normfac",&normfac);
   tsimc->SetBranchAddress("xsn",&xsn);
   tsimc->SetBranchAddress("theta",&theta);
   tsimc->SetBranchAddress("eprime",&eprime);
   tsimc->SetBranchAddress("ecalc",&ecalc);
   tsimc->SetBranchAddress("Q2",&Q2);
   tsimc->SetBranchAddress("W",&W);
   //;ysieve
 TH1F *hsdelta_all = new TH1F("hsdelta_all",";hsdelta ",100,-15.,25.);
 TH1F *hq2_ev3 = new TH1F("hq2_ev3","; Q2 (GeV^2) ",100,0.,2.5);
 TH1F *hw2_ev3 = new TH1F("hw_ev3","; W (GeV) ",100,.5,3.0);
 TH1F *hsdeltadiff_all = new TH1F("hsdeltadiff_all","; (hsdelta-hsdeltai) ",100,-.1,.1);
 TH1F *hsdeltadiff_ev1 = new TH1F("hsdeltadiff_ev1","; (hsdelta-hsdeltai) ",100,-.1,.1);
 TH1F *hsdeltadiff_ev3 = new TH1F("hsdeltadiff_ev3","; (hsdelta-hsdeltai) ",100,-.1,.1);
 TH1F *hsdeltadiff_ev1_xsn = new TH1F("hsdeltadiff_ev1_xsn","; (hsdelta-hsdeltai) ",100,-.1,.1);
 TH1F *hsdeltadiff_ev3_xsn = new TH1F("hsdeltadiff_ev3_xsn","; (hsdelta-hsdeltai) ",100,-.1,.1);
 TH2F *hsdeltaVSdiff_all = new TH2F("hsdeltaVSdiff_all","; hsdelta : (hsdelta-hsdeltai) ",100,-15.,20.,100,-.1,.1);
 TH1F *hsdelta_ev1 = new TH1F("hsdelta_ev1",";hsdelta ",100,-15.,25.);
 TH1F *hsdelta_ev1_xsn = new TH1F("hsdelta_ev1_xsn",";hsdelta ",100,-5.,5.);
 TH1F *hsdelta_ev3_xsn = new TH1F("hsdelta_ev3_xsn",";hsdelta ",100,-5.,5.);
 TH1F *hsdelta_ev3 = new TH1F("hsdelta_ev3",";hsdelta ",100,-15.,25.);
 TH1F *heprime_all = new TH1F("heprime_all","; eprime (GeV) ",100,1.5,2.5);
 TH1F *heprime_ev1 = new TH1F("heprime_ev1","; eprime (GeV) ",100,1.5,2.5);
 TH1F *heprime_ev3 = new TH1F("heprime_ev3","; eprime (GeV) ",100,1.5,2.5);
 TH1F *hediff_all = new TH1F("hediff_all","; (eprime - ecalc_theta)/eprime (%) ",100,-.1,.1);
 TH1F *hediff_ev1 = new TH1F("hediff_ev1","; (eprime - ecalc_theta)/eprime (%) ",100,-.1,.1);
 TH1F *hediff_ev1_xsn = new TH1F("hediff_ev1_xsn","; (eprime - ecalc_theta)/eprime (%) ",100,-.1,.1);
 TH1F *hediff_ev3_xsn = new TH1F("hediff_ev3_xsn","; (eprime - ecalc_theta)/eprime (%) ",100,-.1,.1);
 TH1F *hediff_ev3 = new TH1F("hediff_ev3","; (eprime - ecalc_theta)/eprime (%) ",100,-.1,.1);
 TH2F *heprVStheta_all = new TH2F("heprVStheta_all",";theta (deg); eprime (GeV) ",50,8.,11.,50,1.5,2.5);
 TH2F *heprVStheta_ev1 = new TH2F("heprVStheta_ev1",";theta (deg); eprime (GeV) ",50,8.,11.,50,1.5,2.5);
 TH2F *heprVStheta_ev3 = new TH2F("heprVStheta_ev3",";theta (deg); eprime (GeV) ",50,8.,11.,50,1.5,2.5);
 TH2F *hxsnVStheta_ev1 = new TH2F("hxsnVStheta_ev1",";theta (deg); xsn ",50,8.,11.,50,0.,.0001);
 TH2F *hxsnVStheta_ev3 = new TH2F("hxsnVStheta_ev3",";theta (deg); xsn ",50,8.,11.,50,0.,.0000001);
 //
 int numxsh,numysh;
   Long64_t nentries = tsimc->GetEntries();
   for (int i = 0; i < nentries; i++) {
        tsimc->GetEntry(i);
	hsdelta_all->Fill(hsdelta);
	hsdeltadiff_all->Fill((hsdelta-hsdeltai));
	if (evtype ==3 ) hq2_ev3->Fill(Q2,normfac);
	if (evtype ==3 ) hw_ev3->Fill(W,normfac);
	if (evtype ==1 ) hsdeltadiff_ev1->Fill(hsdelta-hsdeltai);
	if (evtype ==3 ) hsdeltadiff_ev3->Fill(hsdelta-hsdeltai);
	if (evtype ==1 ) hsdeltadiff_ev1_xsn->Fill(hsdelta-hsdeltai,normfac);
	if (evtype ==3 ) hsdeltadiff_ev3_xsn->Fill(hsdelta-hsdeltai,normfac);
	hsdeltaVSdiff_all->Fill(hsdeltai,(hsdelta-hsdeltai));
	if (evtype ==1 ) hsdelta_ev1->Fill(hsdelta);
	if (evtype ==3 ) hsdelta_ev3->Fill(hsdelta);
	if (evtype ==1 ) hsdelta_ev1_xsn->Fill(hsdelta,normfac);
	if (evtype ==3 ) hsdelta_ev3_xsn->Fill(hsdelta,normfac);
	heprime_all->Fill(eprime/1000.);
	if (evtype ==1 ) heprime_ev1->Fill(eprime/1000.);
	if (evtype ==3 ) heprime_ev3->Fill(eprime/1000.);
	heprVStheta_all->Fill(theta*180/3.14159,eprime/1000.);
	if (evtype ==1 ) hxsnVStheta_ev1->Fill(theta*180/3.14159,xsn);
	if (evtype ==3 ) hxsnVStheta_ev3->Fill(theta*180/3.14159,xsn);
	if (evtype ==1 ) heprVStheta_ev1->Fill(theta*180/3.14159,eprime/1000.);
	if (evtype ==3 ) heprVStheta_ev3->Fill(theta*180/3.14159,eprime/1000.);
	hediff_all->Fill(100*(eprime-ecalc)/eprime);
	if (evtype ==1 ) hediff_ev1->Fill(100*(eprime-ecalc)/eprime);
	if (evtype ==1 ) hediff_ev1_xsn->Fill(100*(eprime-ecalc)/eprime,normfac);
	if (evtype ==3 ) hediff_ev3_xsn->Fill(100*(eprime-ecalc)/eprime,normfac);
	if (evtype ==3 ) hediff_ev3->Fill(100*(eprime-ecalc)/eprime);
   }
   //
TCanvas *cfp = new TCanvas("cfp","focal plane",1400,900);
cfp->Divide(2,2);
cfp->cd(1);
 hsdelta_all->Draw("");
 hsdelta_ev1->Draw("same");
 hsdelta_ev1->SetLineColor(2);
 hsdelta_ev3->Draw("same");
 hsdelta_ev3->SetLineColor(3);
cfp->cd(2);
 heprime_all->Draw("");
 heprime_ev1->Draw("same");
 heprime_ev1->SetLineColor(2);
 heprime_ev3->Draw("same");
 heprime_ev3->SetLineColor(3);
cfp->cd(3);
 hediff_all->Draw("");
 hediff_ev1->Draw("same");
 hediff_ev1->SetLineColor(2);
 hediff_ev3->Draw("same");
 hediff_ev3->SetLineColor(3);
cfp->cd(4);
 hsdeltadiff_all->Draw("");
 hsdeltadiff_ev1->Draw("same");
 hsdeltadiff_ev1->SetLineColor(2);
 hsdeltadiff_ev3->Draw("same");
 hsdeltadiff_ev3->SetLineColor(3);
 //
TCanvas *cfp1 = new TCanvas("cfp1","Eprime versus theta",1400,900);
cfp1->Divide(2,2);
cfp1->cd(1);
 heprVStheta_all->Draw("colz") ;
cfp1->cd(2);
 heprVStheta_ev1->Draw("colz") ;
cfp1->cd(3);
 heprVStheta_ev3->Draw("colz") ;
cfp1->cd(4);
 hsdeltaVSdiff_all->Draw("colz") ;
 //
TCanvas *cfp2 = new TCanvas("cfp2","Kinematics",1400,900);
cfp2->Divide(2,2);
cfp2->cd(1);
 hq2_ev3->Draw();
cfp2->cd(2);
 hw_ev3->Draw();
 cout << " Integral rate over all W = " << hw_ev3->Integral() << endl;
 cfp2->cd(3);
   hxsnVStheta_ev1->Draw("colz");
 cfp2->cd(4);
   hxsnVStheta_ev3->Draw("colz");
 //
TCanvas *cfp3 = new TCanvas("cfp3","Rate",1400,900);
cfp3->Divide(2,2);
cfp3->cd(1);
 hsdelta_ev1_xsn->Draw();
 cout << " rate for 20uA = " << hsdelta_ev1_xsn->Integral() << endl;
cfp3->cd(2);
 hsdelta_ev3_xsn->Draw();
 cout << " rate for 20uA = " << hsdelta_ev3_xsn->Integral() << endl;
cfp3->cd(3);
 hediff_ev1_xsn->Draw();
 cout << " rate for 20uA = " << hediff_ev1_xsn->Integral() << endl;
cfp3->cd(4);
 hediff_ev3_xsn->Draw();
 cout << " rate for 20uA = " << hediff_ev3_xsn->Integral() << endl;
//
}
