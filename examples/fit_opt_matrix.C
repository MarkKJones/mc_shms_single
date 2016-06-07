#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <TFile.h>
#include <TNtuple.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TMath.h>
#include <TBox.h>
#include <TPolyLine.h>
#include <TLegend.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TDecompSVD.h>
void fit_opt_matrix() {
  gROOT->Reset();
  gStyle->SetOptStat(0);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.05,"XY");
 gStyle->SetPadLeftMargin(0.17);
 //
 TH1F *hDelta = new TH1F("hDelta","Delta ",20,-.1,.1);
 TH1F *hDeltarecon = new TH1F("hDeltarecon","Delta Recon % ",40,-.12,.12);
 TH1F *hDeltadiff = new TH1F("hDeltadiff","Delta Diff % ",40,-.1,.1);
 TH1F *hDeltanew = new TH1F("hDeltanew","Delta New Recon % ",40,-.12,.12);
 TH1F *hDeltanewdiff = new TH1F("hDeltanewdiff","Delta New Diff % ",40,-.1,.1);
 TH1F *hytar = new TH1F("hytar","ytar (cm)",70,-7.,7.);
 TH1F *hytarrecon = new TH1F("hytarrecon","ytar recon(cm)",70,-7.,7.);
 TH1F *hytarnew = new TH1F("hytarnew","ytar new(cm)",70,-7.,7.);
 TH1F *hytardiff = new TH1F("hytardiff","ytar diff(cm)",70,-2.,2.);
 TH1F *hytarnewdiff = new TH1F("hytarnewdiff","ytar new diff(cm)",70,-2.,2.);
 TH1F *hxptar = new TH1F("hxptar","xptar ",100,-.12,.12);
 TH1F *hxptarrecon = new TH1F("hxptarrecon","xptar recon",100,-.12,.12);
 TH1F *hxptardiff = new TH1F("hxptardiff","xptar diff (mr)",100,-10,10);
 TH1F *hxptarnew = new TH1F("hxptarnew","xptar new recon",100,-.12,.12);
 TH1F *hxptarnewdiff = new TH1F("hxptarnewdiff","xptar new diff (mr)",100,-10,10);
 TH1F *hyptar = new TH1F("hyptar","yptar ",100,-.1,.1);
 TH1F *hyptarrecon = new TH1F("hyptarrecon","yptar recon ",100,-.1,.1);
 TH1F *hyptardiff = new TH1F("hyptardiff","yptar diff (mr) ",100,-10,10);
 TH1F *hyptarnew = new TH1F("hyptarnew","yptar new recon ",100,-.1,.1);
 TH1F *hyptarnewdiff = new TH1F("hyptarnewdiff","yptar new diff (mr) ",100,-10,10);
  //
   string newcoeffsfilename="newfit.dat";
    string oldcoeffsfilename="shms/shms_recon_fit_90deg_1cm_5th_order.dat";
  ifstream oldcoeffsfile(oldcoeffsfilename.c_str());
  ofstream newcoeffsfile(newcoeffsfilename.c_str());

  vector<double> xptarcoeffs_old;
  vector<double> yptarcoeffs_old;
  vector<double> ytarcoeffs_old;
  vector<double> deltacoeffs_old;
  vector<int> xfpexpon_old;
  vector<int> xpfpexpon_old;
  vector<int> yfpexpon_old;
  vector<int> ypfpexpon_old;
  vector<int> xtarexpon_old;

  vector<double> xptarcoeffs_fit;
  vector<double> yptarcoeffs_fit;
  vector<double> ytarcoeffs_fit;
  vector<double> deltacoeffs_fit;
  vector<int> xfpexpon_fit;
  vector<int> xpfpexpon_fit;
  vector<int> yfpexpon_fit;
  vector<int> ypfpexpon_fit;
  vector<int> xtarexpon_fit;

  vector<double> xptarcoeffs_new;
  vector<double> yptarcoeffs_new;
  vector<double> ytarcoeffs_new;
  vector<double> deltacoeffs_new;
  vector<int> xfpexpon_new;
  vector<int> xpfpexpon_new;
  vector<int> yfpexpon_new;
  vector<int> ypfpexpon_new;
  vector<int> xtarexpon_new;
  vector<double> xtartrue,ytartrue,xptartrue,yptartrue,deltatrue;
  vector<double> xfptrue,yfptrue,xpfptrue,ypfptrue;
  TString currentline;
  int num_recon_terms_old;
  int num_recon_terms_fit;

  num_recon_terms_old = 0;
  num_recon_terms_fit = 0;
  int nfit=0,npar,nfit_max=5000,npar_final=0,max_order=6,norder;

  while( currentline.ReadLine(oldcoeffsfile,kFALSE) && !currentline.BeginsWith(" ----") ){

    TString sc1(currentline(1,16));
    TString sc2(currentline(17,16));
    TString sc3(currentline(33,16));
    TString sc4(currentline(49,16));
    
    xptarcoeffs_old.push_back(sc1.Atof());
    ytarcoeffs_old.push_back(sc2.Atof());
    yptarcoeffs_old.push_back(sc3.Atof());
    deltacoeffs_old.push_back(sc4.Atof());
    
    int expontemp[5];

    for(int expon=0; expon<5; expon++){
      TString stemp(currentline(66+expon,1));
      expontemp[expon] = stemp.Atoi();
    }

  

    xfpexpon_old.push_back(expontemp[0]);
    xpfpexpon_old.push_back(expontemp[1]);
    yfpexpon_old.push_back(expontemp[2]);
    ypfpexpon_old.push_back(expontemp[3]);
    xtarexpon_old.push_back(expontemp[4]);

   
    
    num_recon_terms_old++;
    norder= expontemp[0]+expontemp[1]+expontemp[2]+expontemp[3]+expontemp[4];
    if (norder <= max_order) {
    xptarcoeffs_fit.push_back(sc1.Atof());
    ytarcoeffs_fit.push_back(sc2.Atof());
    yptarcoeffs_fit.push_back(sc3.Atof());
    deltacoeffs_fit.push_back(sc4.Atof());
    xfpexpon_fit.push_back(expontemp[0]);
    xpfpexpon_fit.push_back(expontemp[1]);
    yfpexpon_fit.push_back(expontemp[2]);
    ypfpexpon_fit.push_back(expontemp[3]);
    xtarexpon_fit.push_back(expontemp[4]);
    num_recon_terms_fit++;
    }
  }

  cout << "num recon terms in OLD matrix = " << num_recon_terms_old << endl;
  cout << "num recon terms in fit matrix = " << num_recon_terms_fit << endl;
   npar= num_recon_terms_fit ;
   //
  TVectorD b_ytar(npar);
  TVectorD b_yptar(npar);
  TVectorD b_xptar(npar);
  TVectorD b_delta(npar);
  TMatrixD lambda(npar,nfit_max);
  TMatrixD Ay(npar,npar);
  //
   /*
   TVectorD *b_ytar;
  TVectorD *b_yptar;
  TVectorD *b_xptar;
  TVectorD *b_delta;
  TMatrixD *lambda;
  TMatrixD *Ay;
  b_ytar = new TVectorD(npar);
  b_yptar = new TVectorD(npar);
  b_xptar = new TVectorD(npar);
  b_delta = new TVectorD(npar);
  lambda = new TMatrixD(npar,nfit_max);
   Ay = new TMatrixD(npar,nfit_max);
  */
 //
  //   TString basename = "shms_pointtarg_9p5deg_2p2gev";
     TString basename = "shms_20cmtarg_2gev_nowc_nomscat";
   TString inputroot="rootfiles/"+basename+".root";
   TFile *fsimc = new TFile(inputroot);
   TTree *tsimc = (TTree*)fsimc->Get("h1411");
//Declaration of leaves types
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
   //
   Double_t xfp,yfp,xpfp,ypfp,xtar;
   Double_t xptar,yptar,ytar,delta;
   Long64_t nentries = tsimc->GetEntries();
   for (int i = 0; i < nentries; i++) {
        tsimc->GetEntry(i);
	xfp = hsxfp;
	yfp = hsyfp;
	xpfp = hsxpfp;
	ypfp = hsypfp;
        xtar=hsxtari;
	xptar=hsxptari;
	yptar=hsyptari;
	ytar=hsytari;
	delta=hsdeltai/100.;
	//
	if (TMath::Abs(delta) < .5 ) {
	hDelta->Fill(hsdelta);
	hytar->Fill(hsytar);
	hyptar->Fill(hsyptar);
	hxptar->Fill(hsxptar);
	hDeltadiff->Fill(hsdelta-hsdeltai);
	hytardiff->Fill(hsytar-hsytari);
	hyptardiff->Fill(hsyptar-hsyptari);
	hxptardiff->Fill(hsxptar-hsxptari);
          Double_t ytartemp = 0.0,yptartemp=0.0,xptartemp=0.0,deltatemp=0.0;
          Double_t etemp;
          for( int icoeffold=0; icoeffold<num_recon_terms_old; icoeffold++ ){
        	etemp= 
	  pow( xfp / 100.0, xfpexpon_old[icoeffold] ) * 
	  pow( yfp / 100.0, yfpexpon_old[icoeffold] ) * 
	  pow( xpfp, xpfpexpon_old[icoeffold] ) * 
	  pow( ypfp, ypfpexpon_old[icoeffold] ) * 
	  pow( xtar/100., xtarexpon_old[icoeffold] );
        	deltatemp += deltacoeffs_old[icoeffold] * etemp;
        	ytartemp += ytarcoeffs_old[icoeffold] * etemp;
	        yptartemp += yptarcoeffs_old[icoeffold] * etemp;
	         xptartemp += xptarcoeffs_old[icoeffold] *etemp; 
	  } // for icoeffold loop
	  hDeltarecon->Fill(deltatemp*100.);
 	hytarrecon->Fill(ytartemp*100.);
	hyptarrecon->Fill(yptartemp);
	hxptarrecon->Fill(xptartemp);
          for( int icoeff_fit=0; icoeff_fit<num_recon_terms_fit; icoeff_fit++ ){
        	etemp= 
	  pow( xfp / 100.0, xfpexpon_fit[icoeff_fit] ) * 
	  pow( yfp / 100.0, yfpexpon_fit[icoeff_fit] ) * 
	  pow( xpfp, xpfpexpon_fit[icoeff_fit] ) * 
	  pow( ypfp, ypfpexpon_fit[icoeff_fit] ) * 
	  pow( xtar/100., xtarexpon_fit[icoeff_fit] );
 		 if (nfit < nfit_max ) {
              lambda[icoeff_fit][nfit] = etemp;
	      b_xptar[icoeff_fit] += (xptar) * etemp;
	      b_yptar[icoeff_fit] += (yptar) * etemp;
	      b_ytar[icoeff_fit] += (ytar) /100.0 * etemp;
	      b_delta[icoeff_fit] += (delta) * etemp;
              }
	  } // for icoeff_fit loop
	      if (nfit < nfit_max ) {
	      nfit++;
  	    xfptrue.push_back( xfp );
	    yfptrue.push_back( yfp );
	    xpfptrue.push_back( xpfp );
	    ypfptrue.push_back( ypfp );
	    xtartrue.push_back( xtar );
	    xptartrue.push_back( xptar );
	    ytartrue.push_back( ytar  );
	    yptartrue.push_back( yptar  );
	    deltatrue.push_back( delta  );
	      }
        }
   }
   //
   //
   cout << " number to fit = " << nfit << " max = " << nfit_max << endl;
 for(int i=0; i<npar; i++){
    for(int j=0; j<npar; j++){
      Ay[i][j] = 0.0;
    }
 }
   for( int ifit=0; ifit<nfit; ifit++){
    if( ifit % 100 == 0 ) cout << ifit << endl;
    for( int ipar=0; ipar<npar; ipar++){
      for( int jpar=0; jpar<npar; jpar++){
      	Ay[ipar][jpar] += lambda[ipar][ifit] * lambda[jpar][ifit];
       }
    }
  }
   //
   TDecompSVD Ay_svd(Ay);
  bool ok;
  ok = Ay_svd.Solve( b_ytar );
  cout << "ytar solution ok = " << ok << endl;
  //b_ytar.Print();
  ok = Ay_svd.Solve( b_yptar );
  cout << "yptar solution ok = " << ok << endl;
  //b_yptar.Print();
  ok = Ay_svd.Solve( b_xptar );
  cout << "xptar solution ok = " << ok << endl;
  //b_xptar.Print();
  ok = Ay_svd.Solve( b_delta );
  cout << "Delta solution ok = " << ok << endl;
  //b_delta.Print();
  // calculate target quantities with new fit parameter
  for( int ifit=0; ifit<nfit; ifit++){
          Double_t ytarnew = 0.0,yptarnew=0.0,xptarnew=0.0,deltanew=0.0;
	  Double_t etemp;
     for( int ipar=0; ipar<npar; ipar++){
       etemp=lambda[ipar][ifit];
        	deltanew += b_delta[ipar] * etemp;
        	ytarnew += b_ytar[ipar] * etemp;
	        yptarnew += b_yptar[ipar] * etemp;
	         xptarnew += b_xptar[ipar] *etemp;        
    }
	  hytarnew->Fill(ytarnew*100.);
	  hyptarnew->Fill(yptarnew);
	  hxptarnew->Fill(xptarnew);
	  hDeltanew->Fill(deltanew*100);
	  hytarnewdiff->Fill(ytarnew*100.-ytartrue.at(ifit));
	  hyptarnewdiff->Fill(1000*(yptarnew-yptartrue.at(ifit)));
	  hxptarnewdiff->Fill(1000*(xptarnew-xptartrue.at(ifit)));
	  hDeltanewdiff->Fill(100*(deltanew-deltatrue.at(ifit)));    
  }
  // write out coeff
  char coeffstring[100];
  Double_t tt;
  cout << "writing new coeffs file" << endl;
          for( int icoeff_fit=0; icoeff_fit<num_recon_terms_fit; icoeff_fit++ ){
      newcoeffsfile << " ";
      tt=b_xptar[icoeff_fit];
      sprintf( coeffstring, "%16.9g", tt );
      newcoeffsfile << coeffstring; 
      //      newcoeffsfile << " ";
      sprintf( coeffstring, "%16.9g", b_ytar[icoeff_fit] );
      newcoeffsfile << coeffstring;
      sprintf( coeffstring, "%16.9g", b_yptar[icoeff_fit] );
      //newcoeffsfile << " ";
      newcoeffsfile << coeffstring; 
      sprintf( coeffstring, "%16.9g", b_delta[icoeff_fit] );
      //newcoeffsfile << " ";
      newcoeffsfile << coeffstring; 
      newcoeffsfile << " ";
	newcoeffsfile << setw(1) << setprecision(1) << xfpexpon_fit[icoeff_fit]; 
	newcoeffsfile << setw(1) << setprecision(1) << xpfpexpon_fit[icoeff_fit]; 
	newcoeffsfile << setw(1) << setprecision(1) << yfpexpon_fit[icoeff_fit]; 
	newcoeffsfile << setw(1) << setprecision(1) << ypfpexpon_fit[icoeff_fit]; 
	newcoeffsfile << setw(1) << setprecision(1) << xtarexpon_fit[icoeff_fit]; 
      newcoeffsfile << endl;

	  }
  newcoeffsfile << " ---------------------------------------------" << endl;

  newcoeffsfile.close();
  cout << "wrote new coeffs file" << endl;
//
    TCanvas *cdiff = new TCanvas("cdiff","Old Diff target",800,800);
    cdiff->Divide(2,2);
    cdiff->cd(1);
    hytardiff->Draw();
 hytardiff->Fit("gaus");
 TF1 *fitcydiff=hytardiff->GetFunction("gaus");
    cdiff->cd(2);
    hyptardiff->Draw();
 hyptardiff->Fit("gaus");
 TF1 *fitcypdiff=hyptardiff->GetFunction("gaus");
    hyptarnewdiff->Draw("same");
    cdiff->cd(3);
    hxptardiff->Draw();
   hxptardiff->Fit("gaus");
 TF1 *fitcxpdiff=hxptardiff->GetFunction("gaus");
    cdiff->cd(4);
    hDeltadiff->Draw();
 hDeltadiff->Fit("gaus");
 TF1 *fitcdeldiff=hDeltadiff->GetFunction("gaus");
//
    TCanvas *cnewdiff = new TCanvas("cnewdiff","Old Newdiff target",800,800);
    cnewdiff->Divide(2,2);
    cnewdiff->cd(1);
    hytarnewdiff->Draw();
 hytarnewdiff->Fit("gaus");
 TF1 *fitcynewdiff=hytarnewdiff->GetFunction("gaus");
    cnewdiff->cd(2);
    hyptarnewdiff->Draw();
 hyptarnewdiff->Fit("gaus");
 TF1 *fitcypnewdiff=hyptarnewdiff->GetFunction("gaus");
    cnewdiff->cd(3);
    hxptarnewdiff->Draw();
     hxptarnewdiff->Fit("gaus");
 TF1 *fitcxpnewdiff=hxptarnewdiff->GetFunction("gaus");
    cnewdiff->cd(4);
    hDeltanewdiff->Draw();
 hDeltanewdiff->Fit("gaus");
 TF1 *fitcdelnewdiff=hDeltanewdiff->GetFunction("gaus");
  //
}
