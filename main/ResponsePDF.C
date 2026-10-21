
#include "TH2D.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TMath.h"

//FisherZenith(ThetaReco, ThetaTrue, kappa);

double PDFcth(double threco, double thtrue, double kappa);
double PDFth(double threco, double thtrue, double sigmath);

void ResponsePDF()
{ 

  double threcomin = (180-180)*TMath::Pi()/180.0;
  double threcomax = (180-140)*TMath::Pi()/180.0;

  double n = 1000;

  double dth = (threcomax-threcomin)/n;

  double thtrue = (180-169)*TMath::Pi()/180.0;

  //0.04 + 0.18/(sqrt(e))
  //double Ath    = (2)*TMath::Pi()/180.0;
  //double Bth    = (10)*TMath::Pi()/180.0;

  double Ath    = 0.04;
  double Bth    = 0.18;


  double Etrue  = 10;

  double sigmath = Ath + Bth/sqrt(Etrue);
  double kappa = 1/(sigmath*sigmath);

  TGraph * plotcth = new TGraph(n);
  plotcth->SetLineColor(kBlue);
  plotcth->SetLineWidth(4);
  plotcth->SetTitle("VMF");


  TGraph * plotth  = new TGraph(n);
  plotth->SetLineColor(kRed);
  plotth->SetLineStyle(10);
  plotth->SetLineWidth(4);
  plotth->SetTitle("Gauss");

  TGraph * tho = new TGraph(2);
  tho->SetPoint (0, (180.0*thtrue)/TMath::Pi(), 0);
  tho->SetPoint (1, (180.0*thtrue)/TMath::Pi(), 4.5 );




  TMultiGraph *mg = new TMultiGraph();



  for (int i = 0; i <=n; ++i)
  {
    
  double threco = threcomin + i*dth;
  double pdfcth = PDFcth(threco, thtrue, kappa );
  double pdfth = PDFth(threco,thtrue, sigmath );

  plotcth->SetPoint (i, 180.0*threco/TMath::Pi(), pdfcth );
  plotth ->SetPoint (i, 180.0*threco/TMath::Pi(), pdfth );
  std::cout <<  Ath<< " " << Bth << std::endl;
  std::cout << "reference " <<threco<< " " <<thtrue << std::endl;

  }


  TCanvas *c = new TCanvas();
  mg->Add(plotcth,"lp");
  mg->Add(plotth,"cp");
  mg->Add(tho,"cp");
  mg->Draw("a");

}

double PDFcth(double threco, double thtrue, double kappa){

    /// @param sigma - std. dev. in radians. Angular resolution
    /// kappa = 1/sigma^2 - VMF concentrarion parameter. 

  double cos_z = cos(threco);
  double cos_z0 = cos(thtrue);

  double sin_z = sqrt(1 - pow(cos_z, 2));
  double sin_z0 = sqrt(1 - pow(cos_z0, 2));

  double exp_arg = kappa * cos_z * cos_z0;

  double norm = kappa;
  if (kappa < 50)
    norm /= 2 * sinh(kappa);
  else
    exp_arg -= kappa;

  double arg = kappa * sin_z * sin_z0;
  double bess = 1;
  if (arg < 50)
    bess = TMath::BesselI0(arg);
  else {
    exp_arg += arg;
    bess *=
        (1 + 1. / (8 * arg) + 4.5 / pow(8 * arg, 2) + 37.5 / pow(8 * arg, 3));
    bess /= sqrt(2 * M_PI * arg);
  }

  double expo = exp(exp_arg);

  std::cout << "VMF Values x:" << threco << " sigma:" << sqrt(1/kappa) <<  " x0:" << thtrue << std::endl; 

  return norm * expo * bess * sin_z;

}

double PDFth(double threco, double thtrue, double sigmath){

  std::cout << "Gaussian Values x:" << threco << " sigma:" << sigmath <<  " x0:" << thtrue << std::endl; 

  return ROOT::Math::gaussian_pdf   (  threco, sigmath, thtrue  );   



}