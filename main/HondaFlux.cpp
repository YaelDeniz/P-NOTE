#define _USE_MATH_DEFINES
#include <cmath>

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

//Cern ROOT
#include "TApplication.h"
#include "TRootCanvas.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "Math/Interpolator.h"
//My Libraries
#include "MathTools.h"
#include "PhyTools.h"
#include "GetTrueEvents.h"


# define myPi 3.14159265358979323846  /* pi */
// Some Constants i need
# define R_earth 6368.0 //km
# define mN   1.67492749804E-27  // Nucleons mass ~ Neutron mass
# define MTon  1E9  //Metric MegaTon
# define years2sec 3.154E7 // Years in Seconds

using namespace std;

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include "TH2D.h"



void ExportHistogramToCSV(TH2D* hist, const std::string& filename = "histogram2d.csv") {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Could not open file for writing!\n";
        return;
    }

    out << "x_bin_low,x_bin_high,y_bin_low,y_bin_high,count\n";

    for (int ix = 1; ix <= hist->GetNbinsX(); ++ix) {
        double xlow = hist->GetXaxis()->GetBinLowEdge(ix);
        double xhigh = hist->GetXaxis()->GetBinUpEdge(ix);

        for (int iy = 1; iy <= hist->GetNbinsY(); ++iy) {
            double ylow = hist->GetYaxis()->GetBinLowEdge(iy);
            double yhigh = hist->GetYaxis()->GetBinUpEdge(iy);
            double count = hist->GetBinContent(ix, iy);

            out << std::fixed << std::setprecision(6)
                << xlow << "," << xhigh << ","
                << ylow << "," << yhigh << ","
                << count << "\n";
        }
    }

    out.close();
    std::cout << "Histogram exported to " << filename << "\n";
}

//------------------
int main(int argc, char **argv)
{
 
    NuFlux HondaFlux;
    //SouthPole
    std::string SPFolder = "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation/NuFlux/";
    std::string SPFile = "SP_AziAveraged_solmin/spl-nu-20-01-000.d";
    std::string SPTable = SPFolder + SPFile; //Class Assumes Sout Pole flux as default

    //Gran Sasso
    std::string GRNFolder = "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation/NuFlux/";
    std::string GRNFile = "GRN_AziAveraged_solmin/grn-nu-20-01-000.d";
    std::string GRNTable = GRNFolder + GRNFile; //Class Assumes Sout Pole flux as default
 
    std::string FluxTable = SPTable; 
    
    std::vector< std::vector<double> > FluxData = HondaFlux.SetFluxData(FluxTable);

    //Matrix for Histogram & Histogram Draw

    TH2D* muflux =  HondaFlux.GetFluxHist(1,FluxData); //MuFlux
    TH2D* mubflux =  HondaFlux.GetFluxHist(2,FluxData); //MuBarFlux
    TH2D* eflux =  HondaFlux.GetFluxHist(3,FluxData); //EFlux
    TH2D* ebflux =  HondaFlux.GetFluxHist(4,FluxData); //EBarFlux

//---------------------------------------------------------------------------------------------------------------------------------------
muflux->SetStats(0);
muflux-> SetTitle("#nu_{#mu} Flux");
muflux->GetXaxis()->SetTitle("log_{10}(E/GeV)");
muflux->GetYaxis()->SetTitle("Cos(#theta_{z})");
//muflux->GetZaxis()->SetTitle("log_{10}(#Phi_{#nu})");

mubflux->SetStats(0);
mubflux-> SetTitle("#bar{#nu}_{#mu} Flux");
mubflux->GetXaxis()->SetTitle("log_{10}(E/GeV)");
mubflux->GetYaxis()->SetTitle("Cos(#theta_{z})");
//mubflux->GetZaxis()->SetTitle("log_{10}(#Phi_{#nu})");

eflux->SetStats(0);
eflux-> SetTitle("#nu_{e} Flux");
eflux->GetXaxis()->SetTitle("log_{10}(E/GeV)");
eflux->GetYaxis()->SetTitle("Cos(#theta_{z})");
//eflux->GetZaxis()->SetTitle("log_{10}(#Phi_{#nu})");

ebflux->SetStats(0);
ebflux-> SetTitle("#bar{#nu}_{e} Flux");
ebflux->GetXaxis()->SetTitle("log_{10}(E/GeV)");
ebflux->GetYaxis()->SetTitle("Cos(#theta_{z})");
//ebflux->GetZaxis()->SetTitle("log_{10}(#Phi_{#nu})");




//gStyle->SetPadLeftMargin(10000); 
//gStyle->SetPadRightMargin(0.1);




//Interpolation of the flux





    double Emax = 100;

    double Emin = 1;

    double pts = 1000; 

    double hE = (100 - 1)/pts;

    double Eo = 1;

    

    TGraph  *numu = new TGraph(pts);
    TGraph  *numub = new TGraph(pts);
    TGraph  *nue = new TGraph(pts);
    TGraph  *nueb = new TGraph(pts);

 std::cout << "----------------" << std::endl;

    for (int i = 0; i < pts; ++i)
     {
          

 

         double Ei =  Eo + i*hE;

         double logEi = log10(Ei);
         
         double logphimu = muflux->Interpolate(logEi,-0.9);
         double logphimub = mubflux->Interpolate(logEi,-0.9);
         double logphie = eflux->Interpolate(logEi,-0.9);
         double logphieb = ebflux->Interpolate(logEi,-0.9);

         

         double phimu = pow(10,logphimu)*Ei*Ei*Ei;
         double phimub = pow(10,logphimub)*Ei*Ei*Ei;
         double phie = pow(10,logphie)*Ei*Ei*Ei;
         double phieb = pow(10,logphieb)*Ei*Ei*Ei;

         numu -> SetPoint (i,  logEi,  phimu);
         numub -> SetPoint (i,  logEi,  phimub);
         nue -> SetPoint (i,  logEi,  phie);
         nueb -> SetPoint (i,  logEi,  phieb);
     } 


    ExportHistogramToCSV(eflux, "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation/HondaFlux_hists/eflux_hist.csv");
    ExportHistogramToCSV(ebflux, "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation/HondaFlux_hists/ebflux_hist.csv");
    ExportHistogramToCSV(muflux, "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation/HondaFlux_hists/muflux_hist.csv");
    ExportHistogramToCSV(mubflux, "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation/HondaFlux_hists/mubflux_hist.csv");
   
    TMultiGraph * multi = new TMultiGraph();

    numu->SetLineColor(4);
    nue->SetLineColor(2);
    numub->SetLineColor(4);
    nueb->SetLineColor(2);

    numub->SetLineStyle(9);
    nueb->SetLineStyle(9);

    multi->Add(numu);
    multi->Add(numub);
    multi->Add(nue);
    multi->Add(nueb);



    multi->GetXaxis()->SetLimits(0.0,2.0);
    multi->GetXaxis()->SetTitle("log_{10}(E/GeV)");

    multi->GetYaxis()->SetLimits(0.0,0.03);
    multi->GetYaxis()->SetTitle("#Phi_{#nu}*E^{3}");

    multi->SetTitle("Bilinear interpolation of Honda flux");

   
    TApplication app("app", &argc, argv);    


    TCanvas *c = new TCanvas();
    c->Divide(2,2);
    c->cd(1);
    eflux->Draw("COLZ");
    c->cd(2);
    ebflux->Draw("COLZ");
    c->cd(3);
    muflux->Draw("COLZ");
    c->cd(4);
    mubflux->Draw("COLZ");
    gPad->Update();

    c->Modified(); c->Update();
    TRootCanvas *rc = (TRootCanvas *)c->GetCanvasImp();
    rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");

   
    
    TCanvas *c1 = new TCanvas();
    gPad->SetGrid(1,1);
    multi->Draw("A");
    gPad->Update();

    c1->Modified(); c1->Update();
    TRootCanvas *rc1 = (TRootCanvas *)c1->GetCanvasImp();
    rc1->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
    
    app.Run();



    c->SaveAs("/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation/SimulationResults/NeutrinoFlux/FluxHistSP.jpg");
    c1->SaveAs("/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation/SimulationResults/NeutrinoFlux/FluxInterpolSP.jpg");

   return 0;
   
}

/*
// Writes one CSV with columns:
// log10E, cosThetaZ, nu_mu, nu_mu_bar, nu_e, nu_e_bar
// (bin *centers*; excludes under/overflow)
void write_flux_csv(
    const std::string& csv_path,
    const TH2D* muflux,
    const TH2D* mubflux,
    const TH2D* eflux,
    const TH2D* ebflux)
{
    // (Optional) basic sanity checks: same binning for all
    auto same_binning = [&](const TH2D* a, const TH2D* b) {
        return a->GetNbinsX() == b->GetNbinsX() &&
               a->GetNbinsY() == b->GetNbinsY() &&
               std::equal(a->GetXaxis()->GetXbins()->GetArray(),
                          a->GetXaxis()->GetXbins()->GetArray() + (a->GetXaxis()->GetNbins()+1),
                          b->GetXaxis()->GetXbins()->GetArray(),
                          [](double lhs, double rhs){ return std::abs(lhs-rhs) < 1e-12; }) &&
               std::equal(a->GetYaxis()->GetXbins()->GetArray(),
                          a->GetYaxis()->GetXbins()->GetArray() + (a->GetYaxis()->GetNbins()+1),
                          b->GetYaxis()->GetXbins()->GetArray(),
                          [](double lhs, double rhs){ return std::abs(lhs-rhs) < 1e-12; });
    };

    // If variable binning wasn’t used, GetXbins()->GetArray() is null; fall back to counts only
    auto compatible_counts_only = [&](const TH2D* a, const TH2D* b) {
        return a->GetNbinsX() == b->GetNbinsX() && a->GetNbinsY() == b->GetNbinsY();
    };

    bool grids_ok = compatible_counts_only(muflux, mubflux) &&
                    compatible_counts_only(muflux, eflux)   &&
                    compatible_counts_only(muflux, ebflux);

    if (!grids_ok) {
        std::cerr << "[write_flux_csv] Histograms do not have matching bin counts. "
                     "Cannot write a single table reliably.\n";
        return;
    }

    std::ofstream out(csv_path);
    if (!out) {
        std::cerr << "[write_flux_csv] Could not open: " << csv_path << "\n";
        return;
    }

    out << std::setprecision(10) << std::fixed;
    out << "log10E,cosThetaZ,nu_mu,nu_mu_bar,nu_e,nu_e_bar\n";

    const int nx = muflux->GetNbinsX();
    const int ny = muflux->GetNbinsY();

    for (int ix = 1; ix <= nx; ++ix) {
        const double x = muflux->GetXaxis()->GetBinCenter(ix); // log10(E/GeV)
        for (int iy = 1; iy <= ny; ++iy) {
            const double y = muflux->GetYaxis()->GetBinCenter(iy); // cos(theta_z)

            const double v_mu   = muflux->GetBinContent(ix, iy);
            const double v_mub  = mubflux->GetBinContent(ix, iy);
            const double v_e    = eflux->GetBinContent(ix, iy);
            const double v_eb   = ebflux->GetBinContent(ix, iy);

            out << x << ',' << y << ','
                << v_mu  << ',' << v_mub << ','
                << v_e   << ',' << v_eb  << '\n';
        }
    }
    out.close();
    std::cout << "[write_flux_csv] Wrote: " << csv_path << "\n";
}
*/