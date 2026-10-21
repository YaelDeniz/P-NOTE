#include "ExportData.h"
#include <fstream>
#include "TH2D.h"   // actual ROOT header

void SenvHistToCSV(TH2D* stdhist,
                 TH2D* althist,
                 TH2D* hist,
                 const std::string& filename)
{

    std::ofstream outfile(filename);

    for (int i = 1; i <= hist->GetNbinsX(); ++i)
    {
        for (int j = 1; j <= hist->GetNbinsY(); ++j)
        {
            double x = hist->GetXaxis()->GetBinCenter(i);
            double y = hist->GetYaxis()->GetBinCenter(j);

            double nexp = stdhist->GetBinContent(i, j);
            double nobs = althist->GetBinContent(i, j);
            double diff = hist->GetBinContent(i, j);

            outfile << x << "," << y << "," << nexp
                    << " , " << nobs
                    << " , " << diff << '\n';
        }
    }
}

