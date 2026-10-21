#pragma once

#include <string>
#include "TH2D.h"   // actual ROOT header

void SenvHistToCSV(TH2D* stdhist,
                 TH2D* althist,
                 TH2D* hist,
                 const std::string& filename);

