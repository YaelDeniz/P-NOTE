// LayerUpdate.h
#pragma once

struct LayerUpdate
{
    int    tableNumber;      // 1-based PREM CSV row
    double densityContrast;  // relative change in density (%)
    double chemicalContrast; // relative change in Z/A (%)
};
