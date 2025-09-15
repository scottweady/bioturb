#ifndef CONFIG_H_
#define CONFIG_H_

#include "GrowthTable.h"
#include "HydroTable.h"

#include <iostream>
#include <memory>
#include <unordered_map>

struct Growth {
    double tauD = 0;
    double Delta = 0;
    double sigma = 0;
    GrowthTable cellGrowthTable;
    Growth() = delete;
    Growth(double time, double length, double s, const GrowthTable &table) : cellGrowthTable(table) {
        tauD = time;
        Delta = length;
        sigma = s;
    }
    void echo() const;
};

struct Hydro {
    int numberOfQuadraturePoints = 16;
    double forceExternal[3] = {0.0, 0.0, 0.0};
    double torqueExternal[3] = {0.0, 0.0, 0.0};
    HydroTable cellHydroTable;
    Hydro() = delete;
    Hydro(int numQuadPt, const double (&Fext)[3], const double (&Text)[3], const HydroTable &table)
        : cellHydroTable(table) {
        numberOfQuadraturePoints = numQuadPt;
        std::copy(std::begin(Fext), std::end(Fext), std::begin(forceExternal));
        std::copy(std::begin(Text), std::end(Text), std::begin(torqueExternal));
    }
    void echo() const;
};

class Config {
  public:
    // hydro swim setting
    bool hydro = false;                             ///< false = turn off hydro
    double freeStreamVelocity[3] = {0.0, 0.0, 0.0}; ///< global background velocity
    int fmmMultOrder = 10;                          ///< fmm multipole order, should be between [6,16]
    double flowGrid = -1.0;                         ///< negative for no flow dump
    bool wallZ0;

    std::unordered_map<int, Hydro> cellHydro;   ///< hydro settings for all species
    std::unordered_map<int, Growth> cellGrowth; ///< division settings for all species

    Config(std::string);
    ~Config() = default;

    void echo() const;
};

#endif
