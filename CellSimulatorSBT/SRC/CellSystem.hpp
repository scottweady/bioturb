#ifndef CELLSYSTEM_HPP_
#define CELLSYSTEM_HPP_

#include "Config.h"

#include "HydroRodMobility.hpp"
#include "STKFMM/Lib/include/STKFMM/STKFMM.hpp"

#include "SimToolbox/Sylinder/SylinderSystem.hpp"
#include "SimToolbox/Util/QuadInt.hpp"

#include <memory>
#include <vector>

#include <mpi.h>
#include <omp.h>

template <int N>
class CellSystem {

    SylinderSystem<N> rodSystem;
    Config cellConfig;
    Teuchos::RCP<const TCOMM> commRcp;
    Teuchos::RCP<HydroRodMobility<PS::ParticleSystem<Sylinder<N>>>> hydroMobOpRcp;

    // FMM stuff
    std::shared_ptr<stkfmm::STKFMM> fmmPtr;
    std::vector<double> srcCoord;
    std::vector<double> srcValue;
    std::vector<double> trgCoord;
    std::vector<double> trgValue;

    using Quad = QuadInt<N>;
    std::vector<Quad> quads;

    void prepareQuad();

    void prepareRandomSpecies(int speciesNum);

    void initHydro();

    void calcHydroVelocity();

    void calcNonBrownVelocity();

    void initCellGrowth();

    void calcCellGrowth();

    void calcCellDivision();

    void writeFlowGridVTI(std::string prefix);

    void setFMMBox();

    void configCheck() const;

  public:
    CellSystem(const std::string &runConfig, const std::string &cellConfig, const std::string &posFile,
               const std::string &restartFile, int argc, char **argv);

    void output();

    void step();

    bool stop();
};

#include "CellSystem.tpp"
#endif