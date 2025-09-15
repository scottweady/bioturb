#ifndef HYDRORODOPERATOR_HPP_
#define HYDRORODOPERATOR_HPP_

#include "Config.h"

#include <STKFMM/STKFMM.hpp>

#include "Sylinder/SylinderConfig.hpp"
#include "Trilinos/Preconditioner.hpp"
#include "Trilinos/TpetraUtil.hpp"
#include "Util/EigenDef.hpp"
#include "Util/QuadInt.hpp"
#include "Util/SpecialQuadWeights.hpp"
#include "SimToolbox/Util/IOHelper.hpp"

struct HydroOption {
    bool withSwim = false;
    bool withForceTorqueExternal = false;
};

// /****************************************************************************
//  *                                                                          *
//  * Weighted RPY Force,a Vel,lapVel kernel, source: 4, target: 6             *
//  *       fx,fy,fz,a -> ux,uy,uz,lapux,lapuy,lapuz                           *
//  ****************************************************************************/
void rpy_weighted_ulapu(const int nsrc, const double *src_coord, const double *src_value, const double *trg_coord,
                        double *trg_value, const double *w1_all, const double *w3_all, const double *w5_all,
                        double scale) {
    // This kernel is for a single target point - line source pair with known quadrature weights w1, w3, w5 related to
    // the rinv, rinv3, and rinv5 terms, respectively. use scale to control the bounds of integration: scale=1
    // corresponds to an integral from [-1,1]
    constexpr double Pi = M_PI;
    const double FACV = 1.0 / (8 * Pi);
    const double tx = trg_coord[0];
    const double ty = trg_coord[1];
    const double tz = trg_coord[2];

    double vx = 0.0;
    double vy = 0.0;
    double vz = 0.0;
    double lapvx = 0.0;
    double lapvy = 0.0;
    double lapvz = 0.0;

    for (int s = 0; s < nsrc; s++) {
        const double dx = tx - src_coord[3 * s + 0];
        const double dy = ty - src_coord[3 * s + 1];
        const double dz = tz - src_coord[3 * s + 2];

        const double fx = src_value[4 * s + 0];
        const double fy = src_value[4 * s + 1];
        const double fz = src_value[4 * s + 2];
        const double a = src_value[4 * s + 3];
        const double w1 = w1_all[s] * scale;
        const double w3 = w3_all[s] * scale;
        const double w5 = w5_all[s] * scale;

        const double a2_over_three = (1.0 / 3.0) * a * a;
        const double r2 = dx * dx + dy * dy + dz * dz;

        const double rinv = 1.0 / std::sqrt(r2);
        const double rinv3 = rinv * rinv * rinv;
        const double rinv5 = rinv * rinv * rinv3;
        const double w1_rinv = w1 * rinv;
        const double w3_rinv3 = w3 * rinv3;
        const double w5_rinv5 = w5 * rinv5;
        const double fdotr = fx * dx + fy * dy + fz * dz;

        const double three_fdotr_w5_rinv5 = 3 * fdotr * w5_rinv5;
        const double cx = fx * w3_rinv3 - three_fdotr_w5_rinv5 * dx;
        const double cy = fy * w3_rinv3 - three_fdotr_w5_rinv5 * dy;
        const double cz = fz * w3_rinv3 - three_fdotr_w5_rinv5 * dz;

        const double fdotr_w3_rinv3 = fdotr * w3_rinv3;
        vx = vx + fx * w1_rinv + dx * fdotr_w3_rinv3 + a2_over_three * cx;
        vy = vy + fy * w1_rinv + dy * fdotr_w3_rinv3 + a2_over_three * cy;
        vz = vz + fz * w1_rinv + dz * fdotr_w3_rinv3 + a2_over_three * cz;

        lapvx = lapvx + 2 * cx;
        lapvy = lapvy + 2 * cy;
        lapvz = lapvz + 2 * cz;
    }

    vx = vx * FACV + trg_value[0];
    vy = vy * FACV + trg_value[1];
    vz = vz * FACV + trg_value[2];
    lapvx = lapvx * FACV + trg_value[3];
    lapvy = lapvy * FACV + trg_value[4];
    lapvz = lapvz * FACV + trg_value[5];

    trg_value[0] = vx;
    trg_value[1] = vy;
    trg_value[2] = vz;
    trg_value[3] = lapvx;
    trg_value[4] = lapvy;
    trg_value[5] = lapvz;
}

/***
 * terminology:
 * upper case F/T for force/torque per rod
 * lower case f for force density per qudrature point
 * upper case Vel/Omega for Velocity/Omega per rod
 * lower case u for velocity per quadrature point
 *
 * for force/torque per rod
 * 'external'    : force/torque input in CellConfig
 * 'input'   : force/torque input by function parameters
 * for background velocity
 * 'freestreamvelocity' : input in CellConfig
 * for unknown force density
 * 'e': hydro force excluding the fs(s)p
 * hydro force density
 * 'fh': hydro force density, can be either fe, or fs or both
 */

template <class Container>
class HydroRodOperator : public TOP {
  private:
    const Container *const containerPtr;    ///< read-only
    const int nRodLocal;                    ///< local number of rods
    std::shared_ptr<stkfmm::STKFMM> fmmPtr; ///< pointer to fmm, either 3d or above wall
    const Config &cellConfig;               ///< hydro user input
    const SylinderConfig &runConfig;        ///< general user input

    // indices, maps, etc, precomputed in constructor
    // Q: number of quadrature points per rod
    Teuchos::RCP<const TCOMM> commRcp;
    Teuchos::RCP<TMAP> rodMapRcp;         ///< 1 dof per rod
    Teuchos::RCP<TMAP> rodPtsMapRcp;      ///< Q dof per rod
    Teuchos::RCP<TMAP> pointValuesMapRcp; ///< 3Q dof per rod
    std::vector<int> rodPts;              ///< 1 dof per rod, stores Q
    std::vector<int> rodPtsIndex;         ///< beginning of quadrature points per rod in rodPtsMap

    // these are precomputed in constructor
    std::vector<double> FTExt; ///< external force/torque specified in HydroTable, 6 dof per rod
    std::vector<double> Hu;    ///< <us(s),1> integral, 1 dof per rod
    std::vector<double> Hf;    ///< <fs(s),1> integral, 1 dof per rod

    // these two are original [-1,1] data, not scaled by actual length
    std::vector<double> sloc;         ///< Q dof per rod, quadrature points in ([-1,1])
    std::vector<double> weight;       ///< Q dof per rod, quadrature weights
    std::vector<double> slipVelocity; ///< Q dof per rod, slipVelocity for all rods
    std::vector<double> activeStress; ///< Q dof per rod, activeStress for all rods
    std::vector<double> radius;       ///< Q dof per rod, radius for all rods
    Teuchos::RCP<TV> fsRcp;           ///< 3Q dof per rod, fs(s)p of each point
    Teuchos::RCP<TV> fhRcp;           ///< 3Q dof per rod, fs(s)p + fe of each point

    // temporary data for uinf and their integrals
    mutable std::vector<double> uinf;   ///< 3Q dof per rod, uinf due to fs(s)p along other rods
    mutable std::vector<double> Oinf;   ///< 3 dof per rod, integral <uinf,1>
    mutable std::vector<double> Sinf;   ///< 3 dof per rod, integral <uinf,s>
    mutable std::vector<double> uinfOS; ///< 3Q dof per rod, 1/c(Imhalfpp)uinf-1/2cl(Imhalfpp)Oinf-3s/2cl^3(Impp)S

    // these are computed at every runFMM(), thus mutable
    mutable std::vector<double> srcSLCoord; ///< fmm data, 3*number of quadrature points per rod
    mutable std::vector<double> trgCoord;   ///< fmm data, 3*number of quadrature points per rod
    mutable std::vector<double> srcSLValue; ///< fmm data, 4*number of quadrature points per rod
    mutable std::vector<double> trgValue;   ///< fmm data, 6*number of quadrature points per rod

    void setupDOF();

    // initialize FMM tree
    void setupFMM();

    // calc Oinf, Sinf, and 1/c(Imhalfpp)uinf-1/2cl(Imhalfpp)Oinf-3s/2cl^3(Impp)Sinf
    void calcuinfOS() const;

  public:
    /************************************************
     * Interface required by TOP
     *
     *
     *
     ***********************************************/
    HydroRodOperator(const Container *const containerPtr_, const int nRodLocal_,
                     std::shared_ptr<stkfmm::STKFMM> &fmmPtr_, const Config &cellConfig_,
                     const SylinderConfig &runConfig_);

    ~HydroRodOperator() = default;

    Teuchos::RCP<const TMAP> getDomainMap() const { return pointValuesMapRcp; }

    Teuchos::RCP<const TMAP> getRangeMap() const { return pointValuesMapRcp; }

    bool hasTransposeApply() const { return false; }

    // Compute Y := alpha Op X + beta Y.
    void apply(const TMV &X, TMV &Y, Teuchos::ETransp mode = Teuchos::NO_TRANS,
               scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
               scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const;

    /***********************************************
     * Interface for hydro calculation
     *
     *
     *
     ***********************************************/

    void writeFlowGridVTI(const std::string &prefix, const double &flowGrid, const bool &wallz0, const double& viscosity, const int &fileID);

    // calc uinf generated by given hydro force at quadrature points
    void runFMM(const TV &fhVec) const;

    // temporary O(N^2) test of the mixed quadrature method. Call in place of runFMM.
    void testMixedQuadrature(const TV &fhVec) const;

    // calculate Afe with arbitrary given fe on quadrature points
    void calcAfe(const TV &feVec, TV &AfeVec) const;

    // calculate velocity/omega with arbitrary given fe on quadrature points
    void calcVelOmega(const TV &feVec, TV &VelOmegaVec, const HydroOption &options) const;

    // calculate rightside vector b with arbitrary given force/torque on each rod
    void calcb(const TV &FTinputVec, TV &bVec, const HydroOption &options) const;

    // return a TV of total fh on quad points. fh = fs + fe
    Teuchos::RCP<TV> getfhRcp() const { return fhRcp; };

    // return a TV of fs on quad points.
    Teuchos::RCP<TV> getfsRcp() const { return fsRcp; };

    const std::vector<double> &getuinf() const { return uinf; };

    const std::vector<double> &getOinf() const { return Oinf; };

    const std::vector<double> &getSinf() const { return Sinf; };

    const std::vector<int> &getPts() const { return rodPts; };

    const std::vector<int> &getPtsIndex() const { return rodPtsIndex; };
};

template <class Container>
HydroRodOperator<Container>::HydroRodOperator(const Container *const containerPtr_, const int nRodLocal_,
                                              std::shared_ptr<stkfmm::STKFMM> &fmmPtr_, const Config &cellConfig_,
                                              const SylinderConfig &runConfig_)
    : containerPtr(containerPtr_), nRodLocal(nRodLocal_), fmmPtr(fmmPtr_), cellConfig(cellConfig_),
      runConfig(runConfig_) {
    commRcp = getMPIWORLDTCOMM();
    setupDOF();
    setupFMM();
}

template <class Container>
void HydroRodOperator<Container>::setupDOF() {
    // task: initialize all data structures independent of fe data
    rodPts.resize(nRodLocal);
    rodPtsIndex.resize(nRodLocal + 1, 0);
    const auto &rodContainer = *containerPtr;
    for (int i = 0; i < nRodLocal; i++) {
        rodPts[i] = rodContainer[i].numQuadPt;
        rodPtsIndex[i + 1] = rodPtsIndex[i] + rodPts[i];
    }

    rodMapRcp = getTMAPFromLocalSize(nRodLocal, commRcp);
    rodPtsMapRcp = getTMAPFromLocalSize(rodPtsIndex.back(), commRcp);
    pointValuesMapRcp = getTMAPFromLocalSize(3 * rodPtsIndex.back(), commRcp);

    const int nQuadLocal = rodPtsIndex.back();

    // preallocate spaces
    FTExt.resize(6 * nRodLocal, 0); ///< external force/torque specified in HydroTable, 6 dof per rod
    Hu.resize(nRodLocal, 0);        ///< <us(s),1> integral, 1 dof per rod
    Hf.resize(nRodLocal, 0);        ///< <fs(s),1> integral, 1 dof per rod

    // these two are original [-1,1] data, not scaled by actual length
    sloc.resize(nQuadLocal, 0);         ///< Q dof per rod, quadrature points in ([-1,1])
    weight.resize(nQuadLocal, 0);       ///< Q dof per rod, quadrature weights
    slipVelocity.resize(nQuadLocal, 0); ///< Q dof per rod, slipVelocity for all rods
    activeStress.resize(nQuadLocal, 0); ///< Q dof per rod, activeStress for all rods
    radius.resize(nQuadLocal, 0);       ///< Q dof per rod, radius for all rods
    fsRcp = Teuchos::rcp(new TV(pointValuesMapRcp, true));
    fhRcp = Teuchos::rcp(new TV(pointValuesMapRcp, true));
    auto fsPtr = fsRcp->getLocalView<Kokkos::HostSpace>(Tpetra::Access::ReadWrite);
    auto fhPtr = fhRcp->getLocalView<Kokkos::HostSpace>(Tpetra::Access::ReadWrite);

    // temporary data for uinf and their integrals
    uinf.resize(3 * nQuadLocal, 0); ///< 3Q dof per rod,
    Oinf.resize(3 * nQuadLocal, 0); ///< 3 dof per rod,
    Sinf.resize(3 * nQuadLocal, 0); ///< 3 dof per rod,

#pragma omp parallel for
    for (size_t i = 0; i < nRodLocal; i++) {
        const auto &sy = rodContainer[i];
        const double length = sy.length;
        const int numQuadPt = sy.numQuadPt;
        const auto quadPtr = sy.quadPtr;
        assert(numQuadPt == rodPts[i]);
        assert(sy.numQuadPt == quadPtr->getSize());

        const Evec3 direction = ECmapq(sy.orientation) * Evec3(0, 0, 1);
        const double *sQuadPt = quadPtr->getPoints();
        const double *weightQuadPt = quadPtr->getWeights();

        std::vector<double> slipVelocityPt(numQuadPt, 0);
        std::vector<double> activeStressPt(numQuadPt, 0);
        std::vector<double> radiusPt(numQuadPt, 0);
        // setup data for quadrature points according to cell species
        const Hydro *pHydro = &(cellConfig.cellHydro.at(sy.speciesID)); 
        const auto &table = pHydro->cellHydroTable;
        table.getValue(numQuadPt, quadPtr->getPoints(), slipVelocityPt.data(), activeStressPt.data(), radiusPt.data());

        Hu[i] = length * 0.5 * quadPtr->intSamples(slipVelocityPt.data());
        Hf[i] = length * 0.5 * quadPtr->intSamples(activeStressPt.data());
        FTExt[6 * i + 0] = pHydro->forceExternal[0];
        FTExt[6 * i + 1] = pHydro->forceExternal[1];
        FTExt[6 * i + 2] = pHydro->forceExternal[2];
        FTExt[6 * i + 3] = pHydro->torqueExternal[0];
        FTExt[6 * i + 4] = pHydro->torqueExternal[1];
        FTExt[6 * i + 5] = pHydro->torqueExternal[2];

        const int idx = rodPtsIndex[i];
        for (size_t j = 0; j < numQuadPt; j++) {
            sloc[idx + j] = sQuadPt[j];
            weight[idx + j] = weightQuadPt[j];
            slipVelocity[idx + j] = slipVelocityPt[j];
            activeStress[idx + j] = activeStressPt[j];
            radius[idx + j] = radiusPt[j];
            fsPtr(3 * (idx + j) + 0, 0) = activeStressPt[j] * direction[0];
            fsPtr(3 * (idx + j) + 1, 0) = activeStressPt[j] * direction[1];
            fsPtr(3 * (idx + j) + 2, 0) = activeStressPt[j] * direction[2];
            fhPtr(3 * (idx + j) + 0, 0) = sy.forceHydro[3 * j + 0];
            fhPtr(3 * (idx + j) + 1, 0) = sy.forceHydro[3 * j + 1];
            fhPtr(3 * (idx + j) + 2, 0) = sy.forceHydro[3 * j + 2];
        }
    }

    // for (auto &v : slipVelocity) {
    //     printf("%g\n", v);
    // }
}

template <class Container>
void HydroRodOperator<Container>::setupFMM() {
    // task: get point coordinates and setup fmm tree
    using namespace stkfmm;
    const auto &rodContainer = *containerPtr;
    const int nPtsLocal = rodPtsIndex.back(); // number of quadrature points on local
    const int nSrc = nPtsLocal;
    const int nTrg = nPtsLocal;
    srcSLCoord.clear();
    trgCoord.clear();
    uinf.clear();
    srcSLCoord.resize(3 * nSrc, 0);
    trgCoord.resize(3 * nTrg, 0);

#pragma omp parallel for
    for (int i = 0; i < nRodLocal; i++) {
        const auto &sy = rodContainer[i];
        const double length = sy.length;
        const int numQuadPt = sy.numQuadPt;
        assert(sy.numQuadPt == rodPts[i]);

        const Evec3 direction = ECmapq(sy.orientation) * Evec3(0, 0, 1);
        const int idx = rodPtsIndex[i];
        for (int iQuadPt = 0; iQuadPt < numQuadPt; iQuadPt++) {
            Evec3 loc = ECmap3(sy.pos) + (length * 0.5 * sloc[idx + iQuadPt]) * direction;
            srcSLCoord[3 * (idx + iQuadPt) + 0] = loc[0];
            srcSLCoord[3 * (idx + iQuadPt) + 1] = loc[1];
            srcSLCoord[3 * (idx + iQuadPt) + 2] = loc[2];
        }
    }

    trgCoord = srcSLCoord;

    fmmPtr->setPoints(nSrc, srcSLCoord.data(), nTrg, trgCoord.data());
    fmmPtr->setupTree(KERNEL::RPY);
}

template <class Container>
void HydroRodOperator<Container>::runFMM(const TV &fhVec) const {
    // task: compute uinf with given force density fhVec.
    using namespace stkfmm;
    if (!fmmPtr) {
        printf("error: fmm not initialized.\n");
    }

    assert(fhVec.getLocalLength() == 3 * this->rodPtsIndex.back());
    assert(fhVec.getLocalLength() == this->pointValuesMapRcp->getLocalNumElements());
    auto fhPtr = fhVec.getLocalView<Kokkos::HostSpace>(Tpetra::Access::ReadOnly);

    // the FMM bounding box must be set before this step
    const double viscosity = runConfig.viscosity;
    const auto &rodContainer = *containerPtr;

    const int nPtsLocal = rodPtsIndex.back(); // number of quadrature points on local
    const int nSrc = nPtsLocal;
    const int nTrg = nPtsLocal;
    srcSLValue.clear();
    srcSLValue.resize(4 * nSrc, 0); // RPY FMM, fx,fy,fz,a
    trgValue.clear();
    trgValue.resize(6 * nTrg, 0); // RPY FMM, vx,vy,vz,lapvx,lapvy,lapvz

#pragma omp parallel for
    for (size_t i = 0; i < nRodLocal; i++) {
        const auto &sy = rodContainer[i];
        const double length = sy.length;
        const int numQuadPt = sy.numQuadPt;
        // const auto quadPtr = sy.quadPtr;
        assert(sy.numQuadPt == rodPts[i]);

        const Evec3 direction = ECmapq(sy.orientation) * Evec3(0, 0, 1);
        const int idx = rodPtsIndex[i];

        for (size_t j = 0; j < numQuadPt; j++) {
            Evec3 loc = ECmap3(sy.pos) + (length * 0.5 * sloc[idx + j]) * direction;
            const double w = this->weight[idx + j];
            // NOTE: This is the line integral over the rods from -1 to 1,
            // thus we multiply the force by the respective weight scaled by the half length.
            srcSLValue[4 * (idx + j) + 0] = fhPtr(3 * (idx + j) + 0, 0) * w * length * 0.5;
            srcSLValue[4 * (idx + j) + 1] = fhPtr(3 * (idx + j) + 1, 0) * w * length * 0.5;
            srcSLValue[4 * (idx + j) + 2] = fhPtr(3 * (idx + j) + 2, 0) * w * length * 0.5;
            srcSLValue[4 * (idx + j) + 3] = radius[idx + j];
        }
    }

    // run FMM
    fmmPtr->clearFMM(KERNEL::RPY);
    fmmPtr->evaluateFMM(KERNEL::RPY, nSrc, srcSLValue.data(), nTrg, trgValue.data());
    commRcp->barrier();

    uinf.clear();
    uinf.resize(3 * nPtsLocal, 0);

    // subtract off self interactions and compute uinf
#pragma omp parallel for
    for (int i = 0; i < nRodLocal; i++) {
        const auto &sy = rodContainer[i];
        const int numQuadPt = sy.numQuadPt;
        assert(sy.numQuadPt == rodPts[i]);
        const auto quadPtr = sy.quadPtr;
        assert(sy.numQuadPt == quadPtr->getSize());

        const int idx = rodPtsIndex[i];
        std::vector<double> trgValueSelf(6 * numQuadPt, 0);
        fmmPtr->evaluateKernel(KERNEL::RPY, 1, PPKERNEL::SLS2T, numQuadPt, srcSLCoord.data() + 3 * idx,
                               srcSLValue.data() + 4 * idx, numQuadPt, trgCoord.data() + 3 * idx, trgValueSelf.data());
        for (int j = 0; j < 6 * numQuadPt; j++) {
            trgValue[6 * idx + j] -= trgValueSelf[j];
        }

        for (int j = 0; j < numQuadPt; j++) {
            const double radiusPt = radius[idx + j];
            const double rpyfac = radiusPt * radiusPt * (1.0 / 6.0);
            uinf[3 * (idx + j) + 0] =
                (trgValue[6 * (idx + j) + 0] + rpyfac * trgValue[6 * (idx + j) + 3]) * (1 / viscosity);
            uinf[3 * (idx + j) + 1] =
                (trgValue[6 * (idx + j) + 1] + rpyfac * trgValue[6 * (idx + j) + 4]) * (1 / viscosity);
            uinf[3 * (idx + j) + 2] =
                (trgValue[6 * (idx + j) + 2] + rpyfac * trgValue[6 * (idx + j) + 5]) * (1 / viscosity);
        }
    }

    if (commRcp->getRank() == 0) {
        printf("FMM complete\n");
    }
}

template <class Container>
void HydroRodOperator<Container>::calcuinfOS() const {
    // Constants
    constexpr double Pi = M_PI;
    const double mu = runConfig.viscosity;
    const auto &rodContainer = *containerPtr;
    const int nPtsLocal = rodPtsIndex.back();
    Oinf.clear();
    Oinf.resize(3 * nRodLocal, 0);
    Sinf.clear();
    Sinf.resize(3 * nRodLocal, 0);
    uinfOS.clear();
    uinfOS.resize(nPtsLocal * 3, 0);

#pragma omp parallel for
    for (int i = 0; i < nRodLocal; i++) {
        const auto &sy = rodContainer[i];
        const int numQuadPt = sy.numQuadPt;
        assert(sy.numQuadPt == rodPts[i]);
        auto quadPtr = sy.quadPtr;
        assert(quadPtr->getSize() == numQuadPt);

        // Rod Properties
        const double length = sy.length;
        const double ell = length * 0.5;
        const double diameter = sy.radius * 2;
        const double c = std::log(2 * length / diameter) / (4 * Pi * mu);

        const Evec3 p = ECmapq(sy.orientation) * Evec3(0, 0, 1);
        const Emat3 pp = p * p.transpose();
        const Emat3 Impp = Emat3::Identity() - pp;
        const Emat3 Imhalfpp = Emat3::Identity() - 0.5 * pp;
        const int idx = rodPtsIndex[i];

        // compute Oinf and Ssinf
        // calc ux, uy, uz from FMM trgValue
        std::vector<double> ux(numQuadPt, 0);
        std::vector<double> uy(numQuadPt, 0);
        std::vector<double> uz(numQuadPt, 0);
        for (int j = 0; j < numQuadPt; j++) {
            ux[j] = uinf[3 * (idx + j) + 0];
            uy[j] = uinf[3 * (idx + j) + 1];
            uz[j] = uinf[3 * (idx + j) + 2];
        }

        // calc Hx,Hy,Hz
        const double intScaleFac = ell;
        Evec3 OinfRod, SinfRod;
        OinfRod[0] = intScaleFac * quadPtr->intSamples(ux.data());
        OinfRod[1] = intScaleFac * quadPtr->intSamples(uy.data());
        OinfRod[2] = intScaleFac * quadPtr->intSamples(uz.data());
        SinfRod[0] = intScaleFac * intScaleFac * quadPtr->intSSamples(ux.data());
        SinfRod[1] = intScaleFac * intScaleFac * quadPtr->intSSamples(uy.data());
        SinfRod[2] = intScaleFac * intScaleFac * quadPtr->intSSamples(uz.data());

        for (int k = 0; k < 3; k++) {
            Oinf[3 * i + k] = OinfRod[k];
            Sinf[3 * i + k] = SinfRod[k];
        }

        const Evec3 Hpart1 = -(1 / (2 * c * ell)) * (Imhalfpp * OinfRod);
        const Evec3 Hpart2 = -(3 / (2 * c * ell * ell * ell)) * (Impp * SinfRod);

        for (int j = 0; j < numQuadPt; j++) {
            const int idj = 3 * (idx + j);
            const Evec3 uinfPt = ECmap3(uinf.data() + idj);
            const double s = sloc[idx + j] * length * 0.5;

            Evec3 uinfOSPt = (1 / c) * Imhalfpp * uinfPt + Hpart1 + s * Hpart2;
            uinfOS[idj + 0] = uinfOSPt[0];
            uinfOS[idj + 1] = uinfOSPt[1];
            uinfOS[idj + 2] = uinfOSPt[2];
        }
    }
}

template <class Container>
void HydroRodOperator<Container>::calcb(const TV &FTinputVec, TV &bVec, const HydroOption &options) const {

    assert(bVec.getMap()->isSameAs(*pointValuesMapRcp));
    bVec.putScalar(0);
    auto bVecPtr = bVec.getLocalView<Kokkos::HostSpace>(Tpetra::Access::ReadWrite);

    assert(FTinputVec.getLocalLength() == 6 * nRodLocal);
    auto FTinputPtr = FTinputVec.getLocalView<Kokkos::HostSpace>(Tpetra::Access::ReadOnly);

    // Constants
    constexpr double Pi = M_PI;
    const double mu = runConfig.viscosity;
    const auto &rodContainer = *containerPtr;
    const int nPtsLocal = rodPtsIndex.back();

    // step 1, freestreamvelocity should not appear
    // step 2, swim velocity
    if (options.withSwim) {
#pragma omp parallel for
        for (size_t i = 0; i < nRodLocal; i++) {
            const auto &sy = rodContainer[i];
            const int numQuadPt = sy.numQuadPt;
            assert(sy.numQuadPt == rodPts[i]);
            // Rod Properties
            const double length = sy.length;
            const double ell = length * 0.5;
            const double diameter = sy.radius * 2;
            const double c = std::log(2 * length / diameter) / (4 * Pi * mu);
            const Evec3 p = ECmapq(sy.orientation) * Evec3(0, 0, 1);
            const int idx = rodPtsIndex[i];
            for (size_t j = 0; j < numQuadPt; j++) {
                // Hu, Hf term
                const double pfactor = (Hf[i] / (2 * ell) - Hu[i] / (4 * c * ell) + slipVelocity[idx + j] / (2 * c) -
                                        activeStress[idx + j]);
                bVecPtr(3 * (idx + j) + 0, 0) += p[0] * pfactor;
                bVecPtr(3 * (idx + j) + 1, 0) += p[1] * pfactor;
                bVecPtr(3 * (idx + j) + 2, 0) += p[2] * pfactor;
            }
        }
        // calc uinf generated by fh
        runFMM(*fsRcp);
        // testMixedQuadrature(*fsRcp);

        // calc operator
        calcuinfOS();
#pragma omp parallel for
        for (size_t i = 0; i < nPtsLocal; i++) {
            bVecPtr(3 * i + 0, 0) -= uinfOS[3 * i + 0];
            bVecPtr(3 * i + 1, 0) -= uinfOS[3 * i + 1];
            bVecPtr(3 * i + 2, 0) -= uinfOS[3 * i + 2];
        }
    }

    // step 3, FTExt and FTinput
#pragma omp parallel for
    for (size_t i = 0; i < nRodLocal; i++) {
        const auto &sy = rodContainer[i];
        const int numQuadPt = sy.numQuadPt;
        assert(sy.numQuadPt == rodPts[i]);

        // Rod Properties
        const double length = sy.length;
        const double ell = length * 0.5;
        const double diameter = sy.radius * 2;
        const double c = std::log(2 * length / diameter) / (4 * Pi * mu);

        const Evec3 p = ECmapq(sy.orientation) * Evec3(0, 0, 1);
        const Emat3 pp = p * p.transpose();
        const Emat3 Impp = Emat3::Identity() - pp;
        const Emat3 Imhalfpp = Emat3::Identity() - 0.5 * pp;

        const Evec3 FExt = ECmap3(FTExt.data() + 6 * i);
        const Evec3 TExt = ECmap3(FTExt.data() + 6 * i + 3);
        const Evec3 Finput(FTinputPtr(6 * i + 0, 0), FTinputPtr(6 * i + 1, 0), FTinputPtr(6 * i + 2, 0));
        const Evec3 Tinput(FTinputPtr(6 * i + 3, 0), FTinputPtr(6 * i + 4, 0), FTinputPtr(6 * i + 5, 0));
        const Evec3 FTotal = FExt + Finput;
        const Evec3 TTotal = TExt + Tinput;
        const Evec3 Tcrossp = TTotal.cross(p);

        const int idx = rodPtsIndex[i];
        for (size_t j = 0; j < numQuadPt; j++) {
            const int idj = 3 * (idx + j);
            const double s = sloc[idx + j] * length * 0.5;
            for (size_t k = 0; k < 3; k++) {
                bVecPtr(idj + k, 0) += (1 / (2 * ell)) * FTotal[k] + (3 * s / (2 * ell * ell * ell)) * Tcrossp[k];
            }
        }
    }
}

template <class Container>
void HydroRodOperator<Container>::calcVelOmega(const TV &feVec, TV &VelOmegaVec, const HydroOption &options) const {
    Teuchos::RCP<TV> fhRcp = Teuchos::rcp(new TV(pointValuesMapRcp, true));
    fhRcp->update(1.0, feVec, 0.0);
    if (options.withSwim) {
        fhRcp->update(1.0, *fsRcp, 1.0); // fh = fe + fs(s)p
    }

    runFMM(*fhRcp);
    // testMixedQuadrature(*fhRcp);

    calcuinfOS();

    auto VelOmegaPtr = VelOmegaVec.getLocalView<Kokkos::HostSpace>(Tpetra::Access::ReadWrite);

    // Constants
    const double mu = runConfig.viscosity;
    const auto &rodContainer = *containerPtr;
    auto fePtr = feVec.getLocalView<Kokkos::HostSpace>(Tpetra::Access::ReadOnly);

#pragma omp parallel for
    for (size_t i = 0; i < nRodLocal; i++) {
        const auto &sy = rodContainer[i];
        const int numQuadPt = sy.numQuadPt;
        auto quadPtr = sy.quadPtr;
        assert(sy.numQuadPt == rodPts[i]);
        assert(sy.numQuadPt == quadPtr->getSize());

        // Rod Properties
        const double length = sy.length;
        const double ell = length * 0.5;
        const double diameter = sy.radius * 2;
        const double c = std::log(2 * length / diameter) / (4 * Pi * mu);

        const Evec3 p = ECmapq(sy.orientation) * Evec3(0, 0, 1);
        const Emat3 pp = p * p.transpose();
        const Emat3 Impp = Emat3::Identity() - pp;
        const Emat3 Imhalfpp = Emat3::Identity() - 0.5 * pp;
        const double intFactor = length * 0.5;
        const int idx = rodPtsIndex[i];

        // these are force/torque external specified in HydroTable
        const Evec3 FExt = ECmap3(FTExt.data() + 6 * i);
        const Evec3 TExt = ECmap3(FTExt.data() + 6 * i + 3);

        // these are computed from imposed fe along quadrature points
        // no need to add fs because integral of fs is zero and torque is always zero
        Evec3 Finput = Evec3::Zero();
        Evec3 Tinput = Evec3::Zero();
        Evec3 Finput_ints = Evec3::Zero();
        for (int k = 0; k < 3; k++) {
            std::vector<double> f(numQuadPt, 0); // component of force
            for (int j = 0; j < numQuadPt; j++) {
                f[j] = fePtr(3 * (idx + j) + k, 0);
            }
            Finput[k] = quadPtr->intSamples(f.data()) * intFactor;
            Finput_ints[k] = quadPtr->intSSamples(f.data()) * (intFactor * intFactor);
        }
        // std::cout << Finput.transpose();
        Tinput = p.cross(Finput_ints);

        const Evec3 FTotal = (options.withForceTorqueExternal ? FExt : Evec3::Zero()) + Finput;
        const Evec3 TTotal = (options.withForceTorqueExternal ? TExt : Evec3::Zero()) + Tinput;

        const Evec3 OinfRod = ECmap3(Oinf.data() + 3 * i);
        const Evec3 SinfRod = ECmap3(Sinf.data() + 3 * i);
        const double HuRod = options.withSwim ? Hu[i] : 0;
        const double HfRod = options.withSwim ? Hf[i] : 0;
        const Evec3 FTotalPara = (FTotal.dot(p)) * p;
        const Evec3 FTotalPerp = FTotal - FTotalPara;

        const Evec3 Omega = 3 / (2 * ell * ell * ell) * (p.cross(SinfRod) + c * TTotal);
        const Evec3 Vel = OinfRod / (2 * ell) + p * (HfRod * c / ell - HuRod / (2 * ell)) + c * FTotalPerp / (2 * ell) +
                          c * FTotalPara / ell;
        // std::cout << i << " " << Vel.transpose() << " " << Omega.transpose();
        VelOmegaPtr(6 * i + 0, 0) = Vel[0];
        VelOmegaPtr(6 * i + 1, 0) = Vel[1];
        VelOmegaPtr(6 * i + 2, 0) = Vel[2];
        VelOmegaPtr(6 * i + 3, 0) = Omega[0];
        VelOmegaPtr(6 * i + 4, 0) = Omega[1];
        VelOmegaPtr(6 * i + 5, 0) = Omega[2];
    }
}

template <class Container>
void HydroRodOperator<Container>::calcAfe(const TV &feVec, TV &AfeVec) const {
    // 1 FMM call to get u and ulap induced by f. Output is stored in trgValue
    runFMM(feVec);
    // testMixedQuadrature(feVec);
    calcuinfOS();

    const auto &feVecPtr = feVec.getLocalView<Kokkos::HostSpace>(Tpetra::Access::ReadOnly);
    auto AfeVecPtr = AfeVec.getLocalView<Kokkos::HostSpace>(Tpetra::Access::ReadWrite);

    // Constants
    const double Pi = M_PI;
    const double mu = runConfig.viscosity;
    const auto &rodContainer = *containerPtr;

#pragma omp parallel for
    for (size_t i = 0; i < nRodLocal; i++) {
        const auto &sy = rodContainer[i];
        const int numQuadPt = sy.numQuadPt;
        assert(sy.numQuadPt == rodPts[i]);
        const auto quadPtr = sy.quadPtr;
        assert(sy.numQuadPt == quadPtr->getSize());

        // Rod Properties
        const double length = sy.length;
        const double ell = length * 0.5;
        const double diameter = sy.radius * 2;
        const double c = std::log(2 * length / diameter) / (4 * Pi * mu);

        const Evec3 direction = ECmapq(sy.orientation) * Evec3(0, 0, 1);
        const Evec3 q = ECmapq(sy.orientation) * Evec3(0, 0, 1);
        const Emat3 qq = q * q.transpose();
        const Emat3 Imqq = Emat3::Identity() - qq;
        const Emat3 Imhalfqq = Emat3::Identity() - 0.5 * qq;

        const int idx = rodPtsIndex[i];

        Evec3 Hpart1 = -(1 / (2 * c * ell)) * (Imhalfqq * ECmap3(Oinf.data() + 3 * i));
        Evec3 Hpart2 = -(3 / (2 * c * ell * ell * ell)) * (Imqq * ECmap3(Sinf.data() + 3 * i));
        // fill Afe
        for (size_t j = 0; j < numQuadPt; j++) {
            const double s = sloc[idx + j] * length * 0.5;
            const Evec3 f(feVecPtr(3 * (idx + j) + 0, 0), feVecPtr(3 * (idx + j) + 1, 0),
                          feVecPtr(3 * (idx + j) + 2, 0));
            const Evec3 uinfPt = ECmap3(uinf.data() + 3 * (idx + j));
            Evec3 Af = f + (1.0 / c) * (Imhalfqq * uinfPt) + Hpart1 + s * Hpart2;
            AfeVecPtr(3 * (idx + j) + 0, 0) = Af[0];
            AfeVecPtr(3 * (idx + j) + 1, 0) = Af[1];
            AfeVecPtr(3 * (idx + j) + 2, 0) = Af[2];
        }
    }
}

template <class Container>
void HydroRodOperator<Container>::apply(const TMV &X, TMV &Y, Teuchos::ETransp mode, scalar_type alpha,
                                        scalar_type beta) const {
    // compute Y=alpha*Ax+beta*Y;
    assert(X.getMap()->isSameAs(*(Y.getMap())));
    assert(X.getMap()->isSameAs(*pointValuesMapRcp));
    assert(X.getNumVectors() == Y.getNumVectors());

    Teuchos::RCP<TV> YColOld = Teuchos::rcp(new TV(Y.getMap(), true));

    const int nCol = X.getNumVectors();
    for (int c = 0; c < nCol; c++) {
        const auto &XCol = X.getVector(c);
        auto YCol = Y.getVectorNonConst(c);
        YColOld->update(beta, *YCol, Teuchos::ScalarTraits<scalar_type>::zero()); // Yold = beta*Ycol
        calcAfe(*XCol, *YCol);                                                    // Ycol = AXcol
        YCol->update(Teuchos::ScalarTraits<scalar_type>::one(), *YColOld, alpha); // Ycol = alpha*AXcol+beta*Ycol
    }
}

template <class Container>
void HydroRodOperator<Container>::testMixedQuadrature(const TV &fhVec) const {
    // task: compute uinf with given force density fhVec.
    using namespace stkfmm;
    assert(fhVec.getLocalLength() == 3 * this->rodPtsIndex.back());
    assert(fhVec.getLocalLength() == this->pointValuesMapRcp->getLocalNumElements());
    auto fhPtr = fhVec.getLocalView<Kokkos::HostSpace>(Tpetra::Access::ReadOnly);

    // the FMM bounding box must be set before this step
    const double viscosity = runConfig.viscosity;
    const auto &rodContainer = *containerPtr;

    std::vector<double> srcValueTest;
    const int nPtsLocal = rodPtsIndex.back(); // number of quadrature points on local
    const int nSrc = nPtsLocal;
    const int nTrg = nPtsLocal;
    srcSLValue.clear();
    srcSLValue.resize(4 * nSrc, 0); // RPY FMM, fx,fy,fz,a
    srcValueTest.clear();
    srcValueTest.resize(4 * nSrc, 0); // RPY FMM, fx,fy,fz,a
    trgValue.clear();
    trgValue.resize(6 * nTrg, 0); // RPY FMM, vx,vy,vz,lapvx,lapvy,lapvz

#pragma omp parallel for
    for (int i = 0; i < nRodLocal; i++) {
        const auto &sy = rodContainer[i];
        const double length = sy.length;
        const int numQuadPt = sy.numQuadPt;
        // const auto quadPtr = sy.quadPtr;
        assert(sy.numQuadPt == rodPts[i]);

        const Evec3 direction = ECmapq(sy.orientation) * Evec3(0, 0, 1);
        const int idx = rodPtsIndex[i];

        for (int j = 0; j < numQuadPt; j++) {
            Evec3 loc = ECmap3(sy.pos) + (length * 0.5 * sloc[idx + j]) * direction;
            const double w = this->weight[idx + j];
            // NOTE: This is the line integral over the rods from -1 to 1,
            // thus we multiply the force by the respective weight scaled by the half length.
            srcSLValue[4 * (idx + j) + 0] = fhPtr(3 * (idx + j) + 0, 0) * w * length * 0.5;
            srcSLValue[4 * (idx + j) + 1] = fhPtr(3 * (idx + j) + 1, 0) * w * length * 0.5;
            srcSLValue[4 * (idx + j) + 2] = fhPtr(3 * (idx + j) + 2, 0) * w * length * 0.5;
            srcSLValue[4 * (idx + j) + 3] = radius[idx + j];
            srcValueTest[4 * (idx + j) + 0] = fhPtr(3 * (idx + j) + 0, 0);
            srcValueTest[4 * (idx + j) + 1] = fhPtr(3 * (idx + j) + 1, 0);
            srcValueTest[4 * (idx + j) + 2] = fhPtr(3 * (idx + j) + 2, 0);
            srcValueTest[4 * (idx + j) + 3] = radius[idx + j];
        }
    }

    fmmPtr->clearFMM(KERNEL::RPY);
    commRcp->barrier();
    uinf.clear();
    uinf.resize(3 * nPtsLocal, 0);

// compute uinf via direct O(N^2) summation (non-periodic)
#pragma omp parallel for
    for (int i = 0; i < nRodLocal; i++) {
        const auto &sy = rodContainer[i];
        const int numQuadPt = sy.numQuadPt;
        assert(sy.numQuadPt == rodPts[i]);

        const int idx = rodPtsIndex[i];
        const double lineHalfLength = sy.length / 2;
        const Evec3 direction = ECmapq(sy.orientation) * Evec3(0, 0, 1);
        const Evec3 pos = ECmap3(sy.pos);
        for (int t = 0; t < nPtsLocal; t++) {

            // skip the current rod
            bool skip_tf = false;
            for (int n = 0; n < numQuadPt; n++) {
                if ((srcSLCoord[3 * idx + 3 * n + 0] == trgCoord[3 * t + 0]) &&
                    (srcSLCoord[3 * idx + 3 * n + 1] == trgCoord[3 * t + 1]) &&
                    (srcSLCoord[3 * idx + 3 * n + 2] == trgCoord[3 * t + 2])) {
                    skip_tf = true;
                }
            }
            if (skip_tf) {
                continue;
            }

            // check if special quad is necessary.
            const Evec3 tmp(trgCoord.data() + 3 * t);
            double R = (pos - tmp).norm();
            SpecialQuadWeights<32> sqw(numQuadPt); // The recommended maximum number of quadrature points is 32
            if (R > 0.75) {
                const double *w1 = sqw.getGLWeights();
                const double *w3 = sqw.getGLWeights();
                const double *w5 = sqw.getGLWeights();
                rpy_weighted_ulapu(numQuadPt, srcSLCoord.data() + 3 * idx, srcValueTest.data() + 4 * idx,
                                   trgCoord.data() + 3 * t, trgValue.data() + t * 6, w1, w3, w5, lineHalfLength);
            } else {
                sqw.calcWeights(lineHalfLength, pos.data(), trgCoord.data() + 3 * t, direction.data());
                const double *w1 = sqw.getWeights1();
                const double *w3 = sqw.getWeights3();
                const double *w5 = sqw.getWeights5();

                // Uncomment to check against SpecialQuadWeights_test.cpp
                // std::cout << "pos " << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
                // std::cout << "trgCoord " << trgCoord[3 * t + 0] << " " << trgCoord[3 * t + 1] << " " << trgCoord[3 *
                // t + 2]  << std::endl; std::cout << "direction " << direction[0] << " " << direction[1] << " " <<
                // direction[2] << std::endl;
                // for (int i=0; i<numQuadPt; i++) {
                //     std::cout << "w1[i] " << w1[i] << std::endl;
                // }
                // for (int i=0; i<numQuadPt; i++) {
                //     std::cout << "w3[i] " << w3[i] << std::endl;
                // }
                // for (int i=0; i<numQuadPt; i++) {
                //     std::cout << "w5[i] " << w5[i] << std::endl;
                // }
                rpy_weighted_ulapu(numQuadPt, srcSLCoord.data() + 3 * idx, srcValueTest.data() + 4 * idx,
                                   trgCoord.data() + 3 * t, trgValue.data() + t * 6, w1, w3, w5, lineHalfLength);
            }
        }

        for (int j = 0; j < numQuadPt; j++) {
            const double radiusPt = radius[idx + j];
            const double rpyfac = radiusPt * radiusPt * (1.0 / 6.0);
            uinf[3 * (idx + j) + 0] =
                (trgValue[6 * (idx + j) + 0] + rpyfac * trgValue[6 * (idx + j) + 3]) * (1 / viscosity);
            uinf[3 * (idx + j) + 1] =
                (trgValue[6 * (idx + j) + 1] + rpyfac * trgValue[6 * (idx + j) + 4]) * (1 / viscosity);
            uinf[3 * (idx + j) + 2] =
                (trgValue[6 * (idx + j) + 2] + rpyfac * trgValue[6 * (idx + j) + 5]) * (1 / viscosity);
        }
    }

    if (commRcp->getRank() == 0) {
        printf("testMixedQuadrature complete\n");
    }
}

template <class Container>
void HydroRodOperator<Container>::writeFlowGridVTI(const std::string &prefix, const double &flowGrid, const bool &wallz0, const double& viscosity, const int &fileID) {
    // Evaluate the hydrodynamic interaction between the rods and a structured grid.
    // If wallz0 is true, the grid will be in the upper half space.

    // write a VTK rectilinear grid file
    if (!fmmPtr || flowGrid <= 0) {
        return;
    }

    // pre: valid sphere list and fdist
    // post: pointFlowList hold the flow in the structured grid with predifined
    // max size,
    //      flow dumpped to a file
    // generate the fluid velocity with 3D Cartesian mesh size dx, and write it
    // to filename using Lscale, Tscale as the dimension scale of mesh size.
    // adding a shift if necessary.

    // pay attention that the FMM box size is usually larger than the sphere
    // periodic boundary size except for the TP:PXYZ case. Therefore the shift
    // should be carefully set to match the data

    // step 1 generate cartesian mesh on rank0
    // total mesh size

    const auto box = fmmPtr->getBox();
    using std::get;
    double xlow = get<0>(box);
    double xhigh = get<1>(box);
    double ylow = get<2>(box);
    double yhigh = get<3>(box);
    double zlow = get<4>(box);
    double zhigh = get<5>(box);

    if (wallz0) {
        zhigh = 0.5 * (zlow + zhigh);
    }

    const double dxSet = flowGrid;
    const int NX = (xhigh - xlow) / dxSet;
    const int NY = (yhigh - ylow) / dxSet;
    const int NZ = (zhigh - zlow) / dxSet;
    const double dx = (xhigh - xlow) / (NX - 1);
    const double dy = (yhigh - ylow) / (NY - 1);
    const double dz = (zhigh - zlow) / (NZ - 1);

    std::cout << "Box: [" << xlow << ", " << ylow << ", " << zlow << "] to [" << xhigh << ", " << yhigh << ", " << zhigh
              << "]" << std::endl;
    std::cout << "NX: " << NX << " NY: " << NY << " NZ: " << NZ << std::endl;
    std::cout << "dx: " << dx << " dy: " << dy << " dz: " << dz << std::endl;

    if (commRcp->getRank() != 0) {
        trgCoord.clear();
        trgValue.clear();
    } else {
        std::cout << "Dumping flow grid VTI" << std::endl;
        /* TO match VTK point ordering:
            *  # NOTE: VTK expects data in FORTRAN order
            *  The order and number of points must match that specified by the
            * dimensions of the grid. The point order increases in i fastest (from
            * 0<=i<dims[0]), then j (0<=j<dims[1]), then k (0<=k<dims[2]) where
            * dims[] are the dimensions of the grid in the i-j-k topological
            * directions. The number of points is dims[0]*dims[1]*dims[2].
            *
            *  The same is true for the cells of the grid.
            *  The order and number of cells must match that specified by the
            * dimensions of the grid. The cell order increases in i fastest (from
            * 0<=i<(dims[0]-1)), then j (0<=j<(dims[1]-1)), then k
            * (0<=k<(dims[2]-1)) The number of cells is
            * (dims[0]-1)*(dims[1]-1)*(dims[2]-1).
            * */
        trgCoord.clear();
        trgCoord.resize(NX * NY * NZ * 3);
        trgValue.clear();
        trgValue.resize(NX * NY * NZ * 6, 0);
#pragma omp parallel for
        for (int k = 0; k < NZ; k++) {
            for (int j = 0; j < NY; j++) {
                for (int i = 0; i < NX; i++) {
                    trgCoord[(i + j * (NX) + k * (NX) * (NY)) * 3 + 0] = i * dx + xlow;
                    trgCoord[(i + j * (NX) + k * (NX) * (NY)) * 3 + 1] = j * dy + ylow;
                    trgCoord[(i + j * (NX) + k * (NX) * (NY)) * 3 + 2] = k * dz + zlow;
                }
            }
        }

        // shift the lower bound a little
        {
        int i = 0;
    #pragma omp parallel for
            for (int k = 0; k < NZ; k++) {
                for (int j = 0; j < NY; j++) {
                    // shift x
                    trgCoord[(i + j * (NX) + k * (NX) * (NY)) * 3 + 0] += 100 *
                    std::numeric_limits<float>::epsilon();
                }
            }
        }

        {   
        int j = 0;
    #pragma omp parallel for
            for (int k = 0; k < NZ; k++) {
                for (int i = 0; i < NX; i++) {
                    // shift y
                    trgCoord[(i + j * (NX) + k * (NX) * (NY)) * 3 + 1] += 100 *
                    std::numeric_limits<float>::epsilon();
                }
            }
        }

        {
        int k = 0;
    #pragma omp parallel for
            for (int j = 0; j < NY; j++) {
                for (int i = 0; i < NX; i++) {
                    // shift z
                    trgCoord[(i + j * (NX) + k * (NX) * (NY)) * 3 + 2] += 100 *
                    std::numeric_limits<float>::epsilon();
                }
            }
        }


        // shift the upper bound a little
        {
        int i = NX - 1;
    #pragma omp parallel for
            for (int k = 0; k < NZ; k++) {
                for (int j = 0; j < NY; j++) {
                    // shift x
                    trgCoord[(i + j * (NX) + k * (NX) * (NY)) * 3 + 0] -= 100 *
                    std::numeric_limits<float>::epsilon();
                }
            }
        }

        {   
        int j = NY - 1;
    #pragma omp parallel for
            for (int k = 0; k < NZ; k++) {
                for (int i = 0; i < NX; i++) {
                    // shift y
                    trgCoord[(i + j * (NX) + k * (NX) * (NY)) * 3 + 1] -= 100 *
                    std::numeric_limits<float>::epsilon();
                }
            }
        }

        {
        int k = NZ - 1;
    #pragma omp parallel for
            for (int j = 0; j < NY; j++) {
                for (int i = 0; i < NX; i++) {
                    // shift z
                    trgCoord[(i + j * (NX) + k * (NX) * (NY)) * 3 + 2] -= 100 *
                    std::numeric_limits<float>::epsilon();
                }
            }
        }
    }

    commRcp->barrier();

    // step 2 set up FMM tree and calc velocity, on every node
    // srcSLCoord and srcSLValue have been set. No need to change

    const double invMu = 1.0 / viscosity;
    const auto kernel = stkfmm::KERNEL::RPY;
    const int nSL = srcSLCoord.size() / 3;
    const int nTrg = trgCoord.size() / 3;
    fmmPtr->clearFMM(kernel);
    fmmPtr->setPoints(nSL, srcSLCoord.data(), nTrg, trgCoord.data());
    fmmPtr->setupTree(kernel);
    fmmPtr->evaluateFMM(kernel, nSL, srcSLValue.data(), nTrg, trgValue.data());
    commRcp->barrier();

    if (commRcp->getRank() == 0) {
        // const int pointFlowNumber = trgCoord.size() / 3;
        // final step: dump
        std::string filename = prefix + std::string("Flow_") + std::to_string(fileID) + ".vti";

        // apply viscosity
        std::vector<double> vel(3 * nTrg), lap(3 * nTrg);
#pragma omp parallel for
        for (int i = 0; i < nTrg; i++) {
            vel[3 * i] = trgValue[6 * i] * invMu;
            vel[3 * i + 1] = trgValue[6 * i + 1] * invMu;
            vel[3 * i + 2] = trgValue[6 * i + 2] * invMu;
            lap[3 * i] = trgValue[6 * i + 3] * invMu;
            lap[3 * i + 1] = trgValue[6 * i + 4] * invMu;
            lap[3 * i + 2] = trgValue[6 * i + 5] * invMu;
        }

        std::ofstream file(filename, std::ios::out);
        // VTR file header
        int extentLow[3] = {0, 0, 0};
        int extentHigh[3] = {NX - 1, NY - 1, NZ - 1};
        double spacing[3] = {dx, dy, dz};
        double boxLow[3] = {xlow, ylow, zlow};
        IOHelper::writeHeadVTI(file, extentLow, extentHigh, spacing, boxLow);
        file << "<Piece Extent=\"" << extentLow[0] << " " << extentHigh[0] << " " << extentLow[1] << " "
                << extentHigh[1] << " " << extentLow[2] << " " << extentHigh[2] << "\">\n";

        // VTI point data
        file << "<PointData Scalars=\"scalars\">\n";
        IOHelper::writeDataArrayBase64(vel, "vel", 3, file);
        IOHelper::writeDataArrayBase64(lap, "laplacian", 3, file);
        file << "</PointData>\n";

        file << "</Piece>\n";
        IOHelper::writeTailVTI(file);
        file.close();
    }

    commRcp->barrier();
}

// template <class Container>
// Teuchos::RCP<TV> HydroRodOperator<Container>::getfhRcp() const {
//     Teuchos::RCP<TV> fhRcp = Teuchos::rcp(new TV(pointValuesMapRcp, true));
//     auto fhPtr = fhRcp->getLocalView<Kokkos::HostSpace>(Tpetra::Access::ReadWrite);
//     const auto &rodContainer = *containerPtr;
// #pragma omp parallel for
//     for (int i = 0; i < nRodLocal; i++) {
//         const auto &sy = rodContainer[i];
//         const int idx = rodPtsIndex[i];
//         const int numQuadPts = rodPts[i];
//         for (int j = 0; j < numQuadPts; j++) {
//             fhPtr(3 * (idx + j) + 0, 0) = sy.forceHydro[3 * j];
//             fhPtr(3 * (idx + j) + 1, 0) = sy.forceHydro[3 * j + 1];
//             fhPtr(3 * (idx + j) + 2, 0) = sy.forceHydro[3 * j + 2];
//         }
//     }

//     return fhRcp;
// };

#endif