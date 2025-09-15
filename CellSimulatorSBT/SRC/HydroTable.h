#ifndef HYDROTABLE_H_
#define HYDROTABLE_H_

#include <string>
#include <vector>

/**
 * @brief Hydro quadrature point data table
 *
 */
class HydroTable {
  public:
    std::vector<double> sloc;         ///< location [-1,1]
    std::vector<double> slipVelocity; ///< [us0,us1,us2,...]
    std::vector<double> activeStress; ///< [fs1,fs2,fs3,...]
    std::vector<double> radius;       ///< [r0,r1,r2,...]

    /**
     * @brief Construct a new HydroTable object
     *
     */
    HydroTable() = default;
    // HydroTable(const HydroTable &) = default;
    // HydroTable &operator=(const HydroTable &) = default;
    // HydroTable(HydroTable &&) = default;
    // HydroTable &operator=(HydroTable &&) = default;

    /**
     * @brief Construct a new HydroTable object
     *
     * @param csvFile read a csv HydroTable file
     */
    HydroTable(const std::string &csvFile);

    void getValue(const int nPts, const double *sPtr, double *slipVelocityPtr, double *activeStressPtr,
                  double *radiusPtr) const;

    void echo() const;
};

#endif
