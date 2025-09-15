#ifndef GROWTHTABLE_H_
#define GROWTHTABLE_H_

#include <string>
#include <vector>

/**
 * @brief 2D interpolation table f(x,y)
 *
 */
class GrowthTable {
    std::vector<double> xd; // [x0,x1,x2,...]
    std::vector<double> yd; // [y0,y1,y2,...]
    std::vector<double> fd; // [f00,f01,f02,..f10,f11,f12,...]
    std::string xname;
    std::string yname;

    /**
     * @brief Construct a new Table object
     *
     */
  public:
    GrowthTable() = default;
    // GrowthTable(const Table &) = default;
    // GrowthTable &operator=(const GrowthTable &) = default;
    // GrowthTable(GrowthTable &&) = default;
    // GrowthTable &operator=(GrowthTable &&) = default;

    /**
     * @brief Construct a new GrowthTable object
     *
     * @param csvFile read a csv GrowthTable file
     */
    GrowthTable(const std::string &csvFile, const std::string &xname, const std::string &yname);

    std::vector<double> getValue(const std::vector<double> &xi, const std::vector<double> &yi) const;

    void echo() const;
};

#endif
