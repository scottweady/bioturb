#include "HydroTable.h"

#include "mlinterp.hpp"
#include "rapidcsv.h"

#include <cstdio>
#include <iostream>
#include <set>
#include <string>

HydroTable::HydroTable(const std::string &csvFile) {
    // expects a csv file
    // sloc between [-1,1]
    // sloc,slipVelocity,activeStress,radius
    // s0,us0,fs0,r1
    // s1,us1,fs1,r2
    // s2,us2,fs2,r3
    // ...

    rapidcsv::Document doc(csvFile, rapidcsv::LabelParams(0, -1));
    std::vector<std::string> columnNames = doc.GetColumnNames();
    for (auto &c : columnNames) {
        std::cout << c << std::endl;
    }
    sloc = doc.GetColumn<double>("sloc");
    slipVelocity = doc.GetColumn<double>("slipVelocity");
    activeStress = doc.GetColumn<double>("activeStress");
    radius = doc.GetColumn<double>("radius");
}

void HydroTable::echo() const {
    auto printVec = [&](const std::vector<double> &vec, const std::string &name) {
        std::cout << name << ": " << std::endl;
        for (auto &v : vec) {
            printf("%6g,", v);
        }
        printf("\n");
    };
    printVec(sloc, "sloc");
    printVec(slipVelocity, "slipVelocity");
    printVec(activeStress, "activeStress");
    printVec(radius, "radius");
}

void HydroTable::getValue(const int nPts, const double *sPtr, double *slipVelocityPtr, double *activeStressPtr,
                          double *radiusPtr) const {
    int nd[1];
    nd[0] = slipVelocity.size(); // data dimension
    // interpolate slipVelocity
    mlinterp::interp(nd, nPts, slipVelocity.data(), slipVelocityPtr, sloc.data(), sPtr);
    // interpolate activeStress
    mlinterp::interp(nd, nPts, activeStress.data(), activeStressPtr, sloc.data(), sPtr);
    // interpolate radius
    mlinterp::interp(nd, nPts, radius.data(), radiusPtr, sloc.data(), sPtr);
}