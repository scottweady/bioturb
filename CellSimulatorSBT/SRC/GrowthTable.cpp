#include "GrowthTable.h"

#include "mlinterp.hpp"
#include "rapidcsv.h"

#include <cstdio>
#include <iostream>
#include <set>
#include <string>

GrowthTable::GrowthTable(const std::string &csvFile, const std::string &xname_, const std::string &yname_)
    : xname(xname_), yname(yname_) {
    // expect a csv file
    // xname,yname,f
    // x0,y0,f00
    // x0,y1,f01
    // x0,y2,f02
    // ...
    // x1,y0,f10
    // x1,y1,f11
    // x1,y2,f12
    rapidcsv::Document doc(csvFile, rapidcsv::LabelParams(0, -1));
    std::vector<std::string> columnNames = doc.GetColumnNames();
    // for (auto &c : columnNames) {
    //     std::cout << c << std::endl;
    // }
    std::vector<double> x = doc.GetColumn<double>(xname);
    std::vector<double> y = doc.GetColumn<double>(yname);
    std::vector<double> f = doc.GetColumn<double>("f");
    // determine nx and ny according to unique entries
    auto getGrid = [&](const std::vector<double> &col) {
        std::set<double> val;
        for (auto &v : col) {
            val.insert(v);
        }
        std::vector<double> grid;
        for (auto &v : val) {
            grid.push_back(v);
        }
        std::sort(grid.begin(), grid.end());
        return grid;
    };

    xd = getGrid(x);
    yd = getGrid(y);
    const int nxd = xd.size();
    const int nyd = yd.size();
    const int nf = nxd * nyd;
    if (f.size() != nf) {
        printf("table size error: %d\n", nf);
        exit(1);
    }
    fd = std::move(f); // f must be in the correct order (row major) in file
}

std::vector<double> GrowthTable::getValue(const std::vector<double> &xi, const std::vector<double> &yi) const {
    const int ni = xi.size();
    if (ni != yi.size()) {
        printf("interpolation request size error\n");
        exit(1);
    }
    std::vector<double> fi(ni);

    int nd[2];
    nd[0] = xd.size();
    nd[1] = yd.size(); // data dimension
    mlinterp::interp(nd, ni, fd.data(), fi.data(), xd.data(), xi.data(), yd.data(), yi.data());
    return fi;
}

void GrowthTable::echo() const {
    auto printVec = [&](const std::vector<double> &vec, const std::string &name) {
        std::cout << name << std::endl;
        for (auto &v : vec) {
            printf("%6g,", v);
        }
        printf("\n");
    };
    printVec(xd, xname + " grid");
    printVec(yd, yname + " grid");
    printVec(fd, "f grid");
}