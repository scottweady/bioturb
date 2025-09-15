#include "Config.h"

#include "Util/YamlHelper.hpp"

#include <exception>
#include <iostream>

void Growth::echo() const {
    printf("-------------------------------------------\n");
    printf("tauD: %g\n", tauD);
    printf("Delta: %g\n", Delta);
    printf("sigma: %g\n", sigma);
    cellGrowthTable.echo();
    printf("-------------------------------------------\n");
}

void Hydro::echo() const {
    printf("-------------------------------------------\n");
    printf("Number of Quadrature Points: %d\n", numberOfQuadraturePoints);
    printf("forceExternal: [%lf, %lf, %lf]\n", forceExternal[0], forceExternal[1], forceExternal[2]);
    printf("torqueExternal: [%lf, %lf, %lf]\n", torqueExternal[0], torqueExternal[1], torqueExternal[2]);
    cellHydroTable.echo();
    printf("-------------------------------------------\n");
}

Config::Config(std::string fileName) {
    try {
        YAML::Node config = YAML::LoadFile(fileName);
        readConfig(config, VARNAME(hydro), hydro, "");
        readConfig(config, VARNAME(freeStreamVelocity), freeStreamVelocity, 3, "");
        readConfig(config, VARNAME(fmmMultOrder), fmmMultOrder, "");
        readConfig(config, VARNAME(flowGrid), flowGrid, "");

        readConfig(config, VARNAME(wallZ0), wallZ0, "");

        // load cell growth table if any
        if (config[VARNAME(cellGrowth)]) {
            YAML::Node cellGrowthtables = config[VARNAME(cellGrowth)];
            printf("reading cell growth settings\n");
            for (const auto &g : cellGrowthtables) {
                const int group = g["group"].as<int>();
                const double time = g["tauD"].as<double>();
                const double length = g["Delta"].as<double>();
                const double sigma = g["sigma"].as<double>();
                const std::string &csvFileg = g["file"].as<std::string>();
                // printf("%d,%g,%g,%g,\n", group, time, length, sigma);
                // std::cout << csvFile << std::endl;
                GrowthTable tGrowth(csvFileg, "z/H", "t/tauD");
                cellGrowth.insert({0, Growth(time, length, sigma, tGrowth)});
            }
        }

        // load cell hydro table if any
        if (config[VARNAME(cellHydro)]) {
            YAML::Node cellHydrotables = config[VARNAME(cellHydro)];
            printf("reading cell hydro settings\n");
            for (const auto &h : cellHydrotables) {
                const int group = h["group"].as<int>();
                const int numQuadPt = h["numberOfQuadraturePoints"].as<int>();
                const double forceExternal[3] = {h["forceExternal"][0].as<double>(), h["forceExternal"][1].as<double>(),
                                                 h["forceExternal"][2].as<double>()};
                const double torqueExternal[3] = {h["torqueExternal"][0].as<double>(),
                                                  h["torqueExternal"][1].as<double>(),
                                                  h["torqueExternal"][2].as<double>()};
                const std::string &csvFileh = h["file"].as<std::string>();
                HydroTable tHydro(csvFileh);
                cellHydro.insert({group, Hydro(numQuadPt, forceExternal, torqueExternal, tHydro)});
            }
        }

    } catch (const std::exception &exc) {
        std::cout << "Caught an exception during reading the Config file with name '" << fileName << "'." << std::endl;
        std::cout << "The exception was: " << exc.what() << std::endl;
        std::cout << "Ignoring the missing or bad cellConfig.yaml file" << std::endl;
    }

    return;
}

void Config::echo() const {
    // input correctness check
    if (hydro) {
        printf("Hydro Turned ON. Setting: \n");
        printf("Free-stream Velocity: [%lf, %lf, %lf]\n", freeStreamVelocity[0], freeStreamVelocity[1],
               freeStreamVelocity[2]);
        printf("FMM MultOrder: %d\n", fmmMultOrder);
        if (wallZ0) {
            printf("With No-Slip Wall\n");
        }
        if (flowGrid > 0) {
            printf("Dump flow with grid spacing: %lf\n", flowGrid);
        } else {
            printf("No flow dump\n");
        }
    }

    {
        if ((cellGrowth.size() == 1) && (cellHydro.size() == 1)) {
            printf("WARNING: only 1 species supplied, used for all cell groups\n");
        }
        for (auto &v : cellGrowth) {
            std::cout << "cellGrowth group id: " << v.first << std::endl;
            v.second.echo();
        }
        for (auto &v : cellHydro) {
            std::cout << "cellHydro group id: " << v.first << std::endl;
            v.second.echo();
        }
    }
}