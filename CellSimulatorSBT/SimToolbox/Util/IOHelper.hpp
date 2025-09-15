/**
 * @file IOHelper.hpp
 * @author wenyan4work(wenyan4work@gmail.com)
 * @brief
 * @version 0.1
 * @date 2018-12-13
 *
 * @copyright Copyright (c) 2018
 *
 */
#ifndef IOHELPER_HPP_
#define IOHELPER_HPP_

#include "Base64.hpp"

#include <array>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <type_traits>
#include <vector>

#include <errno.h>
#include <sys/stat.h>

/**
 * @brief an abstract class for helping with VTK IO
 *
 */
class IOHelper {
  public:
    enum class IOTYPE { UInt8, Int32, Float32, Float64 }; ///< VTK XML binary data type

    struct FieldVTU {
        int dimension;
        IOTYPE type; // in most cases, choose between Float32 and Float64
        std::string name;

        FieldVTU(int dimension_, IOTYPE type_, std::string name_) : dimension(dimension_), type(type_), name(name_) {}
    };

    static std::string getTypeName(IOTYPE type) {
        std::string name;
        if (type == IOTYPE::UInt8) {
            name = "UInt8";
        } else if (type == IOTYPE::Int32) {
            name = "Int32";
        } else if (type == IOTYPE::Float32) {
            name = "Float32";
        } else if (type == IOTYPE::Float64) {
            name = "Float64";
        }
        return name;
    }

    static bool fileExist(const std::string &name) {
        std::ifstream f(name.c_str());
        return f.good();
    }

    static void makeSubFolder(const std::string folder) {

       const int err = mkdir(folder.c_str(), 0755); // make one folder at a time. parent folder must exist
        if (errno == EEXIST) {
            // printf("Directory already exists.\n");
            return;
        }
        if (err != 0) {
            printf("errno: %d \n", errno);
            printf("Error creating directory!\n");
            std::exit(1);
        }
    }

    /*******************************
     * VTK unstructured grid  *
     ********************************/

    static void writeHeadVTU(std::ofstream &vtkfile) {
        vtkfile << "<?xml version=\"1.0\"?>\n";
        vtkfile << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\"  "
                   "header_type=\"UInt32\">\n";
        vtkfile << "<UnstructuredGrid>\n";
    }

    static void writeTailVTU(std::ofstream &vtkfile) {
        vtkfile << "</UnstructuredGrid>\n";
        vtkfile << "</VTKFile>" << std::endl;
    }

    static void writePVTUFile(std::string filename, const std::vector<FieldVTU> &dataFields,
                              const std::vector<std::string> &pieceNames) {

        std::ofstream pvtufile(filename, std::ios::out);

        pvtufile << "<?xml version=\"1.0\"?>\n";
        pvtufile << "<VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" "
                    "header_type=\"UInt32\"> \n";
        pvtufile << "<PUnstructuredGrid GhostLevel=\"0\"> \n";

        pvtufile << "<PPointData Scalars=\"scalars\">\n";

        for (int i = 0; i < dataFields.size(); i++) {
            auto &data = dataFields[i];
            // data.first = dimension
            // data.second = name
            std::string type = getTypeName(data.type);
            pvtufile << "<PDataArray Name=\"" << data.name << "\" type=\"" << type << "\" NumberOfComponents=\""
                     << data.dimension << "\" format=\"binary\"/>\n";
        }

        pvtufile << "</PPointData>\n";

        pvtufile << "<PPoints> \n";
        pvtufile << "<PDataArray NumberOfComponents=\"3\" type=\"Float64\" format=\"binary\"/>\n";
        pvtufile << "</PPoints> \n";

        for (const auto &piece : pieceNames) {

            pvtufile << "<Piece Source=\"" << piece << "\"/>\n";
        }
        pvtufile << "</PUnstructuredGrid>\n";
        pvtufile << "</VTKFile>\n";
        pvtufile.close();
    }

    // overload, some are point data, some are cell data
    static void writePVTUFile(std::string filename, const std::vector<FieldVTU> &pointDataFields,
                              const std::vector<FieldVTU> &cellDataFields, const std::vector<std::string> &pieceNames) {

        std::ofstream pvtufile(filename, std::ios::out);

        pvtufile << "<?xml version=\"1.0\"?>\n";
        pvtufile << "<VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" "
                    "header_type=\"UInt32\"> \n";
        pvtufile << "<PUnstructuredGrid GhostLevel=\"0\"> \n";

        // point data
        pvtufile << "<PPointData Scalars=\"scalars\">\n";
        for (int i = 0; i < pointDataFields.size(); i++) {
            auto &data = pointDataFields[i];
            // data.first = dimension
            // data.second = name
            std::string type = getTypeName(data.type);
            pvtufile << "<PDataArray Name=\"" << data.name << "\" type=\"" << type << "\" NumberOfComponents=\""
                     << data.dimension << "\" format=\"binary\"/>\n";
        }
        pvtufile << "</PPointData>\n";

        // cell data
        pvtufile << "<PCellData Scalars=\"scalars\">\n";
        for (int i = 0; i < cellDataFields.size(); i++) {
            auto &data = cellDataFields[i];
            // data.first = dimension
            // data.second = name
            std::string type = getTypeName(data.type);
            pvtufile << "<PDataArray Name=\"" << data.name << "\" type=\"" << type << "\" NumberOfComponents=\""
                     << data.dimension << "\" format=\"binary\"/>\n";
        }
        pvtufile << "</PCellData>\n";

        pvtufile << "<PPoints> \n";
        pvtufile << "<PDataArray NumberOfComponents=\"3\" type=\"Float64\" format=\"binary\"/>\n";
        pvtufile << "</PPoints> \n";

        for (const auto &piece : pieceNames) {

            pvtufile << "<Piece Source=\"" << piece << "\"/>\n";
        }
        pvtufile << "</PUnstructuredGrid>\n";
        pvtufile << "</VTKFile>\n";
        pvtufile.close();
    }

    /*******************************
     * VTK poly data  *
     ********************************/

    static void writeHeadVTP(std::ofstream &vtkfile) {
        vtkfile << "<?xml version=\"1.0\"?>\n";
        vtkfile << "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\"  "
                   "header_type=\"UInt32\">\n";
        vtkfile << "<PolyData>\n";
    }

    static void writeTailVTP(std::ofstream &vtkfile) {
        vtkfile << "</PolyData>\n";
        vtkfile << "</VTKFile>" << std::endl;
    }

    // all data are point data
    static void writePVTPFile(std::string filename, const std::vector<FieldVTU> &dataFields,
                              const std::vector<std::string> &pieceNames) {

        std::ofstream pvtufile(filename, std::ios::out);

        pvtufile << "<?xml version=\"1.0\"?>\n";
        pvtufile << "<VTKFile type=\"PPolyData\" version=\"1.0\" byte_order=\"LittleEndian\" "
                    "header_type=\"UInt32\"> \n";
        pvtufile << "<PPolyData GhostLevel=\"0\"> \n";

        pvtufile << "<PPointData Scalars=\"scalars\">\n";

        for (auto &data : dataFields) {
            std::string type = getTypeName(data.type);
            pvtufile << "<PDataArray Name=\"" << data.name << "\" type=\"" << type << "\" NumberOfComponents=\""
                     << data.dimension << "\" format=\"binary\"/>\n";
        }

        pvtufile << "</PPointData>\n";

        pvtufile << "<PPoints> \n";
        pvtufile << "<PDataArray NumberOfComponents=\"3\" type=\"Float64\" format=\"binary\"/>\n";
        pvtufile << "</PPoints> \n";

        for (const auto &piece : pieceNames) {
            pvtufile << "<Piece Source=\"" << piece << "\"/>\n";
        }
        pvtufile << "</PPolyData>\n";
        pvtufile << "</VTKFile>\n";
        pvtufile.close();
    }

    // overload, some are point data, some are cell data
    static void writePVTPFile(std::string filename, const std::vector<FieldVTU> &pointDataFields,
                              const std::vector<FieldVTU> &cellDataFields, const std::vector<std::string> &pieceNames) {

        std::ofstream pvtufile(filename, std::ios::out);

        pvtufile << "<?xml version=\"1.0\"?>\n";
        pvtufile << "<VTKFile type=\"PPolyData\" version=\"1.0\" byte_order=\"LittleEndian\" "
                    "header_type=\"UInt32\"> \n";
        pvtufile << "<PPolyData GhostLevel=\"0\"> \n";

        pvtufile << "<PPointData Scalars=\"scalars\">\n";
        for (auto &data : pointDataFields) {
            std::string type = getTypeName(data.type);
            pvtufile << "<PDataArray Name=\"" << data.name << "\" type=\"" << type << "\" NumberOfComponents=\""
                     << data.dimension << "\" format=\"binary\"/>\n";
        }
        pvtufile << "</PPointData>\n";

        pvtufile << "<PCellData Scalars=\"scalars\">\n";
        for (auto &data : cellDataFields) {
            std::string type = getTypeName(data.type);
            pvtufile << "<PDataArray Name=\"" << data.name << "\" type=\"" << type << "\" NumberOfComponents=\""
                     << data.dimension << "\" format=\"binary\"/>\n";
        }
        pvtufile << "</PCellData>\n";

        pvtufile << "<PPoints> \n";
        pvtufile << "<PDataArray NumberOfComponents=\"3\" type=\"Float64\" format=\"binary\"/>\n";
        pvtufile << "</PPoints> \n";

        for (const auto &piece : pieceNames) {
            pvtufile << "<Piece Source=\"" << piece << "\"/>\n";
        }
        pvtufile << "</PPolyData>\n";
        pvtufile << "</VTKFile>\n";
        pvtufile.close();
    }

    /**
     * @brief write data array into binary XML data field
     *
     * @tparam T
     * @param data
     * @param name
     * @param numComp number of components (1 for scalar, 3 for vector, etc)
     * @param file
     */
    template <class T>
    static void writeDataArrayBase64(std::vector<T> &data, const std::string &name, int numComp, std::ofstream &file) {
        // set type name
        std::string vtktype;
        if (std::is_same<T, int>::value) {
            vtktype = "Int32";
        } else if (std::is_same<T, float>::value) {
            vtktype = "Float32";
        } else if (std::is_same<T, double>::value) {
            vtktype = "Float64";
        } else if (std::is_same<T, uint8_t>::value) {
            vtktype = "UInt8";
        }
        std::string contentB64;
        B64Converter::getBase64FromVector(data, contentB64);

        file << "<DataArray Name=\"" << name << "\" type=\"" << vtktype << "\" NumberOfComponents=\"" << numComp
             << "\" format=\"binary\">\n";
        file << contentB64 << "\n";
        file << "</DataArray>\n";
    }

    /*******************************
     * VTK rectilinear grid  *
     ********************************/
    static void writeHeadVTR(std::ofstream &vtkfile, int boxLow[3], int boxHigh[3]) {
        vtkfile << "<?xml version=\"1.0\"?>\n";
        vtkfile << "<VTKFile type=\"RectilinearGrid\" version=\"1.0\" byte_order=\"LittleEndian\"  "
                   "header_type=\"UInt32\">\n";
        vtkfile << "<RectilinearGrid WholeExtent=\"" << boxLow[0] << " " << boxHigh[0] << " " << boxLow[1] << " "
                << boxHigh[1] << " " << boxLow[2] << " " << boxHigh[2] << "\">\n";
    }

    static void writeTailVTR(std::ofstream &vtkfile) {
        vtkfile << "</RectilinearGrid>\n";
        vtkfile << "</VTKFile>" << std::endl;
    }

    /*******************************
     * VTK ImageData         *
     * equispacing rectilinear grid*
     ********************************/
    static void writeHeadVTI(std::ofstream &vtkfile, int boxLow[3], int boxHigh[3], double spacing[3],
                             double origin[3]) {
        vtkfile << "<?xml version=\"1.0\"?>\n";
        vtkfile << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\"  "
                   "header_type=\"UInt32\">\n";
        vtkfile << "<ImageData WholeExtent=\"" << boxLow[0] << " " << boxHigh[0] << " " << boxLow[1] << " "
                << boxHigh[1] << " " << boxLow[2] << " " << boxHigh[2] << "\""                //
                << " Origin=\" " << origin[0] << " " << origin[1] << " " << origin[2] << "\"" //
                << " Spacing=\" " << spacing[0] << " " << spacing[1] << " " << spacing[2]     //
                << "\">\n";
    }

    static void writeTailVTI(std::ofstream &vtkfile) {
        vtkfile << "</ImageData>\n";
        vtkfile << "</VTKFile>" << std::endl;
    }
};

#endif