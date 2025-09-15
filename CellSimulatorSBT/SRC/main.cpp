#include <string>
#include <vector>

#include <mpi.h>

#include "CellSystem.hpp"

int main(int argc, char **argv) {

    MPI_Init(&argc, &argv);
    PS::Initialize(argc, argv);
    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    // create a system and distribute it to all ranks
    std::string runConfig = "RunConfig.yaml";
    std::string cellConfig = "CellConfig.yaml";
    std::string posFile = "CellInitial.dat";
    std::string restartFile = "TimeStepInfo.txt";
    {
        constexpr int maxQuadPt = 48; 
        CellSystem<maxQuadPt> system(runConfig, cellConfig, posFile, restartFile, argc, argv); // MPI is initialized inside PS::Initialize()
        // main time loop
        while (!system.stop()) {
            system.step();
        }
    }
    // mpi finalize
    // let the root rank wait for other
    PS::Finalize();
    MPI_Finalize();
    return 0;
}
