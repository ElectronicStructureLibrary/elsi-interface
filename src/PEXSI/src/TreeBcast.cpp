#include <vector>
#include <map>
#include <mpi.h>
//#include "pexsi/TreeBcast.hpp"
namespace PEXSI{
std::map< MPI_Comm , std::vector<int> > commGlobRanks;
}
