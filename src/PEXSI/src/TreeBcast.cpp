#include <mpi.h>
#include <vector>
#include <map>
//#include "pexsi/TreeBcast.hpp"
namespace PEXSI{
std::map< MPI_Comm , std::vector<int> > commGlobRanks;
}
