/*
   Copyright (c) 2012 The Regents of the University of California,
   through Lawrence Berkeley National Laboratory.

Authors: Lin Lin and Mathias Jacquelin

This file is part of PEXSI. All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

(1) Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.
(2) Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.
(3) Neither the name of the University of California, Lawrence Berkeley
National Laboratory, U.S. Dept. of Energy nor the names of its contributors may
be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

You are under no obligation whatsoever to provide any bug fixes, patches, or
upgrades to the features, functionality or performance of the source code
("Enhancements") to anyone; however, if you choose to make your Enhancements
available either publicly, or directly to Lawrence Berkeley National
Laboratory, without imposing a separate written license agreement for such
Enhancements, then you hereby grant the following license: a non-exclusive,
royalty-free perpetual license to install, use, modify, prepare derivative
works, incorporate into other computer software, distribute, and sublicense
such enhancements or derivative works thereof, in binary and source code form.
 */
/// @file pselinv.hpp
/// @brief Main file for parallel selected inversion.
/// @date 2013-08-05
#ifndef _PEXSI_PSELINV_HPP_
#define _PEXSI_PSELINV_HPP_

// *********************************************************************
//  Common utilities
// *********************************************************************

#include "pexsi/environment.hpp"
#include "pexsi/NumVec.hpp"
#include "pexsi/NumMat.hpp"
#include "pexsi/sparse_matrix.hpp"

#include "pexsi/superlu_dist_interf.hpp"
#include "pexsi/mpi_interf.hpp"
#include "pexsi/utility.hpp"
#include "pexsi/blas.hpp"
#include "pexsi/lapack.hpp"

#include "pexsi/TreeBcast.hpp"

#include	"pexsi/TreeBcast_v2.hpp"
#include	"pexsi/TreeReduce_v2.hpp"

#include <set>


//#define IDX_TO_TAG(lidx,tag) (SELINV_TAG_COUNT*(lidx)+(tag))
#define sym_IDX_TO_TAG( lidx, tag, numSuper, max)  ((SELINV_TAG_COUNT)*(numSuper)*((lidx)%((max)+1))+(tag))

#define IDX_TO_TAG(lidx,tag,max) ((SELINV_TAG_COUNT*(lidx%(max+1))+(tag)))
#define IDX_TO_TAG2(sidx,lidx,tag) (SELINV_TAG_COUNT*(sidx)+(tag))
#define TAG_TO_IDX(tag,typetag) (((tag)-(typetag))/SELINV_TAG_COUNT)


//#define LIST_BARRIER
//#define ALL_BARRIER

namespace PEXSI{

enum MSGTYPE {LSIZE=0,LROWSIZE,USIZE,UCOLSIZE,LCONTENT,LROWCONTENT,UCONTENT,UCOLCONTENT,MSGCOUNT};

struct ULComparator {
  IntNumVec & lookup;
  ULComparator(IntNumVec & v):lookup(v){}
  bool operator() (int i,int j) { return ((lookup)(i)<(lookup)(j));}
};



typedef std::vector<bool> bitMask;
typedef std::map<bitMask , std::vector<Int> > bitMaskSet;





/// @struct PSelInvOptions
/// @brief A thin interface for passing parameters to set the PSelInv
/// options.
///
struct PSelInvOptions{
  /// @brief The maximum pipeline depth.
  Int              maxPipelineDepth;

  /// @brief Use symmetric storage for the selected inversion or not.
  Int              symmetricStorage;

  // Member functions to setup the default value
  PSelInvOptions(): maxPipelineDepth(-1), symmetricStorage(0) {}
};


/// @struct FactorizationOptions
/// @brief A thin interface for passing parameters to set the Factorization
/// options.
///
struct FactorizationOptions{
  std::string ColPerm;
  std::string RowPerm;
  Int Symmetric;

  // Member functions to setup the default value
  FactorizationOptions(): Symmetric(1) {}
};





/**********************************************************************
 * Basic PSelInv data structure
 **********************************************************************/

/// @struct GridType
///
/// @brief GridType is the PSelInv way of defining the grid.
///
/// GridType should be consistent with the grid used by SuperLU.
///
/// NOTE: It is your responsibility to make sure that the SuperLUGrid
/// and GridType used for SelInv are the same.
struct GridType{
  // Data
  MPI_Comm    comm;
  MPI_Comm    rowComm;
  MPI_Comm    colComm;
  Int         mpirank;
  Int         mpisize;
  Int         numProcRow;
  Int         numProcCol;

  // Member function
  GridType( MPI_Comm Bcomm, int nprow, int npcol );
  ~GridType();
};

/// @struct SuperNodeType
///
/// @brief SuperNodeType describes mapping between supernode and column, the
/// permutation information, and potentially the elimination tree (not
/// implemented here).
///
/// superIdx[i] is the supernode index to which column i belongs.
/// This is the same as supno[i] in SuperLU.
///
/// superPtr[s] is the leading column of the s-th supernode (as in
/// colptr).  This is the same as xsup[s] in SuperLU.
///
///	e.g.   superIdx  0 1 2 2 3 3 3 4 4 4 4 4   (n=12)
///	       superPtr  0 1 2 4 7 12
///
/// This is allocated during symbolic factorization SYMBFACT.
///
/// perm is the permutation vector.  Symmetric permutation is assumed.
/// perm is the same as ScalePermstruct -> perm_c.
///
/// permInv is the inverse of the permutation vector.
///
struct SuperNodeType{
  IntNumVec   perm;
  IntNumVec   permInv;
  IntNumVec   perm_r;
  IntNumVec   permInv_r;
  IntNumVec   superIdx;
  IntNumVec   superPtr;
  IntNumVec   etree;
};


/// @struct LBlock
///
/// @brief LBlock stores a nonzero block in the lower triangular part or
/// the diagonal part in PSelInv.
template<typename T>
struct LBlock{
  // Variables
  /// @brief Block index (supernodal index)
  Int               blockIdx;

  /// @brief Number of nonzero rows.
  Int               numRow;

  /// @brief Number of nonzero columns.
  Int               numCol;

  /// @brief Dimension numRow * 1, index (0-based) for the number of nonzero rows.
  IntNumVec         rows;


  /// @brief Dimension numRow * numCol, nonzero elements.
  NumMat<T>    nzval;

  // Member functions;
  LBlock() {
    blockIdx = -1; numRow = 0; numCol =0;
    nzval.Resize(0,0);
  }
  ~LBlock() {}
  LBlock& operator = (const LBlock& LB) {
    blockIdx    = LB.blockIdx;
    numRow      = LB.numRow;
    numCol      = LB.numCol;
    rows        = LB.rows;
    nzval       = LB.nzval;
    return *this;
  }
  friend std::ostream& operator<<(std::ostream& out, const LBlock& vec) // output
  {
    out << "(" << vec.blockIdx << ", " << vec.numRow << ", " << vec.numCol <<std::endl<< "rows " << vec.rows <<std::endl<< "nzval " <<std::endl<< vec.nzval << ")";
    return out;
  }


};

/// @struct UBlock
///
/// @brief UBlock stores a nonzero block in the upper triangular part in PSelInv.
///
/// In particular, the current version of PSelInv is for sparse
/// symmetric matrices.  All UBlocks, labeled as U(i,j), i<j save the
/// redundant information as saved in L(j,i). The purpose of having
/// UBlocks is to facilitate the communication.
///
/// @see PMatrix::SelInv
template<typename T>
struct UBlock{
  // Variables
  /// @brief Block index (supernodal index)
  Int               blockIdx;

  /// @brief Number of nonzero rows.
  Int               numRow;

  /// @brief Number of nonzero columns.
  Int               numCol;

  /// @brief Dimension numCol * 1, index (0-based) for the number of nonzero rows.
  IntNumVec         cols;

  /// @brief Dimension numRow * numCol, nonzero elements.
  NumMat<T>    nzval;

  // Member functions;
  UBlock() {
    blockIdx = -1; numRow = 0; numCol =0;
  }
  ~UBlock() {}
  UBlock& operator = (const UBlock& UB) {
    blockIdx    = UB.blockIdx;
    numRow      = UB.numRow;
    numCol      = UB.numCol;
    cols        = UB.cols;
    nzval       = UB.nzval;
    return *this;
  }

  friend std::ostream& operator<<(std::ostream& out, const UBlock& vec) // output
  {
    out << "(" << vec.blockIdx << ", " << vec.numRow << ", " << vec.numCol <<std::endl<< "cols " << vec.cols <<std::endl<< "nzval " <<std::endl<< vec.nzval << ")";
    return out;
  }
};

// *********************************************************************
// SuperLU style utility functions
//
// The SuperLU style macros are defined here as inline functions
// so that the code is more portable.
// *********************************************************************

/// @brief MYPROC returns the current processor rank.
inline Int MYPROC( const GridType* g )
{ return g->mpirank; }

/// @brief MYROW returns my processor row
inline Int MYROW( const GridType* g )
{ return g->mpirank / g->numProcCol; }

/// @brief MYCOL returns my processor column
inline Int MYCOL( const GridType* g )
{ return g->mpirank % g->numProcCol; }

/// @brief PROW returns the processor row that the bnum-th block
/// (supernode) belongs to.
inline Int PROW( Int bnum, const GridType* g )
{ return bnum % g->numProcRow; }

/// @brief PCOL returns the processor column that the bnum-th block
/// (supernode) belongs to.
inline Int PCOL( Int bnum, const GridType* g )
{ return bnum % g->numProcCol; }

/// @brief PNUM returns the processor rank that the bnum-th block
/// (supernode) belongs to.
inline Int PNUM( Int i, Int j, const GridType* g )
{ return  (i%g->numProcRow) * g->numProcCol + j%g->numProcCol; }

/// @brief LBi returns the local block number on the processor at
/// processor row PROW( bnum, g ).
inline Int LBi( Int bnum, const GridType* g )
{ return bnum / g->numProcRow; }

/// @brief LBj returns the local block number on the processor at
/// processor column PCOL( bnum, g ).
inline Int LBj( Int bnum, const GridType* g)
{ return bnum / g->numProcCol; }

/// @brief GBi returns the global block number from a local block number
/// in the row direction.
inline Int GBi( Int iLocal, const GridType* g )
{ return iLocal * g->numProcRow + MYROW( g ); }

/// @brief GBj returns the global block number from a local block number
/// in the column direction.
inline Int GBj( Int jLocal, const GridType* g )
{ return jLocal * g->numProcCol + MYCOL( g ); }

/// @brief CEILING is used for computing the storage space for local
/// number of blocks.
inline Int CEILING( Int a, Int b )
{ return (a%b) ? ( a/b + 1 ) : ( a/b ); }

/// @brief BlockIdx returns the block index of a column i.
inline Int BlockIdx( Int i, const SuperNodeType *s )
{ return s->superIdx[i]; }

/// @brief FirstBlockCol returns the first column of a block
/// bnum.
inline Int FirstBlockCol( Int bnum, const SuperNodeType *s )
{ return s->superPtr[bnum]; }


/// @brief FirstBlockRow returns the first column of a block
/// bnum. Note: the functionality of FirstBlockRow is exactly the same
/// as in FirstBlockCol.
inline Int FirstBlockRow( Int bnum, const SuperNodeType *s )
{ return s->superPtr[bnum]; }


/// @brief SuperSize returns the size of the block bnum.
inline Int SuperSize( Int bnum, const SuperNodeType *s )
{ return s->superPtr[bnum+1] - s->superPtr[bnum]; }

/// @brief NumSuper returns the total number of supernodes.
inline Int NumSuper( const SuperNodeType *s )
{ return s->superPtr.m() - 1; }

/// @brief NumCol returns the total number of columns for a supernodal
/// partiiton.
inline Int NumCol( const SuperNodeType *s )
{ return s->superIdx.m(); }


// *********************************************************************
// Serialize / Deserialize
// *********************************************************************

// L part

/// @namespace LBlockMask
///
/// @brief LBlockMask allows one to compress the selected data in
/// LBlock used for communication.
///
/// Example
/// -------
/// std::vector<Int> mask( LBlockMask::TOTAL_NUMBER, 1 );
///
/// assumes all information is to be communicated.
///
///	mask[LBlockMask::NZVAL] = 0;
///
///	then nzval will not be serialized / deserialized.
///
namespace LBlockMask{
enum {
  BLOCKIDX,
  NUMROW,
  NUMCOL,
  ROWS,
  NZVAL,
  TOTAL_NUMBER
};
}


template<typename T>
Int inline serialize(LBlock<T>& val, std::ostream& os, const std::vector<Int>& mask){
  if(mask[LBlockMask::BLOCKIDX]==1) serialize(val.blockIdx, os, mask);
  if(mask[LBlockMask::NUMROW  ]==1) serialize(val.numRow,  os, mask);
  if(mask[LBlockMask::NUMCOL  ]==1) serialize(val.numCol,  os, mask);
  if(mask[LBlockMask::ROWS    ]==1) serialize(val.rows, os, mask);
  if(mask[LBlockMask::NZVAL   ]==1) serialize(val.nzval, os, mask);
  return 0;
}

template<typename T>
Int inline deserialize(LBlock<T>& val, std::istream& is, const std::vector<Int>& mask){
  if(mask[LBlockMask::BLOCKIDX]==1) deserialize(val.blockIdx, is, mask);
  if(mask[LBlockMask::NUMROW  ]==1) deserialize(val.numRow,  is, mask);
  if(mask[LBlockMask::NUMCOL  ]==1) deserialize(val.numCol,  is, mask);
  if(mask[LBlockMask::ROWS    ]==1) deserialize(val.rows,   is, mask);
  if(mask[LBlockMask::NZVAL   ]==1) deserialize(val.nzval,  is, mask);
  return 0;
}

// U part
/// @namespace UBlockMask
///
/// @brief UBlockMask allows one to compress the selected data in
/// UBlock used for communication.
///
/// Example
/// -------
///
/// std::vector<Int> mask( UBlockMask::TOTAL_NUMBER, 1 );
///
/// assumes all information is to be communicated.
///
///	mask[UBlockMask::NZVAL] = 0;
///
///	then nzval will not be serialized / deserialized.
namespace UBlockMask{
enum {
  BLOCKIDX,
  NUMROW,
  NUMCOL,
  COLS,
  NZVAL,
  TOTAL_NUMBER
};
}

template<typename T>
Int inline serialize(UBlock<T>& val, std::ostream& os, const std::vector<Int>& mask){
  Int i = 0;
  if(mask[i]==1) serialize(val.blockIdx, os, mask); i++;
  if(mask[i]==1) serialize(val.numRow,  os, mask); i++;
  if(mask[i]==1) serialize(val.numCol,  os, mask); i++;
  if(mask[i]==1) serialize(val.cols, os, mask);   i++;
  if(mask[i]==1) serialize(val.nzval, os, mask);  i++;
  return 0;
}

template<typename T>
Int inline deserialize(UBlock<T>& val, std::istream& is, const std::vector<Int>& mask){
  Int i = 0;
  if(mask[i]==1) deserialize(val.blockIdx, is, mask); i++;
  if(mask[i]==1) deserialize(val.numRow,  is, mask); i++;
  if(mask[i]==1) deserialize(val.numCol,  is, mask); i++;
  if(mask[i]==1) deserialize(val.cols,   is, mask); i++;
  if(mask[i]==1) deserialize(val.nzval,  is, mask); i++;
  return 0;
}

/**********************************************************************
 * Main data structure in PSelInv: PMatrix
 **********************************************************************/

/// @class PMatrix
///
/// @brief PMatrix contains the main data structure and the
/// computational routine for the parallel selected inversion.
///
/// **NOTE** The following is a bit obsolete.
///
/// Procedure for Selected Inversion
/// --------------------------------
///
/// After factorizing a SuperLUMatrix luMat (See SuperLUMatrix for
/// information on how to perform factorization), perform the following
/// steps for parallel selected inversion.
///
/// - Conversion from SuperLU_DIST.
///
///   Symbolic information
///
///       SuperNodeType super;
///       PMatrix PMloc;
///       luMat.SymbolicToSuperNode( super );
///
///   Numerical information, both L and U.
///
///       luMat.LUstructToPMatrix( PMloc );
///
/// - Preparation.
///
///   Construct the communication pattern for SelInv.
///
///       PMloc.ConstructCommunicationPattern();
///       or PMloc.ConstructCommunicationPattern_P2p();
///       or PMloc.ConstructCommunicationPattern_Collectives();
///
///   Numerical preparation so that SelInv only involves Gemm.
///
///       PMloc.PreSelInv();
///
/// - Selected inversion.
///
///       PMloc.SelInv();
///       or PMloc.SelInv_P2p();
///       or PMloc.SelInv_Collectives();
///
/// - Postprocessing.
///
///   Get the information in DistSparseMatrix format
///
///       DistSparseMatrix<Scalar> Ainv;
///       PMloc.PMatrixToDistSparseMatrix( Ainv );
///
/// Note
/// ----
///
/// - All major operations of PMatrix, including the selected inversion
/// are defined directly as the member function of PMatrix.
///
/// - In the current version of PMatrix, square grid is assumed.  This
/// assumption is only used when sending the information to
/// cross-diagonal blocks, i.e. from L(isup, ksup) to U(ksup, isup).
/// This assumption can be relaxed later.

template<typename T>
class PMatrix{
public:
  /// @brief Create is a factory method which returns a pointer either to a new PMatrix or to a new PMatrixUnsym
  /// depending on the pLuOpt parameter.

  static PMatrix<T> * Create(const GridType * pGridType, const SuperNodeType * pSuper, const PSelInvOptions * pSelInvOpt , const FactorizationOptions * pFactOpt);
  static PMatrix<T> * Create(const FactorizationOptions * pFactOpt);

public:
  // This is the tag used for mpi communication for selinv

  enum{
    /**/        SELINV_TAG_U_SIZE,
    SELINV_TAG_U_CONTENT,
    /**/        SELINV_TAG_L_SIZE,
    SELINV_TAG_L_CONTENT,
    SELINV_TAG_L_REDUCE,
    /**/        SELINV_TAG_D_SIZE,
    /**/        SELINV_TAG_D_CONTENT,
    SELINV_TAG_D_REDUCE,
    SELINV_TAG_L_SIZE_CD,
    SELINV_TAG_L_CONTENT_CD,
    SELINV_TAG_COUNT
  };


  std::vector<std::vector<Int> > ColBlockIdx_;
  std::vector<std::vector<Int> > RowBlockIdx_;
protected:

  virtual void PMatrixToDistSparseMatrix_ ( const NumVec<Int> & AcolptrLocal, const NumVec<Int> & ArowindLocal, const Int Asize, const LongInt Annz, const Int AnnzLocal, DistSparseMatrix<T>& B );

  // *********************************************************************
  // Variables
  // *********************************************************************
  // Data variables

  const GridType*       grid_;
  const SuperNodeType*  super_;
  const PSelInvOptions * options_;
  const FactorizationOptions * optionsFact_;

  Int limIndex_;
  Int maxTag_;

  std::vector<std::vector<LBlock<T> > > L_;
  std::vector<std::vector<UBlock<T> > > U_;

  std::vector<std::vector<Int> > workingSet_;


  // Communication variables
  BolNumMat                       isSendToBelow_;
  BolNumMat                       isSendToRight_;
  BolNumVec                       isSendToDiagonal_;
  BolNumMat                       isSendToCrossDiagonal_;

  BolNumMat                       isRecvFromBelow_;
  BolNumVec                       isRecvFromAbove_;
  BolNumVec                       isRecvFromLeft_;
  BolNumMat                       isRecvFromCrossDiagonal_;


  std::vector<std::shared_ptr<TreeReduce_v2<T> > > redDTree2_;

  std::vector<std::shared_ptr<TreeBcast_v2<char> > > bcastLDataTree_;
  std::vector<std::shared_ptr<TreeReduce_v2<T> > > redLTree2_;

  std::vector<TreeBcast *> fwdToBelowTree_;
  std::vector<TreeBcast *> fwdToRightTree_;
  std::vector<TreeReduce<T> *> redToLeftTree_;
  std::vector<TreeReduce<T> *> redToAboveTree_;

  //This is for the symmetric storage implementation
  std::vector<Int> snodeEtree_;
  std::vector<Int> snodeTreeOffset_;
  std::vector<std::map<Int,Int> > snodeBlkidxToTree_;
  std::vector<std::vector<Int> > snodeTreeToBlkidx_;
  std::list<Int> syncPoints_;

//  double localFlops_;

  struct SuperNodeBufferType{
    //This is for the symmetric storage implementation
    std::vector<NumMat<T> > LUpdateBufBlk;
    std::vector< std::vector<char> > SstrLcolSendBlk;
    std::vector< Int > SizeSstrLcolSendBlk;

    NumMat<T>    LUpdateBuf;
    std::vector<char> SstrLcolSend;
    Int               SizeSstrLcolSend;
    std::vector<char> SstrUrowSend;
    std::vector<char> SstrLcolRecv;
    std::vector<char> SstrUrowRecv;
    Int               SizeSstrUrowSend;
    Int               SizeSstrLcolRecv;
    Int               SizeSstrUrowRecv;

    NumMat<T>    DiagBuf;
    std::vector<Int>  RowLocalPtr;
    std::vector<Int>  BlockIdxLocal;
    Int               Index;
    Int               Rank;
    Int               isReady;


    SuperNodeBufferType():
      SizeSstrLcolSend(0),
      SizeSstrUrowSend(0),
      SizeSstrLcolRecv(0),
      SizeSstrUrowRecv(0),
      Index(0),
      Rank(0),
      isReady(0){}

    SuperNodeBufferType(Int &pIndex) :
      SuperNodeBufferType(){
        Index = pIndex;
      }

  };


  /// @brief SelInvIntra_P2p
  inline void SelInvIntra_P2p(Int lidx,Int & rank);

  /// @brief SelInv_lookup_indexes
  inline void SelInv_lookup_indexes(SuperNodeBufferType & snode, std::vector<LBlock<T> > & LcolRecv, std::vector<UBlock<T> > & UrowRecv, NumMat<T> & AinvBuf,NumMat<T> & UBuf);
  inline void SelInv_lookup_indexes_seq(SuperNodeBufferType & snode, std::vector<LBlock<T> > & LcolRecv, std::vector<UBlock<T> > & UrowRecv, NumMat<T> & AinvBuf,NumMat<T> & UBuf);

  /// @brief GetWorkSet
  inline void GetWorkSet(std::vector<Int> & snodeEtree, std::vector<std::vector<Int> > & WSet);

  /// @brief UnpackData
  inline void UnpackData(SuperNodeBufferType & snode, std::vector<LBlock<T> > & LcolRecv, std::vector<UBlock<T> > & UrowRecv);

  /// @brief ComputeDiagUpdate
  inline void ComputeDiagUpdate(SuperNodeBufferType & snode);

  /// @brief SendRecvCD_UpdateU
  inline void SendRecvCD_UpdateU(std::vector<SuperNodeBufferType > & arrSuperNodes, Int stepSuper);

public:
  // *********************************************************************
  // Public member functions
  // *********************************************************************

//  double GetTotalFlops();

  PMatrix();

  PMatrix( const GridType* g, const SuperNodeType* s, const PEXSI::PSelInvOptions * o, const PEXSI::FactorizationOptions * oFact  );

  void deallocate();


  virtual ~PMatrix();
  PMatrix( const PMatrix & C);
  PMatrix & operator = ( const PMatrix & C);

  void Setup( const GridType* g, const SuperNodeType* s, const PEXSI::PSelInvOptions * o, const PEXSI::FactorizationOptions * oFact  );

  const PSelInvOptions * Options() const { return options_; }

  Int NumCol() const { return super_ -> superIdx.m(); }

  Int NumSuper() const { return super_ ->superPtr.m() - 1; }

  /// @brief NumLocalBlockCol returns the total number of block columns.
  Int NumLocalBlockCol() const { return CEILING( this->NumSuper(), grid_->numProcCol ); }

  /// @brief NumLocalBlockRow returns the total number of block rows.
  Int NumLocalBlockRow() const { return CEILING( this->NumSuper(), grid_->numProcRow); }


  std::vector< std::vector<Int> > & ColBlockIdx() { return ColBlockIdx_; }
  std::vector< std::vector<Int> > & RowBlockIdx() { return RowBlockIdx_; }
  std::vector<Int> & ColBlockIdx(Int jLocal) { return ColBlockIdx_[jLocal]; }
  std::vector<Int> & RowBlockIdx(Int iLocal) { return RowBlockIdx_[iLocal]; }


  /// @brief NumBlockL returns the number of nonzero L blocks for the
  /// local block column jLocal.
  Int NumBlockL( Int jLocal ) const { return L_[jLocal].size(); }

  /// @brief NumBlockU returns the number of nonzero U blocks for the
  /// local block row iLocal.
  Int NumBlockU( Int iLocal ) const { return U_[iLocal].size(); }

  /// @brief Grid returns the GridType structure of the current PMatrix.
  const GridType* Grid() const { return grid_; }

  /// @brief SuperNode returns the supernodal partition of the current
  /// PMatrix.
  const SuperNodeType* SuperNode() const { return super_; }

  /// @brief L returns the vector of nonzero L blocks for the local
  /// block column jLocal.
  std::vector<LBlock<T> >& L( Int jLocal ) { return L_[jLocal]; }

  /// @brief U returns the vector of nonzero U blocks for the local
  /// block row iLocal.
  std::vector<UBlock<T> >& U( Int iLocal ) { return U_[iLocal]; }

  /// @brief WorkingSet returns the ordered list of supernodes which could
  /// be done in parallel.
  std::vector<std::vector<int> >& WorkingSet( ) { return workingSet_; }
  //void WorkingSet(const std::vector<std::vector<int> > * pWSet,const std::vector<std::vector<int> > * pWRanks) { pWSet = &workingSet_; pWRanks = &workingRanks_; }

  /// @brief CountSendToRight returns the number of processors
  /// to the right of current processor with which it has to communicate
  Int CountSendToRight(Int ksup) {  Int count= std::count (isSendToRight_.VecData(ksup), isSendToRight_.VecData(ksup) + grid_->numProcCol, true); return (isSendToRight_(MYCOL(grid_),ksup)?count-1:count); }

  /// @brief CountSendToBelow returns the number of processors
  /// below current processor with which it has to communicate
  Int CountSendToBelow(Int ksup) {  Int count= std::count (isSendToBelow_.VecData(ksup), isSendToBelow_.VecData(ksup) + grid_->numProcRow, true); return (isSendToBelow_(MYROW(grid_),ksup)?count-1:count); }

  /// @brief CountRecvFromBelow returns the number of processors
  /// below the current processor from which it receives data
  Int CountRecvFromBelow(Int ksup) {  Int count= std::count (isRecvFromBelow_.VecData(ksup), isRecvFromBelow_.VecData(ksup) + grid_->numProcRow, true); return (isRecvFromBelow_(MYROW(grid_),ksup)?count-1:count); }

  /// @brief CountSendToCrossDiagonal returns the number of cross diagonal
  /// processors with which current processor has to communicate
  Int CountSendToCrossDiagonal(Int ksup) {  Int count= std::count (isSendToCrossDiagonal_.VecData(ksup), isSendToCrossDiagonal_.VecData(ksup) + grid_->numProcCol, true);  return ((isSendToCrossDiagonal_(MYCOL(grid_),ksup) && MYROW(grid_)==PROW(ksup,grid_))?count-1:count); }

  /// @brief CountRecvFromCrossDiagonal returns the number of cross diagonal
  /// processors with which current processor has to communicate
  Int CountRecvFromCrossDiagonal(Int ksup) {  Int count= std::count (isRecvFromCrossDiagonal_.VecData(ksup), isRecvFromCrossDiagonal_.VecData(ksup) + grid_->numProcRow, true);  return ((isRecvFromCrossDiagonal_(MYROW(grid_),ksup) && MYCOL(grid_)==PCOL(ksup,grid_))?count-1:count); }




  /// @brief GetEtree computes the supernodal elimination tree
  /// to be used later in the pipelined selected inversion stage.
  void GetEtree(std::vector<Int> & etree_supno );


  /// @brief ConstructCommunicationPattern constructs the communication
  /// pattern to be used later in the selected inversion stage.
  /// The supernodal elimination tree is used to add an additional level of parallelism between supernodes.
  /// [ConstructCommunicationPattern_P2p](@ref PEXSI::PMatrix::ConstructCommunicationPattern_P2p) is called by default.
  virtual void ConstructCommunicationPattern( );


  /// @brief ConstructCommunicationPattern_P2p constructs the communication
  /// pattern to be used later in the selected inversion stage.
  /// The supernodal elimination tree is used to add an additional level of parallelism between supernodes.
  void ConstructCommunicationPattern_P2p( );



  /// @brief PreSelInv prepares the structure in L_ and U_ so that
  /// SelInv only involves matrix-matrix multiplication.
  ///
  /// @todo Move documentation to a more proper place and update the
  /// information.
  ///
  /// Procedure
  /// ---------
  /// PreSelInv performs
  ///
  /// - Compute the inverse of the diagonal blocks
  ///
  ///   L_{kk} <- (L_{kk} U_{kk})^{-1}
  ///
  /// - Update the lower triangular L blocks
  ///
  ///   L_{ik} <- L_{ik} L_{kk}^{-1}
  ///
  /// - Update the upper triangular U blocks which saves redundant
  /// information as in L
  ///
  ///   U_{kj} <- L_{ik}
  ///
  /// Note
  /// ----
  ///
  /// PreSelInv assumes that
  /// PEXSI::PMatrix::ConstructCommunicationPattern has been executed.
  virtual void PreSelInv( );

  /// @brief SelInv is the main function for the selected inversion.
  ///
  /// @todo Move documentation to a more proper place and update the
  /// information.
  ///
  /// Procedure
  /// ---------
  ///
  /// PSelInv is a right-looking based parallel selected inversion
  /// subroutine for sparse symmetric matrices.  Static pivoting is
  /// used in this version.
  ///
  /// Although the matrix is symmetric, the key idea of the current
  /// implementation of PSelInv is that the upper-triangular matrix is
  /// saved (in the form of UBlock).  Such redundant information is
  /// effective for reducing the complexity for designing the
  /// communication pattern.
  ///
  /// At each supernode ksup, the lower triangular part Ainv(isup, ksup)
  /// (isup > ksup) are first updated.  The blocks in the processor row
  /// of ksup first sends the nonzero blocks in U(ksup, jsup) (which is
  /// L(isup, ksup)^T) to the Schur complements Ainv(isup, jsup).  At
  /// the same time the blocks in the processor column of ksup sends the
  /// nonzero blocks (only nonzero row indices) to the Schur complement
  /// Ainv(isup, jsup).  Then
  ///
  /// sum_{jsup} Ainv(isup, jsup) U^{T}(ksup, jsup)
  ///
  /// is performed.  In this procedure, only processors with
  /// isRecvFromAbove[ksup] == true && isRecvFromLeft[ksup] == true
  /// participate in the computation.
  ///
  ///
  /// The result is reduced to the processor column ksup within the same
  /// processor row.  The diagonal block Ainv(ksup, ksup) is simply updated
  /// by a reduce procedure within the column processor group of ksup.
  ///
  /// Then we update the Ainv(ksup, isup) blocks, simply via the update
  /// from the cross diagonal processors.
  ///
  /// <b> NOTE </b>: The cross diagonal processor is only well here
  /// defined for square grids.  For a P x P square grid, (ip,
  /// jp) is the cross diagonal processor of (jp, ip) if ip != jp.  The
  /// current version of SelInv only works for square processor grids.
  ///
  ///
  /// Communication pattern
  /// ---------------------
  ///
  /// The communication is controlled by 3 sending varaibles and 3
  /// receiving variables. The first dimension of all the sending and
  /// receiving variables are numSuper.  The information contains
  /// redundancy since not all processors have access to all the
  /// supernodes.  However, this increases the readability of the output
  /// significantly and only increases a small amount of memory cost for
  /// indexing.  This set of sending / receiving mechanism avoids the
  /// double indexing of the supernodes and can scale to matrices of
  /// large size.
  ///
  /// - isSendToBelow:
  ///
  ///   Dimension: numSuper x numProcRow
  ///
  ///   Role     : At supernode ksup, if isSendToBelow(ksup, ip) == true, send
  ///   all local blocks {U(ksup, jsup) | jsup > ksup} to the processor row ip.
  ///
  /// - isRecvFromAbove:
  ///
  ///   Dimension: numSuper
  ///
  ///   Role     :
  ///
  ///     * At supernode ksup, if isRecvFromAbove(ksup) == true,
  ///       receive blocks from the processor owning the block row of ksup
  ///       within the same column processor group.
  ///
  ///     * If isRecvFromAbove(ksup) == true && isRecvFromLeft(ksup) ==
  ///     true, the ucrrent processor participate in updating Ainv(isup,
  ///     ksup).
  ///
  ///
  /// - isSendToRight:
  ///
  ///   Dimension: numSuper x numProcCol
  ///
  ///   Role     : At supernode ksup, if isSendToRight(ksup, jp) == true, send
  ///   all local blocks (mainly the nonzero row indicies, without the
  ///   values to save the communication cost) {L(isup, ksup) | isup >
  ///   ksup} to the processor column jp.
  ///
  /// - isRecvFromLeft:
  ///
  ///   Dimension: numSuper
  ///
  ///   Role     :
  ///
  ///     * At supernode ksup, if isRecvFromLeft(ksup) == true, receive
  ///     blocks from the processor owning the block column of ksup
  ///     within the same row processor group.
  ///
  ///     * If isRecvFromAbove(ksup) == true && isRecvFromLeft(ksup) ==
  ///     true, the ucrrent processor participate in updating Ainv(isup,
  ///     ksup).
  ///
  /// - isSendToCrossDiagonal:
  ///
  ///   Dimension: numSuper
  ///
  ///   Role     : At supernode ksup, if isSendToCrossDiagonal(ksup) ==
  ///   true, send all local blocks {(isup, ksup) | isup > ksup} to the
  ///   cross-diagonal processor.  <b> NOTE </b>: This requires a square
  ///   processor grid.
  ///
  /// - isRecvCrossDiagonal:
  ///
  ///   Dimension: numSuper
  ///
  ///   Role     : At supernode ksup, if isRecvFromCrossDiagonal(ksup) ==
  ///   true, receive from the cross-diagonal processor.  <b> NOTE </b>:
  ///   This requires a square processor grid.
  ///
  ///
  ///
  ///
  virtual void SelInv( );

  /// @brief Point-to-point version of the selected inversion.
  void SelInv_P2p( );


  /// @brief GetDiagonal extracts the diagonal elements of the PMatrix.
  ///
  /// 1) diag is permuted back to the natural order
  ///
  /// 2) diag is shared by all processors in grid_->comm through a
  /// Allreduce procedure.
  void GetDiagonal( NumVec<T>& diag );
  void GetColumn	( Int colIdx,  NumVec<T>& col );


  /// @brief PMatrixToDistSparseMatrix converts the PMatrix into a
  /// distributed compressed sparse column matrix format.
  /// The DistSparseMatrix follows the natural order.
  ///
  /// @param[out] A Output sparse matrix.
  void PMatrixToDistSparseMatrix( DistSparseMatrix<T>& A );

  /// @brief PMatrixToDistSparseMatrix_OLD converts the PMatrix into a
  /// distributed compressed sparse column matrix format B, which has
  /// the same sparsity pattern as the DistSparseMatrix A.
  ///
  /// The DistSparseMatrix follows the natural order.
  ///
  /// @param[in]  A Input sparse matrix to provide the sparsity pattern.
  ///
  /// @param[out] B Output sparse matrix.
  void PMatrixToDistSparseMatrix_OLD(
      const DistSparseMatrix<T>& A,
      DistSparseMatrix<T>& B	);


  /// @brief PMatrixToDistSparseMatrix converts the PMatrix into a
  /// distributed compressed sparse column matrix format B, which has
  /// the same sparsity pattern as the DistSparseMatrix A.
  ///
  /// The DistSparseMatrix follows the natural order.
  ///
  /// @param[in]  A Input sparse matrix to provide the sparsity pattern.
  ///
  /// @param[out] B Output sparse matrix.
  template<typename T1 = T>
  void PMatrixToDistSparseMatrix(
      const DistSparseMatrix<T1>& A, DistSparseMatrix<T>& B );



  /// @brief NnzLocal computes the number of nonzero elements (L and U)
  /// saved locally.
  Int  NnzLocal();

  /// @brief Nnz computes the total number of nonzero elements in the
  /// PMatrix.
  LongInt  Nnz();

  /// @brief GetNegativeInertia computes the negative inertia of a
  /// PMatrix.  This can be used to estimate e.g. the number of
  /// eigenvalues of a matrix below a certain threshold.
  void GetNegativeInertia	( Real& inertia );

  //      inline int IdxToTag(Int lidx, Int tag) { return SELINV_TAG_COUNT*(lidx)+(tag);}

  void DumpLU();

  void CopyLU( const PMatrix & C);
  inline int IdxToTag(Int lidx, Int tag) { return SELINV_TAG_COUNT*(lidx)+(tag);}



};

template<typename T>  class PMatrixUnsym;


} // namespace PEXSI


#include "pexsi/pselinv_impl.hpp"


#endif //_PEXSI_PSELINV_HPP_
