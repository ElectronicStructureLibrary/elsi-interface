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
/// @file pselinv_unsym.hpp
/// @brief Main file for parallel selected inversion on unsymmetric matrices.
/// @date 2013-08-05
#ifndef _PEXSI_PSELINV_UNSYM_HPP_
#define _PEXSI_PSELINV_UNSYM_HPP_

// *********************************************************************
//  Common utilities
// *********************************************************************

#include  "pexsi/environment.hpp"
#include	"pexsi/NumVec.hpp"
#include	"pexsi/NumMat.hpp"
#include  "pexsi/sparse_matrix.hpp"

#include  "pexsi/superlu_dist_interf.hpp"
#include  "pexsi/mpi_interf.hpp"
#include	"pexsi/utility.hpp"
#include	"pexsi/blas.hpp"
#include	"pexsi/lapack.hpp"
#include	"pexsi/pselinv.hpp"


#include <memory>
#include <set>



namespace PEXSI{

/// @brief CDROW returns cross diagonal processor row
inline Int CDROW(Int pnum, const GridType* g )
{ return (pnum % g->numProcCol) % g->numProcRow; }

/// @brief CDCOL returns cross diagonal processor column
inline Int CDCOL(Int pnum, const GridType* g )
{ return (pnum / g->numProcCol) % g->numProcCol; }

/// @brief CDPNUM returns the cross diagonal processor rank.
inline Int CDPNUM(Int pr, Int pc, const GridType* g )
{ return PNUM(pc,pr,g); }


/// @brief MYCDROW returns cross diagonal processor row
inline Int MYCDROW( const GridType* g )
{ return MYCOL(g) % g->numProcRow; }

/// @brief MYCDCOL returns cross diagonal processor column
inline Int MYCDCOL( const GridType* g )
{ return MYROW(g) % g->numProcCol; }

/// @brief MYCDPROC returns the cross diagonal processor rank.
inline Int MYCDPROC( const GridType* g )
{ return PNUM(MYCDROW(g),MYCDCOL(g),g); }




struct CDBuffers;

template<typename T>
  bool LBlockEqualComparator(const LBlock<T> & a,const LBlock<T> & b){
    return a.blockIdx==b.blockIdx;
  }
template<typename T>
  bool LBlockComparator(const LBlock<T> & a,const LBlock<T> & b){
    return a.blockIdx<b.blockIdx;
  }

template<typename T>
  bool UBlockComparator(const UBlock<T> & a,const UBlock<T> & b){
    return a.blockIdx<b.blockIdx;
  }



/**********************************************************************
 * Main data structure in PSelInv: PMatrixUnsym
 **********************************************************************/

/// @class PMatrixUnsym
///
/// @brief PMatrixUnsym contains the main data structure and the
/// computational routine for the parallel selected inversion.
///
/// **NOTE**
///
///  For general non-symmetric matrices, the selected elements are
///  :math:`\{B_{i,j}\vert A_{j,i}\ne 0\}`.  This means that the matrix
///  elements computed corresponding to the sparsity pattern of
///  :math:`A^T` (for more detailed explanation of the mathematical
///  reason, please see the paper :ref:`PC2018 <citeNonSymPSelInv>`
///  in the PEXSI documentation).  However, storing the matrix elements
///  :math:`\{A^{-1}_{i,j}\vert A_{j,i}\ne 0\}` is practically
///  cumbersome, especially in the context of distributed computing.
///  Hence we choose to store the selected elements for :math:`A^{-T}`,
///  i.e. :math:`\{A^{-T}_{i,j}\vert A_{i,j}\ne 0\}`.  These are the
///  values obtained from the non-symmetric version of PSelInv.


template<typename T>
  class PMatrixUnsym: public PMatrix<T>{

  protected:
    virtual void PMatrixToDistSparseMatrix_( const NumVec<Int> & AcolptrLocal, const NumVec<Int> & ArowindLocal, const Int Asize, const LongInt Annz, const Int AnnzLocal, DistSparseMatrix<T>& B );
    // NOTE: By default the returned matrix should store the transpose of the inverse matrix.

  std::vector<std::vector<Int> > isSendToCD_;
  std::vector<std::vector<Int> > isRecvFromCD_;

  std::vector<std::shared_ptr<TreeBcast_v2<char> > > bcastLStructTree_;

  std::vector<std::vector< std::shared_ptr<TreeBcast_v2<char> > > > bcastUDataTrees_;
  std::vector<std::shared_ptr<TreeBcast_v2<char> > > bcastUDataTree_;
  std::vector<std::shared_ptr<TreeBcast_v2<char> > > bcastUStructTree_;

  std::vector<std::shared_ptr<TreeReduce_v2<T> > > redUTree2_;

  std::vector<TreeReduce<T> *> redLTree_;
  std::vector<TreeReduce<T> *> redDTree_;
  std::vector<TreeReduce<T> *> redUTree_;



    // *********************************************************************
    // Variables
    // *********************************************************************
    // Data variables

    std::vector<std::vector<UBlock<T> > > Ucol_;
    std::vector<std::vector<LBlock<T> > > Lrow_;

    std::vector<Int > UcolSize_;
    std::vector<Int > LrowSize_;

    // Communication variables
    // This is the tag used for mpi communication for selinv

    enum{
      SELINV_TAG_U_SIZE,
      SELINV_TAG_U_CONTENT,
      SELINV_TAG_L_SIZE,
      SELINV_TAG_L_CONTENT,
      SELINV_TAG_UCOL_SIZE,
      SELINV_TAG_UCOL_CONTENT,
      SELINV_TAG_LROW_SIZE,
      SELINV_TAG_LROW_CONTENT,
      SELINV_TAG_L_REDUCE,
      SELINV_TAG_U_REDUCE,
      SELINV_TAG_D_SIZE,
      SELINV_TAG_D_CONTENT,
      SELINV_TAG_D_REDUCE,
      SELINV_TAG_U_SIZE_CD,
      SELINV_TAG_U_CONTENT_CD,
      SELINV_TAG_L_SIZE_CD,
      SELINV_TAG_L_CONTENT_CD,
      SELINV_TAG_COUNT
    };


    struct SuperNodeBufferTypeUnsym:public PMatrix<T>::SuperNodeBufferType {
      NumMat<T>    UUpdateBuf;
      std::vector<Int>  ColLocalPtr;
      std::vector<Int>  BlockIdxLocalU;

      //TODO this will become unecessary
      std::vector<char> SstrLrowSend;
      std::vector<char> SstrUcolSend;
      std::vector<char> SstrLrowRecv;
      std::vector<char> SstrUcolRecv;
      Int               SizeSstrLrowSend;
      Int               SizeSstrUcolSend;
      Int               SizeSstrLrowRecv;
      Int               SizeSstrUcolRecv;


      Int isReadyL;
      Int isReadyU;
      bool updatesLU;

      SuperNodeBufferTypeUnsym():PMatrix<T>::SuperNodeBufferType(){
       isReadyL = 0;
       isReadyU = 0;
       updatesLU = false;
      }

    };

    /// @brief SelInvIntra_P2p
    inline void SelInvIntra_P2p(Int lidx);
    inline void SelInvIntra_New(Int lidx, Int & rank);

    /// @brief SelInv_lookup_indexes
    inline void SelInv_lookup_indexes(SuperNodeBufferTypeUnsym & snode,
        std::vector<LBlock<T> > & LcolRecv,
        std::vector<LBlock<T> > & LrowRecv,
        std::vector<UBlock<T> > & UcolRecv,
        std::vector<UBlock<T> > & UrowRecv,
        NumMat<T> & AinvBuf,
        NumMat<T> & LBuf,
        NumMat<T> & UBuf);
    inline void SelInv_lookup_indexes_New(SuperNodeBufferTypeUnsym & snode,
        std::vector<LBlock<T> > & LRecv,
        std::vector<UBlock<T> > & URecv,
        NumMat<T> & AinvBuf,
        NumMat<T> & LBuf,
        NumMat<T> & UBuf);

    /// @brief UnpackData
    inline void UnpackData( SuperNodeBufferTypeUnsym & snode,
        std::vector<LBlock<T> > & LcolRecv,
        std::vector<LBlock<T> > & LrowRecv,
        std::vector<UBlock<T> > & UcolRecv,
        std::vector<UBlock<T> > & UrowRecv
        );
    inline void UnpackData_New( SuperNodeBufferTypeUnsym & snode,
        std::vector<LBlock<T> > & LRecv,
        std::vector<UBlock<T> > & URecv
        );


    /// @brief ComputeDiagUpdate
    inline void ComputeDiagUpdate(SuperNodeBufferTypeUnsym & snode);
    inline void ComputeDiagUpdate_New(SuperNodeBufferTypeUnsym & snode,bool fromU);

    /// @brief SendRecvCD_UpdateU
    inline void SendRecvCD(
        std::vector<SuperNodeBufferTypeUnsym > & arrSuperNodes,
        Int stepSuper
        );


    inline void SendRecvSizesCD(std::vector<Int > & arrSuperNodes, Int stepSuper, CDBuffers & buffers);
    inline void IRecvContentCD( std::vector<Int > & arrSuperNodes, Int stepSuper, CDBuffers & buffers);
    inline void WaitContentLCD( std::vector<Int > & arrSuperNodes, Int stepSuper, CDBuffers & buffers);
    inline void WaitContentUCD( std::vector<Int > & arrSuperNodes, Int stepSuper, CDBuffers & buffers);






  public:
    // *********************************************************************
    // Public member functions
    // *********************************************************************

    PMatrixUnsym():PMatrix<T>() {}

    PMatrixUnsym( const GridType* g, const SuperNodeType* s, const PSelInvOptions * o, const FactorizationOptions * oFact  );

    //virtual ~PMatrixUnsym() { statusOFS<<"DESTRUCTOR UNSYM CALLED"<<std::endl;    }

    void Setup( const GridType* g, const SuperNodeType* s, const PSelInvOptions * o, const FactorizationOptions * oFact  );

    /// @brief NumBlockL returns the number of nonzero L blocks for the
    /// local block column jLocal.
    //Int NumBlockL( Int jLocal ) const { return L_[jLocal].size(); }

    /// @brief NumBlockU returns the number of nonzero U blocks for the
    /// local block row iLocal.
    //Int NumBlockU( Int iLocal ) const { return U_[iLocal].size(); }


    /// @brief Lrow returns the vector of nonzero L blocks for the local
    /// block row iLocal.
    std::vector<LBlock<T> >& Lrow( Int iLocal ) { return Lrow_[iLocal]; }

    /// @brief Ucol returns the vector of nonzero U blocks for the local
    /// block col jLocal.
    std::vector<UBlock<T> >& Ucol( Int jLocal ) { return Ucol_[jLocal]; }


    /// @brief ConstructCommunicationPattern constructs the communication
    /// pattern to be used later in the selected inversion stage.
    /// The supernodal elimination tree is used to add an additional level of parallelism between supernodes.
    /// [ConstructCommunicationPattern_P2p](@ref PEXSI::PMatrix::ConstructCommunicationPattern_P2p) is called by default.
    virtual void ConstructCommunicationPattern( );


    /// @brief ConstructCommunicationPattern_P2p constructs the communication
    /// pattern to be used later in the selected inversion stage.
    /// The supernodal elimination tree is used to add an additional level of parallelism between supernodes.
    void ConstructCommunicationPattern_P2p( );
    void ConstructCommunicationPattern_New( );


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
    virtual void PreSelInv_New( );

    /// @brief SelInv is the main function for the selected inversion.
    ///
    /// @todo Move documentation to a more proper place and update the
    /// information.
    ///
    /// Procedure
    /// ---------
    ///
    /// PSelInv is a right-looking based parallel selected inversion
    /// subroutine for sparse matrices.  Static pivoting is
    /// used in this version.
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
    virtual void SelInv( );

    /// @brief Point-to-point version of the selected inversion.
    void SelInv_P2p( );

  };



} // namespace PEXSI


#include "pexsi/pselinv_unsym_impl.hpp"

#endif //_PEXSI_PSELINV_UNSYM_HPP_
