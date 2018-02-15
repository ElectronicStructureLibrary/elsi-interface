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
/// @file pselinv_unsym_impl.hpp
/// @brief Implementation of the parallel SelInv.
/// @date 2013-08-05
#ifndef _PEXSI_PSELINV_UNSYM_IMPL_HPP_
#define _PEXSI_PSELINV_UNSYM_IMPL_HPP_

#include "pexsi/timer.h"
#include "pexsi/superlu_dist_interf.hpp"
#include <algorithm>




namespace PEXSI{
template<typename T>
  PMatrixUnsym<T>::PMatrixUnsym ( 
      const GridType* g, 
      const SuperNodeType* s, 
      const PSelInvOptions * o, 
      const FactorizationOptions * oFact
      )
  {

    this->Setup( g, s, o, oFact );

    return ;
  } 		// -----  end of method PMatrixUnsym::PMatrixUnsym  ----- 


template<typename T>
  void PMatrixUnsym<T>::Setup( 
      const GridType* g, 
      const SuperNodeType* s, 
      const PSelInvOptions * o, 
      const FactorizationOptions * oFact  
      ) {

    PMatrix<T>::Setup(g,s,o,oFact);



    Lrow_.clear();
    Ucol_.clear();
    Lrow_.resize( this->NumLocalBlockRow() );
    Ucol_.resize( this->NumLocalBlockCol() );

    LrowSize_.clear();
    UcolSize_.clear();
    LrowSize_.resize( this->NumLocalBlockRow() );
    UcolSize_.resize( this->NumLocalBlockCol() );

    return ;
  } 		// -----  end of method PMatrixUnsym::Setup   ----- 


///////////// Utility functions ///////////////////
template<typename T>
  inline  void PMatrixUnsym<T>::SelInv_lookup_indexes(
      SuperNodeBufferTypeUnsym & snode, 
      std::vector<LBlock<T> > & LcolRecv, 
      std::vector<LBlock<T> > & LrowRecv, 
      std::vector<UBlock<T> > & UcolRecv, 
      std::vector<UBlock<T> > & UrowRecv, /*useless so far*/ 
      NumMat<T> & AinvBuf,
      NumMat<T> & LBuf,
      NumMat<T> & UBuf )
  {
    TIMER_START(Compute_Sinv_LT_Lookup_Indexes);





    TIMER_START(Build_colptr_rowptr);

#if ( _DEBUGlevel_ >= 2 )
    statusOFS<<"UrowRecv blockIdx: "<<std::endl;
    for( Int jb = 0; jb < UrowRecv.size(); jb++ ){
      statusOFS<<UrowRecv[jb].blockIdx<<" ";
    }
    statusOFS<<std::endl;

    statusOFS<<"UcolRecv blockIdx: "<<std::endl;
    for( Int jb = 0; jb < UcolRecv.size(); jb++ ){
      statusOFS<<UcolRecv[jb].blockIdx<<" ";
    }
    statusOFS<<std::endl;

    statusOFS<<"LcolRecv blockIdx: "<<std::endl;
    for( Int jb = 0; jb < LcolRecv.size(); jb++ ){
      statusOFS<<LcolRecv[jb].blockIdx<<" ";
    }
    statusOFS<<std::endl;


    statusOFS<<"LrowRecv blockIdx: "<<std::endl;
    for( Int jb = 0; jb < LrowRecv.size(); jb++ ){
      statusOFS<<LrowRecv[jb].blockIdx<<" ";
    }
    statusOFS<<std::endl;
#endif

    TIMER_START(Sort_LrowRecv_UcolRecv);
    //LrowRecv should be sorted by blockIdx, same for UcolRecv
    std::sort(LrowRecv.begin(),LrowRecv.end(),LBlockComparator<T>);
    std::sort(UcolRecv.begin(),UcolRecv.end(),UBlockComparator<T>);
    TIMER_STOP(Sort_LrowRecv_UcolRecv);

    // rowPtrL[ib] gives the row index in snode.LUpdateBuf for the first
    // nonzero row in LcolRecv[ib]. The total number of rows in
    // snode.LUpdateBuf is given by rowPtr[end]-1
    std::vector<Int> rowPtrL(LcolRecv.size() + 1);
    rowPtrL[0] = 0;
    for( Int ib = 0; ib < LcolRecv.size(); ib++ ){
      LBlock<T> & LB = LcolRecv[ib];
      rowPtrL[ib+1] = rowPtrL[ib] + LB.numRow;
    }

    // colPtrL[jb] gives the column index in LBuf for the first
    // nonzero column in UrowRecv[jb]. The total number of rows in
    // LBuf is given by colPtr[end]-1

    std::vector<Int> colPtrL(LrowRecv.size() + 1);
    std::vector<Int> colPtrU(UcolRecv.size() + 1,-1);
    colPtrL[0] = 0;
    for( Int jb = 0; jb < LrowRecv.size(); jb++ ){
      LBlock<T> & LB = LrowRecv[jb];
      colPtrL[jb+1] = colPtrL[jb] + LB.numCol;
    }

    colPtrU[0] = 0;
    for( Int jb = 0; jb < UcolRecv.size(); jb++ ){
      UBlock<T> & UB = UcolRecv[jb];
      colPtrU[jb+1] = colPtrU[jb] + UB.numCol;
    }

#if ( _DEBUGlevel_ >= 1 )
    statusOFS<<rowPtrL<<std::endl;
    statusOFS<<colPtrL<<std::endl;
    statusOFS<<colPtrU<<std::endl;
#endif


    Int numRowAinvBuf = *colPtrU.rbegin();
    //Int numRowAinvBuf = *rowPtrL.rbegin();
    Int numColAinvBuf = *colPtrL.rbegin();
    TIMER_STOP(Build_colptr_rowptr);

    TIMER_START(Allocate_lookup);
    // Allocate for the computational storage
    AinvBuf.Resize( numRowAinvBuf, numColAinvBuf );
    LBuf.Resize( SuperSize( snode.Index, this->super_ ), numColAinvBuf );
    UBuf.Resize( SuperSize( snode.Index, this->super_ ), numRowAinvBuf );
    TIMER_STOP(Allocate_lookup);



    TIMER_START(Fill_LBuf);
    // Fill LBuf first. Make the transpose later in the Gemm phase.
    for( Int jb = 0; jb < LrowRecv.size(); jb++ ){
      LBlock<T>& LB = LrowRecv[jb];
      if(1 || (LB.numRow>0 && LB.numCol>0)){
        if( LB.numRow != SuperSize(snode.Index, this->super_) ){
          ErrorHandling( 
              "The size of LB is not right. Something is seriously wrong." );
        }
        lapack::Lacpy( 'A', LB.numRow, LB.numCol, LB.nzval.Data(),
            LB.numRow, LBuf.VecData( colPtrL[jb] ), LBuf.m() );
      }
    }
    TIMER_STOP(Fill_LBuf);

    TIMER_START(Fill_UBuf);
    // Fill UBuf first. 
    SetValue(UBuf, ZERO<T>());
    for( Int jb = 0; jb < UcolRecv.size(); jb++ ){
      UBlock<T>& UB = UcolRecv[jb];

      if(1||(UB.numRow>0 && UB.numCol>0)){
        if( UB.numRow != SuperSize(snode.Index, this->super_) ){
          ErrorHandling( 
              "The size of UB is not right. Something is seriously wrong." );
        }


        lapack::Lacpy( 'A', UB.numRow, UB.numCol, UB.nzval.Data(),
            UB.numRow, UBuf.VecData( colPtrU[jb] ), UBuf.m() );
      }
    }
    TIMER_STOP(Fill_UBuf);

    // Calculate the relative indices for (isup, jsup)
    // Fill AinvBuf with the information in L or U block.
    TIMER_START(JB_Loop);
    //    for( Int jb = 0; jb < UrowRecv.size(); jb++ ){
    for( Int jb = 0; jb < LrowRecv.size(); jb++ ){
      LBlock<T>& LrowB = LrowRecv[jb];
      Int jsup = LrowB.blockIdx;

      Int SinvColsSta = FirstBlockCol( jsup, this->super_ );


      for( Int ib = 0; ib < LcolRecv.size(); ib++ ){
        LBlock<T>& LB = LcolRecv[ib];
        Int isup = LB.blockIdx;
        Int SinvRowsSta = FirstBlockCol( isup, this->super_ );
        T* nzvalAinv = &AinvBuf( rowPtrL[ib], colPtrL[jb] );
        Int     ldAinv    = numRowAinvBuf;

        // Pin down the corresponding block in the part of Sinv.
        if( isup >= jsup ){
          // Column relative indicies
          std::vector<Int> relCols( LrowB.numCol );
          for( Int j = 0; j < LrowB.numCol; j++ ){
            relCols[j] = LrowB.rows[j] - SinvColsSta;
          }

          std::vector<LBlock<T> >& LcolSinv = this->L( LBj(jsup, this->grid_ ));
          bool isBlockFound = false;
          TIMER_START(PARSING_ROW_BLOCKIDX);
          for( Int ibSinv = 0; ibSinv < LcolSinv.size(); ibSinv++ ){
            // Found the (isup, jsup) block in Sinv
            if( LcolSinv[ibSinv].blockIdx == isup ){
              LBlock<T>& SinvB = LcolSinv[ibSinv];

              // Row relative indices
              std::vector<Int> relRows( LB.numRow );
              Int* rowsLBPtr    = LB.rows.Data();
              Int* rowsSinvBPtr = SinvB.rows.Data();

              TIMER_START(STDFIND_ROW);
              Int * pos =&rowsSinvBPtr[0];
              Int * last =&rowsSinvBPtr[SinvB.numRow];
              for( Int i = 0; i < LB.numRow; i++ ){
                pos = std::upper_bound(pos, last, rowsLBPtr[i])-1;
                if(pos != last){
                  relRows[i] =  (Int)(pos - rowsSinvBPtr);
                }
                else{
                  std::ostringstream msg;
                  msg << "Row " << rowsLBPtr[i] <<
                    " in LB cannot find the corresponding row in SinvB" 
                    << std::endl
                    << "LB.rows    = " << LB.rows << std::endl
                    << "SinvB.rows = " << SinvB.rows << std::endl;
                  ErrorHandling( msg.str().c_str() );
                }
              }
              TIMER_STOP(STDFIND_ROW);

              TIMER_START(Copy_Sinv_to_Ainv);
              // Transfer the values from Sinv to AinvBlock
              T* nzvalSinv = SinvB.nzval.Data();
              Int     ldSinv    = SinvB.numRow;
              for( Int j = 0; j < LrowB.numCol; j++ ){
                for( Int i = 0; i < LB.numRow; i++ ){
                  //                    AinvBuf( rowPtrL[ib] + i, colPtrL[jb] + j ) =
                  //                                  SinvB.nzval(relRows[i],relCols[j]);
                  nzvalAinv[i + j*ldAinv] =
                    nzvalSinv[relRows[i]+relCols[j]*ldSinv];
                }
              }
              TIMER_STOP(Copy_Sinv_to_Ainv);

              isBlockFound = true;
              break;
            }
          } // for (ibSinv )
          TIMER_STOP(PARSING_ROW_BLOCKIDX);
          if( isBlockFound == false ){
            std::ostringstream msg;
            msg << "["<<snode.Index<<"] "<<"Block(" << isup << ", " << jsup
              << ") did not find a matching block in Sinv." << std::endl;
            //#if ( _DEBUGlevel_ >= 1 )
            statusOFS<<msg.str();
            //#endif
            //            ErrorHandling( msg.str().c_str() );
          }
        } // if (isup, jsup) is in L
        else{
          // Row relative indices
          std::vector<Int> relRows( LB.numRow );
          for( Int i = 0; i < LB.numRow; i++ ){
            relRows[i] = LB.rows[i] - SinvRowsSta;
          }

          std::vector<UBlock<T> >& UrowSinv = 
            this->U( LBi( isup, this->grid_ ) );
          bool isBlockFound = false;
          TIMER_START(PARSING_COL_BLOCKIDX);
          for( Int jbSinv = 0; jbSinv < UrowSinv.size(); jbSinv++ ){
            // Found the (isup, jsup) block in Sinv
            if( UrowSinv[jbSinv].blockIdx == jsup ){
              UBlock<T>& SinvB = UrowSinv[jbSinv];

              // Column relative indices
              std::vector<Int> relCols( LrowB.numCol );
              Int* rowsLrowBPtr    = LrowB.rows.Data();
              Int* colsSinvBPtr = SinvB.cols.Data();

              TIMER_START(STDFIND_COL);
              Int * pos =&colsSinvBPtr[0];
              Int * last =&colsSinvBPtr[SinvB.numCol];
              for( Int j = 0; j < LrowB.numCol; j++ ){
                //rowsLrowB is sorted
                pos = std::upper_bound(pos, last, rowsLrowBPtr[j])-1;
                if(pos !=last){
                  relCols[j] = (Int)(pos - colsSinvBPtr);
                }
                else{
                  std::ostringstream msg;
                  msg << "Col " << rowsLrowBPtr[j] <<
                    " in UB cannot find the corresponding row in SinvB"
                    << std::endl
                    << "LrowB.rows    = " << LrowB.rows << std::endl
                    << "UinvB.cols = " << SinvB.cols << std::endl;
                  ErrorHandling( msg.str().c_str() );
                }
              }
              TIMER_STOP(STDFIND_COL);

              TIMER_START(Copy_Sinv_to_Ainv);
              // Transfer the values from Sinv to AinvBlock
              T* nzvalSinv = SinvB.nzval.Data();
              Int     ldSinv    = SinvB.numRow;
              for( Int j = 0; j < LrowB.numCol; j++ ){
                for( Int i = 0; i < LB.numRow; i++ ){
                  //                  AinvBuf( rowPtrL[ib] + i, colPtrL[jb] + j ) =
                  //                                SinvB.nzval(relRows[i],relCols[j]);

                  nzvalAinv[i + j*ldAinv] =
                    nzvalSinv[relRows[i]+relCols[j]*ldSinv];
                }
              }
              TIMER_STOP(Copy_Sinv_to_Ainv);

              isBlockFound = true;
              break;
            }
          } // for (jbSinv)
          TIMER_STOP(PARSING_COL_BLOCKIDX);
          if( isBlockFound == false ){
            std::ostringstream msg;
            msg << "["<<snode.Index<<"] "<< "Block(" << isup << ", " << jsup
              << ") did not find a matching block in Sinv." << std::endl;
            statusOFS<<msg.str();
            //            ErrorHandling( msg.str().c_str() );
          }
        } // if (isup, jsup) is in U

      } // for( ib )
    } // for ( jb )
    TIMER_STOP(JB_Loop);

    TIMER_STOP(Compute_Sinv_LT_Lookup_Indexes);
  } // End of method PMatrixUnsym::Selinv_lookup_indexes



  struct CDBuffers{
    //Buffers for L
    std::vector<MPI_Request > arrMpiReqsSendLCD;
    std::vector<MPI_Request > arrMpiReqsSizeSendLCD;
    std::vector<MPI_Request > arrMpiReqsRecvLCD;
    std::vector<MPI_Request > arrMpiReqsSizeRecvLCD;
    std::vector<std::vector<char> > arrSstrLcolSendCD;
    std::vector<int > arrSstrLcolSizeSendCD;
    std::vector<std::vector<char> > arrSstrLrowRecvCD;
    std::vector<int > arrSstrLrowSizeRecvCD;

    //Buffers for U
    std::vector<MPI_Request > arrMpiReqsSendUCD;
    std::vector<MPI_Request > arrMpiReqsSizeSendUCD;
    std::vector<MPI_Request > arrMpiReqsRecvUCD;
    std::vector<MPI_Request > arrMpiReqsSizeRecvUCD;
    std::vector<std::vector<char> > arrSstrUrowSendCD;
    std::vector<int > arrSstrUrowSizeSendCD;
    std::vector<std::vector<char> > arrSstrUcolRecvCD;
    std::vector<int > arrSstrUcolSizeRecvCD;

    std::vector<Int>sendOffset;
    std::vector<Int>recvOffset;

    void resize(Int sendCount, Int recvCount){
      //Buffers for L
      arrMpiReqsSendLCD.assign(sendCount, MPI_REQUEST_NULL );
      arrMpiReqsSizeSendLCD.assign(sendCount, MPI_REQUEST_NULL );
      arrMpiReqsRecvLCD.assign(recvCount, MPI_REQUEST_NULL );
      arrMpiReqsSizeRecvLCD.assign(recvCount, MPI_REQUEST_NULL );
      arrSstrLcolSendCD.resize(sendCount);
      arrSstrLcolSizeSendCD.resize(sendCount);
      arrSstrLrowRecvCD.resize(recvCount);
      arrSstrLrowSizeRecvCD.resize(recvCount);

      //Buffers for U
      arrMpiReqsSendUCD.assign(recvCount, MPI_REQUEST_NULL );
      arrMpiReqsSizeSendUCD.assign(recvCount, MPI_REQUEST_NULL );
      arrMpiReqsRecvUCD.assign(sendCount, MPI_REQUEST_NULL );
      arrMpiReqsSizeRecvUCD.assign(sendCount, MPI_REQUEST_NULL );
      arrSstrUrowSendCD.resize(recvCount);
      arrSstrUrowSizeSendCD.resize(recvCount);
      arrSstrUcolRecvCD.resize(sendCount);
      arrSstrUcolSizeRecvCD.resize(sendCount);
    }

    void WaitAllSend(){
      mpi::Waitall(arrMpiReqsSizeSendLCD);
      mpi::Waitall(arrMpiReqsSendLCD);
      mpi::Waitall(arrMpiReqsSizeSendUCD);
      mpi::Waitall(arrMpiReqsSendUCD);

      //        MPI_Waitall(arrMpiReqsSizeSendLCD.size(), &arrMpiReqsSizeSendLCD[0], MPI_STATUSES_IGNORE); 
      //        MPI_Waitall(arrMpiReqsSendLCD.size(), &arrMpiReqsSendLCD[0], MPI_STATUSES_IGNORE); 
      //        MPI_Waitall(arrMpiReqsSizeSendUCD.size(), &arrMpiReqsSizeSendUCD[0], MPI_STATUSES_IGNORE); 
      //        MPI_Waitall(arrMpiReqsSendUCD.size(), &arrMpiReqsSendUCD[0], MPI_STATUSES_IGNORE); 
    }

    void WaitAllRecv(){
      //        mpi::Waitall(arrMpiReqsSizeRecvLCD);
      //        mpi::Waitall(arrMpiReqsRecvLCD);
      //        mpi::Waitall(arrMpiReqsSizeRecvUCD);
      //        mpi::Waitall(arrMpiReqsRecvUCD);
    }


  };


  template<typename T>
    inline void PMatrixUnsym<T>::SendRecvSizesCD(
        std::vector<Int > & arrSuperNodes, 
        Int stepSuper, CDBuffers & buffers) 
    {
      TIMER_START(SendRecvSizesCD);
      //compute the number of requests
      Int sendCount = 0;
      Int recvCount = 0;
      buffers.sendOffset.resize(stepSuper);
      buffers.recvOffset.resize(stepSuper);
      Int recvIdxL=0;
      Int recvIdxU=0;
      for (Int supidx=0; supidx<stepSuper; supidx++){
        Int snode_index = arrSuperNodes[supidx];
        buffers.sendOffset[supidx]=sendCount;
        buffers.recvOffset[supidx]=recvCount;
        sendCount+= this->CountSendToCrossDiagonal(snode_index);
        recvCount+= this->CountRecvFromCrossDiagonal(snode_index);
      }

      buffers.resize(sendCount,recvCount);


      //Do Isend for size and content of L and Irecv for sizes of U 
      for (Int supidx=0; supidx<stepSuper; supidx++){
        Int snode_index = arrSuperNodes[supidx];

        TIMER_START(Send_L_Recv_Size_U_CrossDiag);

        if( MYCOL( this->grid_ ) == PCOL( snode_index, this->grid_ ) 
            && this->isSendToCrossDiagonal_(this->grid_->numProcCol, snode_index ) ){

          Int sendIdxL = 0;
          Int recvIdxU = 0;
          for(Int dstCol = 0; dstCol<this->grid_->numProcCol; dstCol++){
            if(this->isSendToCrossDiagonal_(dstCol,snode_index) ){
              Int dest = PNUM(PROW(snode_index,this->grid_),dstCol,this->grid_);

              if( MYPROC( this->grid_ ) != dest	){
                //Send size and content of L
                MPI_Request & mpiReqSizeSendL = 
                  buffers.arrMpiReqsSizeSendLCD[buffers.sendOffset[supidx]+sendIdxL];
                MPI_Request & mpiReqSendL = 
                  buffers.arrMpiReqsSendLCD[buffers.sendOffset[supidx]+sendIdxL];
                std::vector<char> & sstrLcolSend = 
                  buffers.arrSstrLcolSendCD[buffers.sendOffset[supidx]+sendIdxL];
                Int & sstrSizeL = 
                  buffers.arrSstrLcolSizeSendCD[buffers.sendOffset[supidx]+sendIdxL];

                // Pack L data
                std::stringstream sstm;
                std::vector<Int> mask( LBlockMask::TOTAL_NUMBER, 1 );
                std::vector<LBlock<T> >&  Lcol = this->L( LBj(snode_index, this->grid_) );
                // All blocks except for the diagonal block are to be sent right
                //TODO not true > this is a scatter operation ! Can we know the destination ?

                //Skip the diagonal block if necessary 
                Int startIdx = ( MYROW( this->grid_ ) == PROW( snode_index, this->grid_ ) ) ? 1:0;
                Int count = 0;
                for( Int ib = startIdx; ib < Lcol.size(); ib++ ){
                  if( Lcol[ib].blockIdx > snode_index && 
                      (Lcol[ib].blockIdx % this->grid_->numProcCol) == dstCol  ){
                    count++;
                  }
                }

                serialize( (Int)count, sstm, NO_MASK );

                for( Int ib = startIdx; ib < Lcol.size(); ib++ ){
                  if( Lcol[ib].blockIdx > snode_index &&  
                      (Lcol[ib].blockIdx % this->grid_->numProcCol) == dstCol  ){ 
#if ( _DEBUGlevel_ >= 1 )
                    statusOFS << "["<<snode_index<<"] SEND contains "<<Lcol[ib].blockIdx
                      << " which corresponds to "<<GBj(ib,this->grid_)
                      << " to P"<<dest
                      << std::endl;
#endif
                    serialize( Lcol[ib], sstm, mask );
                  }
                }

                sstrLcolSend.resize( Size(sstm) );
                sstm.read( &sstrLcolSend[0], sstrLcolSend.size() );
                sstrSizeL = sstrLcolSend.size();



                MPI_Isend( &sstrSizeL, 1, MPI_INT, dest, 
                    IDX_TO_TAG2(snode_index,supidx,SELINV_TAG_L_SIZE_CD),
                    this->grid_->comm, &mpiReqSizeSendL );
                MPI_Isend( (void*)&sstrLcolSend[0], sstrSizeL, MPI_BYTE, dest,
                    IDX_TO_TAG2(snode_index,supidx,SELINV_TAG_L_CONTENT_CD),
                    this->grid_->comm, &mpiReqSendL );
                sendIdxL++;

                //Recv for U size
                Int & sstrSizeU = 
                  buffers.arrSstrUcolSizeRecvCD[buffers.sendOffset[supidx]+recvIdxU];
                MPI_Request & mpiReqSizeRecvU = 
                  buffers.arrMpiReqsSizeRecvUCD[buffers.sendOffset[supidx]+recvIdxU];

                MPI_Irecv( &sstrSizeU, 1, MPI_INT, dest, 
                    IDX_TO_TAG2(snode_index,supidx,SELINV_TAG_U_SIZE_CD),
                    this->grid_->comm, &mpiReqSizeRecvU );
                recvIdxU++;
              }
            }
          }
        } // sender
        TIMER_STOP(Send_L_Recv_Size_U_CrossDiag);
      }


      //Do Irecv for sizes of L and Isend for size and content of U
      for (Int supidx=0; supidx<stepSuper; supidx++){
        Int snode_index = arrSuperNodes[supidx];
        TIMER_START(Send_U_Recv_Size_L_CrossDiag);
        //If I'm a receiver
        if( MYROW( this->grid_ ) == PROW( snode_index, this->grid_ ) 
            && this->isRecvFromCrossDiagonal_(this->grid_->numProcRow, snode_index ) ){
          Int recvIdxL=0;
          Int sendIdxU=0;
          for(Int srcRow = 0; srcRow<this->grid_->numProcRow; srcRow++){
            if(this->isRecvFromCrossDiagonal_(srcRow,snode_index) ){
              Int src = PNUM(srcRow,PCOL(snode_index,this->grid_),this->grid_);
              if( MYPROC( this->grid_ ) != src ){
                //Recv size of L
                Int & sstrSizeL = 
                  buffers.arrSstrLrowSizeRecvCD[buffers.recvOffset[supidx]+recvIdxL];
                MPI_Request & mpiReqSizeRecvL = 
                  buffers.arrMpiReqsSizeRecvLCD[buffers.recvOffset[supidx]+recvIdxL];

                MPI_Irecv( &sstrSizeL, 1, MPI_INT, src, 
                    IDX_TO_TAG2(snode_index,supidx,SELINV_TAG_L_SIZE_CD),
                    this->grid_->comm, &mpiReqSizeRecvL );
                recvIdxL++;


                //Send size and content of U
                MPI_Request & mpiReqSizeSendU = 
                  buffers.arrMpiReqsSizeSendUCD[buffers.recvOffset[supidx]+sendIdxU];
                MPI_Request & mpiReqSendU = 
                  buffers.arrMpiReqsSendUCD[buffers.recvOffset[supidx]+sendIdxU];
                std::vector<char> & sstrUrowSend = 
                  buffers.arrSstrUrowSendCD[buffers.recvOffset[supidx]+sendIdxU];
                Int & sstrSizeU = 
                  buffers.arrSstrUrowSizeSendCD[buffers.recvOffset[supidx]+sendIdxU];

                // Pack U data
                std::stringstream sstm;
                std::vector<Int> mask( UBlockMask::TOTAL_NUMBER, 1 );
                std::vector<UBlock<T> >&  Urow = this->U( LBi(snode_index, this->grid_) );
                // All blocks except for the diagonal block are to be sent right
                //TODO not true > this is a scatter operation ! Can we know the destination ?

                Int count = 0;
                for( Int jb = 0; jb < Urow.size(); jb++ ){
                  if( Urow[jb].blockIdx > snode_index 
                      && (Urow[jb].blockIdx % this->grid_->numProcRow) == srcRow ){
                    count++;
                  }
                }

                serialize( (Int)count, sstm, NO_MASK );

                for( Int jb = 0; jb < Urow.size(); jb++ ){
                  if( Urow[jb].blockIdx > snode_index 
                      && (Urow[jb].blockIdx % this->grid_->numProcRow) == srcRow ){ 
#if ( _DEBUGlevel_ >= 1 )
                    statusOFS << "["<<snode_index<<"] SEND contains "<<Urow[jb].blockIdx
                      << " which corresponds to "<<GBi(jb,this->grid_)
                      << " to P"<<src
                      << std::endl;
#endif
                    serialize( Urow[jb], sstm, mask );
                  }
                }

                sstrUrowSend.resize( Size(sstm) );
                sstm.read( &sstrUrowSend[0], sstrUrowSend.size() );
                sstrSizeU = sstrUrowSend.size();


                MPI_Isend( &sstrSizeU, 1, MPI_INT, src, 
                    IDX_TO_TAG2(snode_index,supidx,SELINV_TAG_U_SIZE_CD), 
                    this->grid_->comm, &mpiReqSizeSendU );
                MPI_Isend( (void*)&sstrUrowSend[0], sstrSizeU, MPI_BYTE, src,
                    IDX_TO_TAG2(snode_index,supidx,SELINV_TAG_U_CONTENT_CD),
                    this->grid_->comm, &mpiReqSendU );
                sendIdxU++;
              }
            }
          }
        }//end if I'm a receiver
        TIMER_STOP(Send_U_Recv_Size_L_CrossDiag);
      }
      TIMER_STOP(SendRecvSizesCD);
    }

  template<typename T>
    inline void PMatrixUnsym<T>::IRecvContentCD(
        std::vector<Int > & arrSuperNodes, 
        Int stepSuper, CDBuffers & buffers) {
      TIMER_START(IrecvContentCD);


      //waitall sizes of L
      TIMER_START(Wait_Size_L_CrossDiag);
      mpi::Waitall(buffers.arrMpiReqsSizeRecvLCD);
      TIMER_STOP(Wait_Size_L_CrossDiag);

      //if(this->grid_->mpirank==0){gdb_lock();}

      //Allocate content and do Irecv for content of L
      for (Int supidx=0; supidx<stepSuper; supidx++){
        Int snode_index = arrSuperNodes[supidx];
        //If I'm a receiver
        TIMER_START(Recv_L_CrossDiag);
        if( MYROW( this->grid_ ) == PROW( snode_index, this->grid_ ) 
            && this->isRecvFromCrossDiagonal_(this->grid_->numProcRow, snode_index ) ){
          Int recvIdxL=0;
          for(Int srcRow = 0; srcRow<this->grid_->numProcRow; srcRow++){
            if(this->isRecvFromCrossDiagonal_(srcRow,snode_index) ){
              Int src = PNUM(srcRow,PCOL(snode_index,this->grid_),this->grid_);
              if( MYPROC( this->grid_ ) != src ){
                Int & sstrSizeL = 
                  buffers.arrSstrLrowSizeRecvCD[buffers.recvOffset[supidx]+recvIdxL];
                std::vector<char> & sstrLrowRecv = 
                  buffers.arrSstrLrowRecvCD[buffers.recvOffset[supidx]+recvIdxL];
                MPI_Request & mpiReqRecvL = 
                  buffers.arrMpiReqsRecvLCD[buffers.recvOffset[supidx]+recvIdxL];
                sstrLrowRecv.resize( sstrSizeL);

                MPI_Irecv( (void*)&sstrLrowRecv[0], sstrSizeL, MPI_BYTE, src,
                    IDX_TO_TAG2(snode_index,supidx,SELINV_TAG_L_CONTENT_CD),
                    this->grid_->comm, &mpiReqRecvL );
                recvIdxL++;

              }
            }
          }
        }//end if I'm a receiver
        TIMER_STOP(Recv_L_CrossDiag);
      }

      //waitall sizes of U
      TIMER_START(Wait_Size_U_CrossDiag);
      mpi::Waitall(buffers.arrMpiReqsSizeRecvUCD);
      TIMER_STOP(Wait_Size_U_CrossDiag);

      //Allocate content and do Irecv for content of U
      for (Int supidx=0; supidx<stepSuper; supidx++){
        Int snode_index = arrSuperNodes[supidx];
        TIMER_START(Recv_U_CrossDiag);
        if( MYCOL( this->grid_ ) == PCOL( snode_index, this->grid_ ) 
            && this->isSendToCrossDiagonal_(this->grid_->numProcCol, snode_index ) ){
          Int recvIdxU = 0;
          for(Int srcCol = 0; srcCol<this->grid_->numProcCol; srcCol++){
            if(this->isSendToCrossDiagonal_(srcCol,snode_index) ){
              Int src = PNUM(PROW(snode_index,this->grid_),srcCol,this->grid_);
              if( MYPROC( this->grid_ ) != src ){
                Int & sstrSizeU = 
                  buffers.arrSstrUcolSizeRecvCD[buffers.sendOffset[supidx]+recvIdxU];
                std::vector<char> & sstrUcolRecv = 
                  buffers.arrSstrUcolRecvCD[buffers.sendOffset[supidx]+recvIdxU];
                MPI_Request & mpiReqRecvU = 
                  buffers.arrMpiReqsRecvUCD[buffers.sendOffset[supidx]+recvIdxU];
                sstrUcolRecv.resize( sstrSizeU );

                MPI_Irecv( (void*)&sstrUcolRecv[0], sstrSizeU, MPI_BYTE, src,
                    IDX_TO_TAG2(snode_index,supidx,SELINV_TAG_U_CONTENT_CD),
                    this->grid_->comm, &mpiReqRecvU );
                recvIdxU++;
              }
            }
          }
        }//end if I'm a receiver
        TIMER_STOP(Recv_U_CrossDiag);
      }

      TIMER_START(IrecvContentCD);
    }


  template<typename T>
    inline void PMatrixUnsym<T>::WaitContentLCD(
        std::vector<Int > & arrSuperNodes, 
        Int stepSuper, CDBuffers & buffers) {
      TIMER_START(WaitContentLCD);
      //waitall content of L
      TIMER_START(Wait_L_CrossDiag);
      mpi::Waitall(buffers.arrMpiReqsRecvLCD);
      TIMER_STOP(Wait_L_CrossDiag);



      //Do the work for Lrow
      NumMat<T> Ltmp;
      for (Int supidx=0; supidx<stepSuper; supidx++){
        Int snode_index = arrSuperNodes[supidx];
        Int ksupProcRow = PROW( snode_index, this->grid_ );
        Int ksupProcCol = PCOL( snode_index, this->grid_ );

        if( MYROW( this->grid_ ) == PROW( snode_index, this->grid_ ) &&
            this->isRecvFromCrossDiagonal_(this->grid_->numProcRow, snode_index ) ){

#if ( _DEBUGlevel_ >= 1 )
          statusOFS << std::endl << " ["<<snode_index<<"] "
            << "Update the upper triangular block" 
            << std::endl << std::endl; 
          statusOFS << std::endl << " ["<<snode_index<<"] "
            //            << "blockIdxLocal:" << snode.BlockIdxLocal
            << std::endl << std::endl; 
          statusOFS << std::endl << " ["<<snode_index<<"] "
            //            << "rowLocalPtr:" << snode.RowLocalPtr
            << std::endl << std::endl; 
#endif

          std::vector<LBlock<T> >& Lcol = this->L( LBj(snode_index, this->grid_) );
          std::vector<LBlock<T> >& Lrow = this->Lrow( LBi( snode_index, this->grid_ ) );

          //          if(Lrow.size() == 0){
          //            Int& LrowSize = this->LrowSize_[ LBi( snode_index, this->grid_ ) ];
          //            Lrow.resize(LrowSize);
          //          }

          std::vector<Int> isBlockFound(Lrow.size(),false);
          Int recvIdxL=0;

          for(Int srcRow = 0; srcRow<this->grid_->numProcRow; srcRow++){
            if(this->isRecvFromCrossDiagonal_(srcRow,snode_index) ){
              Int src = PNUM(srcRow,PCOL(snode_index,this->grid_),this->grid_);
              TIMER_START(Recv_L_CrossDiag);

              //Declaring some pointers to avoid copy if data is local
              std::vector<LBlock<T> > * pLcol;
              std::vector<LBlock<T> > LcolRecv;
              if(MYPROC(this->grid_)!= src){
                Int & sstrSizeL = 
                  buffers.arrSstrLrowSizeRecvCD[buffers.recvOffset[supidx]+recvIdxL];
                std::vector<char> & sstrLrowRecv = 
                  buffers.arrSstrLrowRecvCD[buffers.recvOffset[supidx]+recvIdxL];

                std::stringstream sstm;

#if ( _DEBUGlevel_ >= 1 )
                statusOFS << "["<<snode_index<<"] P"<<MYPROC(this->grid_)<<" ("<<MYROW(this->grid_)
                  << ","<<MYCOL(this->grid_)<<") <--- LBj("<<snode_index<<") <--- P"
                  << src <<" ("<<srcRow<<","<<ksupProcCol<<")"
                  << std::endl;
#endif
                sstm.write( &sstrLrowRecv[0], sstrSizeL );

                // Unpack L data.  
                Int numLBlock;
                std::vector<Int> mask( LBlockMask::TOTAL_NUMBER, 1 );
                deserialize( numLBlock, sstm, NO_MASK );
                LcolRecv.resize(numLBlock);
                for( Int ib = 0; ib < numLBlock; ib++ ){
                  deserialize( LcolRecv[ib], sstm, mask );
#if ( _DEBUGlevel_ >= 1 )
                  statusOFS << "["<<snode_index<<"] RECV contains "
                    << LcolRecv[ib].blockIdx<< " which corresponds to "
                    << ib * this->grid_->numProcRow + srcRow;
                  statusOFS << " L is on row "<< srcRow 
                    << " whereas Lrow is on col "
                    << LcolRecv[ib].blockIdx % this->grid_->numProcCol 
                    << std::endl;
#endif
                }
                recvIdxL++;
                pLcol = &LcolRecv;
              } // sender is not the same as receiver
              else{
                // L is obtained locally, just make a copy. Do not include the diagonal block
                pLcol = &Lcol;
              } // sender is the same as receiver


              TIMER_STOP(Recv_L_CrossDiag);

              //We can directly put LcolRecv in Lrow (stil need to transpose)

              //Skip the diagonal block if necessary 
              Int startIdx = ( MYPROC( this->grid_ ) == src && MYPROC(this->grid_) == PNUM(ksupProcRow,ksupProcCol, this->grid_) ) ? 1:0;
              for( Int ib = startIdx; ib < pLcol->size(); ib++ ){
                LBlock<T> & LB = (*pLcol)[ib];
                if( LB.blockIdx <= snode_index ){
                  ErrorHandling( "LcolRecv contains the wrong blocks." );
                }

                //check that LB would be on this proc if it was a UB
                if(PCOL(LB.blockIdx, this->grid_) != MYCOL(this->grid_)){
                  continue;
                }



                //std::vector<LBlock<T> > &  LrowCD = this->Lrow( LBi( ksup, super_ ) );
                for( Int iib = 0; iib < Lrow.size(); iib++ ){
                  LBlock<T> &  LrowB = Lrow[ iib ];

                  //If not initialized, set it to be the one
                  if( LrowB.blockIdx == -1){
                    LrowB.blockIdx = LB.blockIdx;
                  }

                  if( LB.blockIdx == LrowB.blockIdx ){
                    // Note that the order of the column indices of the U
                    // block may not follow the order of the row indices,
                    // overwrite the information in U.
                    //                    LrowB = LB;
                    LrowB.rows = LB.rows;
                    //Store in "row major" format / i.e transpose is in col-major
                    Transpose(LB.nzval, LrowB.nzval);
                    LrowB.numCol = LB.numRow;
                    LrowB.numRow = LB.numCol;

                    //                    SetValue(LrowB.nzval,ZERO<T>());

#if ( _DEBUGlevel_ >= 1 )
                    statusOFS<<"["<<snode_index<<"] USING LB "<<LB.blockIdx<< std::endl;
#endif
                    isBlockFound[iib] = 1;//true;
                    break;
                  } // if( LB.blockIdx == LrowB.blockIdx )
                } // for (iib)
              } // for (ib)
            }
          }

          for( Int ib = 0; ib < Lrow.size(); ib++ ){
            LBlock<T> &  LB = Lrow[ib];
            if( !isBlockFound[ib] ){
              ErrorHandling( 
                  "LBlock cannot find its update. Something is seriously wrong."
                  );
            }
          }
        } // Done copying L
      }

      TIMER_STOP(WaitContentLCD);
    }



  template<typename T>
    inline void PMatrixUnsym<T>::WaitContentUCD(
        std::vector<Int > & arrSuperNodes, 
        Int stepSuper, CDBuffers & buffers) {

      TIMER_START(WaitContentUCD);
      //waitall content of U
      TIMER_START(Wait_U_CrossDiag);
      mpi::Waitall(buffers.arrMpiReqsRecvUCD);
      TIMER_STOP(Wait_U_CrossDiag);


      //Do the work for U
      for (Int supidx=0; supidx<stepSuper; supidx++){
        Int snode_index = arrSuperNodes[supidx];
        Int ksupProcRow = PROW( snode_index, this->grid_ );
        Int ksupProcCol = PCOL( snode_index, this->grid_ );


        if( MYCOL( this->grid_ ) == PCOL( snode_index, this->grid_ ) &&
            this->isSendToCrossDiagonal_(this->grid_->numProcCol, snode_index ) ){

          std::vector<UBlock<T> >& Ucol = this->Ucol( LBj( snode_index, this->grid_ ) );

          //          if(Ucol.size() == 0){
          //            Int& UcolSize = this->UcolSize_[ LBj( snode_index, this->grid_ ) ];
          //            Ucol.resize(UcolSize);
          //          }

          std::vector<Int> isBlockFound(Ucol.size(),false);
          Int recvIdxU=0;

          for(Int srcCol = 0; srcCol<this->grid_->numProcCol; srcCol++){
            if(this->isSendToCrossDiagonal_(srcCol,snode_index) ){
              Int src = PNUM(PROW(snode_index,this->grid_),srcCol,this->grid_);
              TIMER_START(Recv_U_CrossDiag);

              //Declaring some pointers to avoid copy if data is local
              std::vector<UBlock<T> > * pUrow;
              std::vector<UBlock<T> > UrowRecv;
              if(MYPROC(this->grid_)!= src){
                Int & sstrSizeU = 
                  buffers.arrSstrUcolSizeRecvCD[buffers.sendOffset[supidx]+recvIdxU];
                std::vector<char> & sstrUcolRecv = 
                  buffers.arrSstrUcolRecvCD[buffers.sendOffset[supidx]+recvIdxU];


                std::stringstream sstm;

#if ( _DEBUGlevel_ >= 1 )
                statusOFS << "["<<snode_index<<"] P"<<MYPROC(this->grid_)<<" ("<<MYROW(this->grid_)
                  << ","<<MYCOL(this->grid_)<<") <--- LBi("<<snode_index<<") <--- P"
                  << src<<" ("<<ksupProcCol<<","<<srcCol<<")"<<std::endl;
#endif
                sstm.write( &sstrUcolRecv[0], sstrSizeU );

                // Unpack U data.  
                Int numUBlock;
                std::vector<Int> mask( UBlockMask::TOTAL_NUMBER, 1 );
                deserialize( numUBlock, sstm, NO_MASK );
                UrowRecv.resize(numUBlock);
                for( Int jb = 0; jb < numUBlock; jb++ ){
                  deserialize( UrowRecv[jb], sstm, mask );
#if ( _DEBUGlevel_ >= 1 )
                  statusOFS << "["<<snode_index<<"] RECV contains "
                    << UrowRecv[jb].blockIdx<< " which corresponds to "
                    << jb * this->grid_->numProcCol + srcCol;
                  statusOFS << " Urow is on col "<< srcCol 
                    << " whereas U is on row "
                    << UrowRecv[jb].blockIdx % this->grid_->numProcRow 
                    << std::endl;
#endif
                }
                pUrow = &UrowRecv;
                recvIdxU++;
              } // sender is not the same as receiver
              else{
                // U is obtained locally, just make a copy. 
                std::vector<UBlock<T> >& Urow = this->U( LBi( snode_index, this->grid_ ) );
                pUrow = &Urow;
              } // sender is the same as receiver


              TIMER_STOP(Recv_U_CrossDiag);

              //We can directly put UrowRecv in Ucol
              for( Int jb = 0; jb < pUrow->size(); jb++ ){
                UBlock<T> & UB = (*pUrow)[jb];
                if( UB.blockIdx <= snode_index ){
                  ErrorHandling( "UrowRecv contains the wrong blocks." );
                }

                //check that UB would be on this proc if it was a LB
                if(PROW(UB.blockIdx, this->grid_) != MYROW(this->grid_)){
                  continue;
                }


                //std::vector<UBlock<T> > &  UcolCD = this->Ucol( LBi( ksup, super_ ) );
                for( Int jjb = 0; jjb < Ucol.size(); jjb++ ){
                  UBlock<T> &  UcolB = Ucol[jjb];


                  //If not initialized, set it to be the one
                  if( UcolB.blockIdx == -1){
                    UcolB.blockIdx = UB.blockIdx;
                  }


                  if( UB.blockIdx == UcolB.blockIdx ){
                    // Compare size
                    // Note that the order of the column indices of the U
                    // block may not follow the order of the row indices,
                    // overwrite the information in U.
                    UcolB = UB;

                    //SetValue(UcolB.nzval,ZERO<T>());
#if ( _DEBUGlevel_ >= 1 )
                    statusOFS<<"["<<snode_index<<"] USING UB "<<UB.blockIdx<< std::endl;
#endif
                    isBlockFound[jjb] = 1;
                    break;
                  } // if( UB.blockIdx == UcolB.blockIdx )
                } // for (jjb)
              } // for (jb)

            }
          }

          for( Int jb = 0; jb < Ucol.size(); jb++ ){
            UBlock<T> &  UB = Ucol[jb];
            if( !isBlockFound[jb] ){

              ErrorHandling( 
                  "UBlock cannot find its update. Something is seriously wrong."
                  );
            }
          } // end for (jb)
        } // Done copying U

      }

      TIMER_STOP(WaitContentUCD);
    }




  template<typename T>
    inline void PMatrixUnsym<T>::SendRecvCD(
        std::vector<SuperNodeBufferTypeUnsym > & arrSuperNodes, 
        Int stepSuper)
    {

      TIMER_START(Send_CD);
      //compute the number of requests
      Int sendCount = 0;
      Int recvCount = 0;
      Int sendOffset[stepSuper];
      Int recvOffset[stepSuper];
      Int recvIdxL=0;
      Int recvIdxU=0;
      for (Int supidx=0; supidx<stepSuper; supidx++){
        SuperNodeBufferTypeUnsym & snode = arrSuperNodes[supidx];
        sendOffset[supidx]=sendCount;
        recvOffset[supidx]=recvCount;
        sendCount+= this->CountSendToCrossDiagonal(snode.Index);
        recvCount+= this->CountRecvFromCrossDiagonal(snode.Index);
      }

      //Buffers for L
      std::vector<MPI_Request > arrMpiReqsSendLCD(sendCount, MPI_REQUEST_NULL );
      std::vector<MPI_Request > arrMpiReqsSizeSendLCD(sendCount, MPI_REQUEST_NULL );
      std::vector<MPI_Request > arrMpiReqsRecvLCD(recvCount, MPI_REQUEST_NULL );
      std::vector<MPI_Request > arrMpiReqsSizeRecvLCD(recvCount, MPI_REQUEST_NULL );
      std::vector<std::vector<char> > arrSstrLcolSendCD(sendCount);
      std::vector<int > arrSstrLcolSizeSendCD(sendCount);
      std::vector<std::vector<char> > arrSstrLrowRecvCD(recvCount);
      std::vector<int > arrSstrLrowSizeRecvCD(recvCount);

      //Buffers for U
      std::vector<MPI_Request > arrMpiReqsSendUCD(recvCount, MPI_REQUEST_NULL );
      std::vector<MPI_Request > arrMpiReqsSizeSendUCD(recvCount, MPI_REQUEST_NULL );
      std::vector<MPI_Request > arrMpiReqsRecvUCD(sendCount, MPI_REQUEST_NULL );
      std::vector<MPI_Request > arrMpiReqsSizeRecvUCD(sendCount, MPI_REQUEST_NULL );
      std::vector<std::vector<char> > arrSstrUrowSendCD(recvCount);
      std::vector<int > arrSstrUrowSizeSendCD(recvCount);
      std::vector<std::vector<char> > arrSstrUcolRecvCD(sendCount);
      std::vector<int > arrSstrUcolSizeRecvCD(sendCount);



      //Do Isend for size and content of L and Irecv for sizes of U 
      for (Int supidx=0; supidx<stepSuper; supidx++){
        SuperNodeBufferTypeUnsym & snode = arrSuperNodes[supidx];

        TIMER_START(Send_L_Recv_Size_U_CrossDiag);

        if( MYCOL( this->grid_ ) == PCOL( snode.Index, this->grid_ ) 
            && this->isSendToCrossDiagonal_(this->grid_->numProcCol, snode.Index ) ){

          Int sendIdxL = 0;
          Int recvIdxU = 0;
          for(Int dstCol = 0; dstCol<this->grid_->numProcCol; dstCol++){
            if(this->isSendToCrossDiagonal_(dstCol,snode.Index) ){
              Int dest = PNUM(PROW(snode.Index,this->grid_),dstCol,this->grid_);

              if( MYPROC( this->grid_ ) != dest	){
                //Send size and content of L
                MPI_Request & mpiReqSizeSendL = 
                  arrMpiReqsSizeSendLCD[sendOffset[supidx]+sendIdxL];
                MPI_Request & mpiReqSendL = 
                  arrMpiReqsSendLCD[sendOffset[supidx]+sendIdxL];
                std::vector<char> & sstrLcolSend = 
                  arrSstrLcolSendCD[sendOffset[supidx]+sendIdxL];
                Int & sstrSizeL = 
                  arrSstrLcolSizeSendCD[sendOffset[supidx]+sendIdxL];

                // Pack L data
                std::stringstream sstm;
                std::vector<Int> mask( LBlockMask::TOTAL_NUMBER, 1 );
                std::vector<LBlock<T> >&  Lcol = this->L( LBj(snode.Index, this->grid_) );
                // All blocks except for the diagonal block are to be sent right
                //TODO not true > this is a scatter operation ! Can we know the destination ?

                //Skip the diagonal block if necessary 
                Int startIdx = ( MYROW( this->grid_ ) == PROW( snode.Index, this->grid_ ) ) ? 1:0;
                Int count = 0;
                for( Int ib = startIdx; ib < Lcol.size(); ib++ ){
                  if( Lcol[ib].blockIdx > snode.Index && 
                      (Lcol[ib].blockIdx % this->grid_->numProcCol) == dstCol  ){
                    count++;
                  }
                }

                serialize( (Int)count, sstm, NO_MASK );

                for( Int ib = startIdx; ib < Lcol.size(); ib++ ){
                  if( Lcol[ib].blockIdx > snode.Index &&  
                      (Lcol[ib].blockIdx % this->grid_->numProcCol) == dstCol  ){ 
#if ( _DEBUGlevel_ >= 1 )
                    statusOFS << "["<<snode.Index<<"] SEND contains "<<Lcol[ib].blockIdx
                      << " which corresponds to "<<GBj(ib,this->grid_)
                      << " to P"<<dest
                      << std::endl;
#endif
                    serialize( Lcol[ib], sstm, mask );
                  }
                }

                sstrLcolSend.resize( Size(sstm) );
                sstm.read( &sstrLcolSend[0], sstrLcolSend.size() );
                sstrSizeL = sstrLcolSend.size();



                MPI_Isend( &sstrSizeL, 1, MPI_INT, dest, 
                    IDX_TO_TAG2(snode.Index,supidx,SELINV_TAG_L_SIZE_CD),
                    this->grid_->comm, &mpiReqSizeSendL );
                MPI_Isend( (void*)&sstrLcolSend[0], sstrSizeL, MPI_BYTE, dest,
                    IDX_TO_TAG2(snode.Index,supidx,SELINV_TAG_L_CONTENT_CD),
                    this->grid_->comm, &mpiReqSendL );
                sendIdxL++;

                //Recv for U size
                Int & sstrSizeU = 
                  arrSstrUcolSizeRecvCD[sendOffset[supidx]+recvIdxU];
                MPI_Request & mpiReqSizeRecvU = 
                  arrMpiReqsSizeRecvUCD[sendOffset[supidx]+recvIdxU];

                MPI_Irecv( &sstrSizeU, 1, MPI_INT, dest, 
                    IDX_TO_TAG2(snode.Index,supidx,SELINV_TAG_U_SIZE_CD),
                    this->grid_->comm, &mpiReqSizeRecvU );
                recvIdxU++;
              }
            }
          }
        } // sender
        TIMER_STOP(Send_L_Recv_Size_U_CrossDiag);
      }


      //Do Irecv for sizes of L and Isend for size and content of U
      for (Int supidx=0; supidx<stepSuper; supidx++){
        SuperNodeBufferTypeUnsym & snode = arrSuperNodes[supidx];
        TIMER_START(Send_U_Recv_Size_L_CrossDiag);
        //If I'm a receiver
        if( MYROW( this->grid_ ) == PROW( snode.Index, this->grid_ ) 
            && this->isRecvFromCrossDiagonal_(this->grid_->numProcRow, snode.Index ) ){
          Int recvIdxL=0;
          Int sendIdxU=0;
          for(Int srcRow = 0; srcRow<this->grid_->numProcRow; srcRow++){
            if(this->isRecvFromCrossDiagonal_(srcRow,snode.Index) ){
              Int src = PNUM(srcRow,PCOL(snode.Index,this->grid_),this->grid_);
              if( MYPROC( this->grid_ ) != src ){
                //Recv size of L
                Int & sstrSizeL = 
                  arrSstrLrowSizeRecvCD[recvOffset[supidx]+recvIdxL];
                MPI_Request & mpiReqSizeRecvL = 
                  arrMpiReqsSizeRecvLCD[recvOffset[supidx]+recvIdxL];

                MPI_Irecv( &sstrSizeL, 1, MPI_INT, src, 
                    IDX_TO_TAG2(snode.Index,supidx,SELINV_TAG_L_SIZE_CD),
                    this->grid_->comm, &mpiReqSizeRecvL );
                recvIdxL++;


                //Send size and content of U
                MPI_Request & mpiReqSizeSendU = 
                  arrMpiReqsSizeSendUCD[recvOffset[supidx]+sendIdxU];
                MPI_Request & mpiReqSendU = 
                  arrMpiReqsSendUCD[recvOffset[supidx]+sendIdxU];
                std::vector<char> & sstrUrowSend = 
                  arrSstrUrowSendCD[recvOffset[supidx]+sendIdxU];
                Int & sstrSizeU = 
                  arrSstrUrowSizeSendCD[recvOffset[supidx]+sendIdxU];

                // Pack U data
                std::stringstream sstm;
                std::vector<Int> mask( UBlockMask::TOTAL_NUMBER, 1 );
                std::vector<UBlock<T> >&  Urow = this->U( LBi(snode.Index, this->grid_) );
                // All blocks except for the diagonal block are to be sent right
                //TODO not true > this is a scatter operation ! Can we know the destination ?

                Int count = 0;
                for( Int jb = 0; jb < Urow.size(); jb++ ){
                  if( Urow[jb].blockIdx > snode.Index 
                      && (Urow[jb].blockIdx % this->grid_->numProcRow) == srcRow ){
                    count++;
                  }
                }

                serialize( (Int)count, sstm, NO_MASK );

                for( Int jb = 0; jb < Urow.size(); jb++ ){
                  if( Urow[jb].blockIdx > snode.Index 
                      && (Urow[jb].blockIdx % this->grid_->numProcRow) == srcRow ){ 
#if ( _DEBUGlevel_ >= 1 )
                    statusOFS << "["<<snode.Index<<"] SEND contains "<<Urow[jb].blockIdx
                      << " which corresponds to "<<GBi(jb,this->grid_)
                      << " to P"<<src
                      << std::endl;
#endif
                    serialize( Urow[jb], sstm, mask );
                  }
                }

                sstrUrowSend.resize( Size(sstm) );
                sstm.read( &sstrUrowSend[0], sstrUrowSend.size() );
                sstrSizeU = sstrUrowSend.size();


                MPI_Isend( &sstrSizeU, 1, MPI_INT, src, 
                    IDX_TO_TAG2(snode.Index,supidx,SELINV_TAG_U_SIZE_CD), 
                    this->grid_->comm, &mpiReqSizeSendU );
                MPI_Isend( (void*)&sstrUrowSend[0], sstrSizeU, MPI_BYTE, src,
                    IDX_TO_TAG2(snode.Index,supidx,SELINV_TAG_U_CONTENT_CD),
                    this->grid_->comm, &mpiReqSendU );
                sendIdxU++;
              }
            }
          }
        }//end if I'm a receiver
        TIMER_STOP(Send_U_Recv_Size_L_CrossDiag);
      }

      //waitall sizes of L
      TIMER_START(Wait_Size_L_CrossDiag);
      mpi::Waitall(arrMpiReqsSizeRecvLCD);
      TIMER_STOP(Wait_Size_L_CrossDiag);

      //if(this->grid_->mpirank==0){gdb_lock();}

      //Allocate content and do Irecv for content of L
      for (Int supidx=0; supidx<stepSuper; supidx++){
        SuperNodeBufferTypeUnsym & snode = arrSuperNodes[supidx];
        //If I'm a receiver
        TIMER_START(Recv_L_CrossDiag);
        if( MYROW( this->grid_ ) == PROW( snode.Index, this->grid_ ) 
            && this->isRecvFromCrossDiagonal_(this->grid_->numProcRow, snode.Index ) ){
          Int recvIdxL=0;
          for(Int srcRow = 0; srcRow<this->grid_->numProcRow; srcRow++){
            if(this->isRecvFromCrossDiagonal_(srcRow,snode.Index) ){
              Int src = PNUM(srcRow,PCOL(snode.Index,this->grid_),this->grid_);
              if( MYPROC( this->grid_ ) != src ){
                Int & sstrSizeL = 
                  arrSstrLrowSizeRecvCD[recvOffset[supidx]+recvIdxL];
                std::vector<char> & sstrLrowRecv = 
                  arrSstrLrowRecvCD[recvOffset[supidx]+recvIdxL];
                MPI_Request & mpiReqRecvL = 
                  arrMpiReqsRecvLCD[recvOffset[supidx]+recvIdxL];
                sstrLrowRecv.resize( sstrSizeL);

                MPI_Irecv( (void*)&sstrLrowRecv[0], sstrSizeL, MPI_BYTE, src,
                    IDX_TO_TAG2(snode.Index,supidx,SELINV_TAG_L_CONTENT_CD),
                    this->grid_->comm, &mpiReqRecvL );
                recvIdxL++;

              }
            }
          }
        }//end if I'm a receiver
        TIMER_STOP(Recv_L_CrossDiag);
      }

      //waitall sizes of U
      TIMER_START(Wait_Size_U_CrossDiag);
      mpi::Waitall(arrMpiReqsSizeRecvUCD);
      TIMER_STOP(Wait_Size_U_CrossDiag);

      //Allocate content and do Irecv for content of U
      for (Int supidx=0; supidx<stepSuper; supidx++){
        SuperNodeBufferTypeUnsym & snode = arrSuperNodes[supidx];
        TIMER_START(Recv_U_CrossDiag);
        if( MYCOL( this->grid_ ) == PCOL( snode.Index, this->grid_ ) 
            && this->isSendToCrossDiagonal_(this->grid_->numProcCol, snode.Index ) ){
          Int recvIdxU = 0;
          for(Int srcCol = 0; srcCol<this->grid_->numProcCol; srcCol++){
            if(this->isSendToCrossDiagonal_(srcCol,snode.Index) ){
              Int src = PNUM(PROW(snode.Index,this->grid_),srcCol,this->grid_);
              if( MYPROC( this->grid_ ) != src ){
                Int & sstrSizeU = 
                  arrSstrUcolSizeRecvCD[sendOffset[supidx]+recvIdxU];
                std::vector<char> & sstrUcolRecv = 
                  arrSstrUcolRecvCD[sendOffset[supidx]+recvIdxU];
                MPI_Request & mpiReqRecvU = 
                  arrMpiReqsRecvUCD[sendOffset[supidx]+recvIdxU];
                sstrUcolRecv.resize( sstrSizeU );

                MPI_Irecv( (void*)&sstrUcolRecv[0], sstrSizeU, MPI_BYTE, src,
                    IDX_TO_TAG2(snode.Index,supidx,SELINV_TAG_U_CONTENT_CD),
                    this->grid_->comm, &mpiReqRecvU );
                recvIdxU++;
              }
            }
          }
        }//end if I'm a receiver
        TIMER_STOP(Recv_U_CrossDiag);
      }

      //waitall content of L
      TIMER_START(Wait_L_CrossDiag);
      mpi::Waitall(arrMpiReqsRecvLCD);
      TIMER_STOP(Wait_L_CrossDiag);



      //Do the work for Lrow
      NumMat<T> Ltmp;
      for (Int supidx=0; supidx<stepSuper; supidx++){
        SuperNodeBufferTypeUnsym & snode = arrSuperNodes[supidx];
        Int ksupProcRow = PROW( snode.Index, this->grid_ );
        Int ksupProcCol = PCOL( snode.Index, this->grid_ );

        if( MYROW( this->grid_ ) == PROW( snode.Index, this->grid_ ) &&
            this->isRecvFromCrossDiagonal_(this->grid_->numProcRow, snode.Index ) ){

#if ( _DEBUGlevel_ >= 1 )
          statusOFS << std::endl << " ["<<snode.Index<<"] "
            << "Update the upper triangular block" 
            << std::endl << std::endl; 
          statusOFS << std::endl << " ["<<snode.Index<<"] "
            << "blockIdxLocal:" << snode.BlockIdxLocal
            << std::endl << std::endl; 
          statusOFS << std::endl << " ["<<snode.Index<<"] "
            << "rowLocalPtr:" << snode.RowLocalPtr
            << std::endl << std::endl; 
#endif

          std::vector<LBlock<T> >& Lcol = this->L( LBj(snode.Index, this->grid_) );
          std::vector<LBlock<T> >& Lrow = this->Lrow( LBi( snode.Index, this->grid_ ) );

          //          if(Lrow.size() == 0){
          //            Int& LrowSize = this->LrowSize_[ LBi( snode.Index, this->grid_ ) ];
          //            Lrow.resize(LrowSize);
          //          }

          std::vector<Int> isBlockFound(Lrow.size(),false);
          Int recvIdxL=0;

          for(Int srcRow = 0; srcRow<this->grid_->numProcRow; srcRow++){
            if(this->isRecvFromCrossDiagonal_(srcRow,snode.Index) ){
              Int src = PNUM(srcRow,PCOL(snode.Index,this->grid_),this->grid_);
              TIMER_START(Recv_L_CrossDiag);

              //Declaring some pointers to avoid copy if data is local
              std::vector<LBlock<T> > * pLcol;
              std::vector<LBlock<T> > LcolRecv;
              if(MYPROC(this->grid_)!= src){
                Int & sstrSizeL = 
                  arrSstrLrowSizeRecvCD[recvOffset[supidx]+recvIdxL];
                std::vector<char> & sstrLrowRecv = 
                  arrSstrLrowRecvCD[recvOffset[supidx]+recvIdxL];

                std::stringstream sstm;

#if ( _DEBUGlevel_ >= 1 )
                statusOFS << "["<<snode.Index<<"] P"<<MYPROC(this->grid_)<<" ("<<MYROW(this->grid_)
                  << ","<<MYCOL(this->grid_)<<") <--- LBj("<<snode.Index<<") <--- P"
                  << src <<" ("<<srcRow<<","<<ksupProcCol<<")"
                  << std::endl;
#endif
                sstm.write( &sstrLrowRecv[0], sstrSizeL );

                // Unpack L data.  
                Int numLBlock;
                std::vector<Int> mask( LBlockMask::TOTAL_NUMBER, 1 );
                deserialize( numLBlock, sstm, NO_MASK );
                LcolRecv.resize(numLBlock);
                for( Int ib = 0; ib < numLBlock; ib++ ){
                  deserialize( LcolRecv[ib], sstm, mask );
#if ( _DEBUGlevel_ >= 1 )
                  statusOFS << "["<<snode.Index<<"] RECV contains "
                    << LcolRecv[ib].blockIdx<< " which corresponds to "
                    << ib * this->grid_->numProcRow + srcRow;
                  statusOFS << " L is on row "<< srcRow 
                    << " whereas Lrow is on col "
                    << LcolRecv[ib].blockIdx % this->grid_->numProcCol 
                    << std::endl;
#endif
                }
                recvIdxL++;
                pLcol = &LcolRecv;
              } // sender is not the same as receiver
              else{
                // L is obtained locally, just make a copy. Do not include the diagonal block
                pLcol = &Lcol;
              } // sender is the same as receiver


              TIMER_STOP(Recv_L_CrossDiag);

              //We can directly put LcolRecv in Lrow (stil need to transpose)

              //Skip the diagonal block if necessary 
              Int startIdx = ( MYPROC( this->grid_ ) == src && MYPROC(this->grid_) == PNUM(ksupProcRow,ksupProcCol, this->grid_) ) ? 1:0;
              for( Int ib = startIdx; ib < pLcol->size(); ib++ ){
                LBlock<T> & LB = (*pLcol)[ib];
                if( LB.blockIdx <= snode.Index ){
                  ErrorHandling( "LcolRecv contains the wrong blocks." );
                }

                //check that LB would be on this proc if it was a UB
                if(PCOL(LB.blockIdx, this->grid_) != MYCOL(this->grid_)){
                  continue;
                }



                //std::vector<LBlock<T> > &  LrowCD = this->Lrow( LBi( ksup, super_ ) );
                for( Int iib = 0; iib < Lrow.size(); iib++ ){
                  LBlock<T> &  LrowB = Lrow[ iib ];

                  //If not initialized, set it to be the one
                  if( LrowB.blockIdx == -1){
                    LrowB.blockIdx = LB.blockIdx;
                  }

                  if( LB.blockIdx == LrowB.blockIdx ){
                    // Note that the order of the column indices of the U
                    // block may not follow the order of the row indices,
                    // overwrite the information in U.
                    //                    LrowB = LB;
                    LrowB.rows = LB.rows;
                    //Store in "row major" format / i.e transpose is in col-major
                    Transpose(LB.nzval, LrowB.nzval);
                    LrowB.numCol = LB.numRow;
                    LrowB.numRow = LB.numCol;

                    //                    SetValue(LrowB.nzval,ZERO<T>());

#if ( _DEBUGlevel_ >= 1 )
                    statusOFS<<"["<<snode.Index<<"] USING LB "<<LB.blockIdx<< std::endl;
#endif
                    isBlockFound[iib] = 1;//true;
                    break;
                  } // if( LB.blockIdx == LrowB.blockIdx )
                } // for (iib)
              } // for (ib)
            }
          }

          for( Int ib = 0; ib < Lrow.size(); ib++ ){
            LBlock<T> &  LB = Lrow[ib];
            if( !isBlockFound[ib] ){
              ErrorHandling( 
                  "LBlock cannot find its update. Something is seriously wrong."
                  );
            }
          }
        } // Done copying L
      }

      //waitall content of U
      TIMER_START(Wait_U_CrossDiag);
      mpi::Waitall(arrMpiReqsRecvUCD);
      TIMER_STOP(Wait_U_CrossDiag);


      //Do the work for U
      for (Int supidx=0; supidx<stepSuper; supidx++){
        SuperNodeBufferTypeUnsym & snode = arrSuperNodes[supidx];
        Int ksupProcRow = PROW( snode.Index, this->grid_ );
        Int ksupProcCol = PCOL( snode.Index, this->grid_ );


        if( MYCOL( this->grid_ ) == PCOL( snode.Index, this->grid_ ) &&
            this->isSendToCrossDiagonal_(this->grid_->numProcCol, snode.Index ) ){

          std::vector<UBlock<T> >& Ucol = this->Ucol( LBj( snode.Index, this->grid_ ) );

          //          if(Ucol.size() == 0){
          //            Int& UcolSize = this->UcolSize_[ LBj( snode.Index, this->grid_ ) ];
          //            Ucol.resize(UcolSize);
          //          }

          std::vector<Int> isBlockFound(Ucol.size(),false);
          Int recvIdxU=0;

          for(Int srcCol = 0; srcCol<this->grid_->numProcCol; srcCol++){
            if(this->isSendToCrossDiagonal_(srcCol,snode.Index) ){
              Int src = PNUM(PROW(snode.Index,this->grid_),srcCol,this->grid_);
              TIMER_START(Recv_U_CrossDiag);

              //Declaring some pointers to avoid copy if data is local
              std::vector<UBlock<T> > * pUrow;
              std::vector<UBlock<T> > UrowRecv;
              if(MYPROC(this->grid_)!= src){
                Int & sstrSizeU = 
                  arrSstrUcolSizeRecvCD[sendOffset[supidx]+recvIdxU];
                std::vector<char> & sstrUcolRecv = 
                  arrSstrUcolRecvCD[sendOffset[supidx]+recvIdxU];


                std::stringstream sstm;

#if ( _DEBUGlevel_ >= 1 )
                statusOFS << "["<<snode.Index<<"] P"<<MYPROC(this->grid_)<<" ("<<MYROW(this->grid_)
                  << ","<<MYCOL(this->grid_)<<") <--- LBi("<<snode.Index<<") <--- P"
                  << src<<" ("<<ksupProcCol<<","<<srcCol<<")"<<std::endl;
#endif
                sstm.write( &sstrUcolRecv[0], sstrSizeU );

                // Unpack U data.  
                Int numUBlock;
                std::vector<Int> mask( UBlockMask::TOTAL_NUMBER, 1 );
                deserialize( numUBlock, sstm, NO_MASK );
                UrowRecv.resize(numUBlock);
                for( Int jb = 0; jb < numUBlock; jb++ ){
                  deserialize( UrowRecv[jb], sstm, mask );
#if ( _DEBUGlevel_ >= 1 )
                  statusOFS << "["<<snode.Index<<"] RECV contains "
                    << UrowRecv[jb].blockIdx<< " which corresponds to "
                    << jb * this->grid_->numProcCol + srcCol;
                  statusOFS << " Urow is on col "<< srcCol 
                    << " whereas U is on row "
                    << UrowRecv[jb].blockIdx % this->grid_->numProcRow 
                    << std::endl;
#endif
                }
                pUrow = &UrowRecv;
                recvIdxU++;
              } // sender is not the same as receiver
              else{
                // U is obtained locally, just make a copy. 
                std::vector<UBlock<T> >& Urow = this->U( LBi( snode.Index, this->grid_ ) );
                pUrow = &Urow;
              } // sender is the same as receiver


              TIMER_STOP(Recv_U_CrossDiag);

              //We can directly put UrowRecv in Ucol
              for( Int jb = 0; jb < pUrow->size(); jb++ ){
                UBlock<T> & UB = (*pUrow)[jb];
                if( UB.blockIdx <= snode.Index ){
                  ErrorHandling( "UrowRecv contains the wrong blocks." );
                }

                //check that UB would be on this proc if it was a LB
                if(PROW(UB.blockIdx, this->grid_) != MYROW(this->grid_)){
                  continue;
                }


                //std::vector<UBlock<T> > &  UcolCD = this->Ucol( LBi( ksup, super_ ) );
                for( Int jjb = 0; jjb < Ucol.size(); jjb++ ){
                  UBlock<T> &  UcolB = Ucol[jjb];


                  //If not initialized, set it to be the one
                  if( UcolB.blockIdx == -1){
                    UcolB.blockIdx = UB.blockIdx;
                  }


                  if( UB.blockIdx == UcolB.blockIdx ){
                    // Compare size
                    // Note that the order of the column indices of the U
                    // block may not follow the order of the row indices,
                    // overwrite the information in U.
                    UcolB = UB;

                    //SetValue(UcolB.nzval,ZERO<T>());
#if ( _DEBUGlevel_ >= 1 )
                    statusOFS<<"["<<snode.Index<<"] USING UB "<<UB.blockIdx<< std::endl;
#endif
                    isBlockFound[jjb] = 1;
                    break;
                  } // if( UB.blockIdx == UcolB.blockIdx )
                } // for (jjb)
              } // for (jb)

            }
          }

          for( Int jb = 0; jb < Ucol.size(); jb++ ){
            UBlock<T> &  UB = Ucol[jb];
            if( !isBlockFound[jb] ){
              ErrorHandling( 
                  "UBlock cannot find its update. Something is seriously wrong."
                  );
            }
          } // end for (jb)
        } // Done copying U

      }

      TIMER_STOP(Send_CD);

      mpi::Waitall(arrMpiReqsSizeSendLCD);
      mpi::Waitall(arrMpiReqsSendLCD);

      mpi::Waitall(arrMpiReqsSizeSendUCD);
      mpi::Waitall(arrMpiReqsSendUCD);


    } // End of method PMatrixUnsym::SendRecvCD

  template<typename T>
    inline void PMatrixUnsym<T>::UnpackData(
        SuperNodeBufferTypeUnsym & snode, 
        std::vector<LBlock<T> > & LcolRecv, 
        std::vector<LBlock<T> > & LrowRecv,
        std::vector<UBlock<T> > & UcolRecv, 
        std::vector<UBlock<T> > & UrowRecv
        )
    {

      TIMER_START(Unpack_data);




#if ( _DEBUGlevel_ >= 1 )
      statusOFS << std::endl << "["<<snode.Index<<"] "<<  "Unpack the received data for processors participate in Gemm. " << std::endl << std::endl; 
#endif
      // Lrow part
      if( MYROW( this->grid_ ) != PROW( snode.Index, this->grid_ ) ){
        std::stringstream sstm;
        sstm.write( &snode.SstrLrowRecv[0], snode.SstrLrowRecv.size() );
        std::vector<Int> mask( LBlockMask::TOTAL_NUMBER, 1 );
        Int numLBlock;
        deserialize( numLBlock, sstm, NO_MASK );
        LrowRecv.resize( numLBlock );
        for( Int ib = 0; ib < numLBlock; ib++ ){
          deserialize( LrowRecv[ib], sstm, mask );
        }

#if ( _DEBUGlevel_ >= 1 )
        statusOFS<< "["<<snode.Index<<"] "<<"Lrow RECV "<<std::endl;
        for(Int ib=0;ib<LrowRecv.size();++ib){statusOFS<<LrowRecv[ib].blockIdx<<" ";}
        statusOFS<<std::endl;
#endif
      } // sender is not the same as receiver
      else{
        // U is obtained locally, just make a copy. Include everything
        // (there is no diagonal block)
        // Is it a copy ?  LL: YES. Maybe we should replace the copy by
        // something more efficient especially for mpisize == 1
        LrowRecv.resize(this->Lrow( LBi( snode.Index, this->grid_ ) ).size());
        std::copy(this->Lrow( LBi( snode.Index, this->grid_ ) ).begin(),this->Lrow( LBi( snode.Index, this->grid_ )).end(),LrowRecv.begin());
      } // sender is the same as receiver


      //L part
      if( MYCOL( this->grid_ ) != PCOL( snode.Index, this->grid_ ) ){
        std::stringstream     sstm;
        sstm.write( &snode.SstrLcolRecv[0], snode.SstrLcolRecv.size() );
        std::vector<Int> mask( LBlockMask::TOTAL_NUMBER, 1 );
        mask[LBlockMask::NZVAL] = 0; // nzval is excluded
        Int numLBlock;
        deserialize( numLBlock, sstm, NO_MASK );
        LcolRecv.resize( numLBlock );
        for( Int ib = 0; ib < numLBlock; ib++ ){
          deserialize( LcolRecv[ib], sstm, mask );
        }


#if ( _DEBUGlevel_ >= 1 )
        statusOFS<< "["<<snode.Index<<"] "<<"Lcol RECV "<<std::endl;
        for(Int ib=0;ib<LcolRecv.size();++ib){statusOFS<<LcolRecv[ib].blockIdx<<" ";}
        statusOFS<<std::endl;
#endif

      } // sender is not the same as receiver
      else{
        // L is obtained locally, just make a copy. 
        // Do not include the diagonal block
        std::vector<LBlock<T> >& Lcol =  this->L( LBj( snode.Index, this->grid_ ) );
        Int startIdx = ( MYROW( this->grid_ ) == PROW( snode.Index, this->grid_ ) )?1:0;
        LcolRecv.resize( Lcol.size() - startIdx );
        std::copy(Lcol.begin()+startIdx,Lcol.end(),LcolRecv.begin());
      } // sender is the same as receiver

      // U part
      if( MYROW( this->grid_ ) != PROW( snode.Index, this->grid_ ) ){
        std::stringstream sstm;
        sstm.write( &snode.SstrUrowRecv[0], snode.SstrUrowRecv.size() );
        std::vector<Int> mask( UBlockMask::TOTAL_NUMBER, 1 );
        mask[UBlockMask::NZVAL] = 0; // nzval is excluded
        Int numUBlock;
        deserialize( numUBlock, sstm, NO_MASK );
        UrowRecv.resize( numUBlock );
        for( Int jb = 0; jb < numUBlock; jb++ ){
          deserialize( UrowRecv[jb], sstm, mask );
        } 

#if ( _DEBUGlevel_ >= 1 )
        statusOFS<< "["<<snode.Index<<"] "<<"Urow RECV "<<std::endl;
        for(Int ib=0;ib<UrowRecv.size();++ib){statusOFS<<UrowRecv[ib].blockIdx<<" ";}
        statusOFS<<std::endl;
#endif
      } // sender is not the same as receiver
      else{
        // U is obtained locally, just make a copy. Include everything
        // (there is no diagonal block)
        // Is it a copy ?  LL: YES. Maybe we should replace the copy by
        // something more efficient especially for mpisize == 1
        UrowRecv.resize(this->U( LBi( snode.Index, this->grid_ ) ).size());
        std::copy(this->U( LBi( snode.Index, this->grid_ ) ).begin(),this->U( LBi( snode.Index, this->grid_ )).end(),UrowRecv.begin());
      } // sender is the same as receiver


      // Ucol part
      if( MYCOL( this->grid_ ) != PCOL( snode.Index, this->grid_ ) ){
        std::stringstream     sstm;
        sstm.write( &snode.SstrUcolRecv[0], snode.SstrUcolRecv.size() );
        std::vector<Int> mask( UBlockMask::TOTAL_NUMBER, 1 );
        Int numUBlock;
        deserialize( numUBlock, sstm, NO_MASK );
        UcolRecv.resize( numUBlock );
        for( Int jb = 0; jb < numUBlock; jb++ ){
          deserialize( UcolRecv[jb], sstm, mask );
        }

#if ( _DEBUGlevel_ >= 1 )
        statusOFS<< "["<<snode.Index<<"] "<<"Ucol RECV "<<std::endl;
        for(Int ib=0;ib<UcolRecv.size();++ib){statusOFS<<UcolRecv[ib].blockIdx<<" ";}
        statusOFS<<std::endl;
#endif
      } // sender is not the same as receiver
      else{
        // Ucol is obtained locally, just make a copy. 
        // Do not include the diagonal block
        std::vector<UBlock<T> >& Ucol =  this->Ucol( LBj( snode.Index, this->grid_ ) );
        UcolRecv.resize( Ucol.size() );
        std::copy(Ucol.begin(),Ucol.end(),UcolRecv.begin());
      } // sender is the same as receiver

      TIMER_STOP(Unpack_data);

    } // End of method PMatrixUnsym<T>::UnpackData

  template<typename T>
    inline void PMatrixUnsym<T>::ComputeDiagUpdate(SuperNodeBufferTypeUnsym & snode)
    {

      TIMER_START(ComputeDiagUpdate);
      //--------- Computing  Diagonal block, all processors in the column
      //--------- are participating to all pipelined supernodes
      if( MYCOL( this->grid_ ) == PCOL( snode.Index, this->grid_ ) ){
#if ( _DEBUGlevel_ >= 1 )
        statusOFS << std::endl << "["<<snode.Index<<"] "
          << "Updating the diagonal block" << std::endl << std::endl;
#endif
        std::vector<LBlock<T> >&  Lcol = this->L( LBj( snode.Index, this->grid_ ) );
        std::vector<UBlock<T> >&  Ucol = this->Ucol( LBj( snode.Index, this->grid_ ) );

        //Allocate DiagBuf even if Lcol.size() == 0
        snode.DiagBuf.Resize(SuperSize( snode.Index, this->super_ ),
            SuperSize( snode.Index, this->super_ ));
        SetValue(snode.DiagBuf, ZERO<T>());

        // Do I own the diagonal block ?
        Int startIb = (MYROW( this->grid_ ) == PROW( snode.Index, this->grid_ ))?1:0;
        for( Int ib = startIb; ib < Lcol.size(); ib++ ){
          LBlock<T> & LcolB = Lcol[ib];
          //find corresponding U block
          Int jb = 0;
          for(jb = 0; jb < Ucol.size(); ++jb){
            if(Ucol[jb].blockIdx == LcolB.blockIdx){
              break;
            }
          }


          assert(jb < Ucol.size());

          UBlock<T> & UcolB = Ucol[jb];

          if(1 || (LcolB.numRow>0 && LcolB.numCol>0 && UcolB.numRow>0 && UcolB.numCol>0)){
            //Compute U S-1 L
            blas::Gemm( 'N', 'N', snode.DiagBuf.m(), snode.DiagBuf.n(), 
                LcolB.numRow, MINUS_ONE<T>(),
                UcolB.nzval.Data(), UcolB.nzval.m(), 
                &snode.LUpdateBuf( snode.RowLocalPtr[ib-startIb], 0 ),
                snode.LUpdateBuf.m(), 
                ONE<T>(), snode.DiagBuf.Data(), snode.DiagBuf.m() );
          }
        } 

#if ( _DEBUGlevel_ >= 1 )
        statusOFS << std::endl << "["<<snode.Index<<"] "
          << "Updated the diagonal block" << std::endl << std::endl; 
#endif
      }
      TIMER_STOP(ComputeDiagUpdate);
    } // End of method PMatrixUnsym<T>::ComputeDiagUpdate 








  template<typename T>
    inline void PMatrixUnsym<T>::SelInvIntra_P2p(Int lidx)
    {

#if defined (PROFILE) || defined(PMPI) || defined(USE_TAU)
      Real begin_SendULWaitContentFirst, end_SendULWaitContentFirst, time_SendULWaitContentFirst = 0;
#endif
      Int numSuper = this->NumSuper(); 
      std::vector<std::vector<Int> > & superList = this->WorkingSet();
      Int numSteps = superList.size();
      Int stepSuper = superList[lidx].size(); 

      TIMER_START(AllocateBuffer);

      //This is required to send the size and content of U/L
      std::vector<std::vector<MPI_Request> >  arrMpireqsSendLToBelow;
      arrMpireqsSendLToBelow.resize( stepSuper, std::vector<MPI_Request>( 2 * this->grid_->numProcRow, MPI_REQUEST_NULL ));
      std::vector<std::vector<MPI_Request> >  arrMpireqsSendLToRight;
      arrMpireqsSendLToRight.resize(stepSuper, std::vector<MPI_Request>( 2 * this->grid_->numProcCol, MPI_REQUEST_NULL ));

      std::vector<std::vector<MPI_Request> >  arrMpireqsSendUToBelow;
      arrMpireqsSendUToBelow.resize( stepSuper, std::vector<MPI_Request>( 2 * this->grid_->numProcRow, MPI_REQUEST_NULL ));
      std::vector<std::vector<MPI_Request> >  arrMpireqsSendUToRight;
      arrMpireqsSendUToRight.resize(stepSuper, std::vector<MPI_Request>( 2 * this->grid_->numProcCol, MPI_REQUEST_NULL ));

      //This is required to reduce L
      std::vector<MPI_Request>  arrMpireqsSendLToLeft;
      arrMpireqsSendLToLeft.resize(stepSuper, MPI_REQUEST_NULL );

      //This is required to reduce U
      std::vector<MPI_Request>  arrMpireqsSendUToAbove;
      arrMpireqsSendUToAbove.resize(stepSuper, MPI_REQUEST_NULL );

      //This is required to reduce D
      std::vector<MPI_Request>  arrMpireqsSendToAbove;
      arrMpireqsSendToAbove.resize(stepSuper, MPI_REQUEST_NULL );

      //------------------------------------------------------------------------

      //This is required to receive the size and content of U/L
      std::vector<MPI_Request>   arrMpireqsRecvLSizeFromAny;
      arrMpireqsRecvLSizeFromAny.resize(stepSuper*2 , MPI_REQUEST_NULL);
      std::vector<MPI_Request>   arrMpireqsRecvLContentFromAny;
      arrMpireqsRecvLContentFromAny.resize(stepSuper*2 , MPI_REQUEST_NULL);

      std::vector<MPI_Request>   arrMpireqsRecvUSizeFromAny;
      arrMpireqsRecvUSizeFromAny.resize(stepSuper*2 , MPI_REQUEST_NULL);
      std::vector<MPI_Request>   arrMpireqsRecvUContentFromAny;
      arrMpireqsRecvUContentFromAny.resize(stepSuper*2 , MPI_REQUEST_NULL);




      //allocate the buffers for this supernode
      std::vector<SuperNodeBufferTypeUnsym> arrSuperNodes(stepSuper);
      for (Int supidx=0; supidx<stepSuper; supidx++){ 
        arrSuperNodes[supidx].Index = superList[lidx][supidx];  
      }

      NumMat<T> AinvBuf, UBuf, LBuf;

      TIMER_STOP(AllocateBuffer);



      Int next_lidx = lidx+1;
      CDBuffers nextCDBuffers;


      //Perhaps this should be done for the next step super in a non blocking way
      if(lidx==0){
        CDBuffers buffers;
        Int next_lidx = lidx;

        //Resize L and U first
        bool doSend = false;
        for (Int supidx=0; supidx<superList[next_lidx].size(); supidx++){ 
          Int snode_index = superList[next_lidx][supidx];
          if( MYROW( this->grid_ ) == PROW( snode_index, this->grid_ ) ){
            std::vector<LBlock<T> >&  Lrow = this->Lrow( LBi(snode_index, this->grid_) );
            Int&  LrowSize = this->LrowSize_[ LBi(snode_index, this->grid_) ];
            if(Lrow.size()!=LrowSize){
              doSend= true;
              Lrow.resize(LrowSize,LBlock<T>());
            }
          }

          if( MYCOL( this->grid_ ) == PCOL( snode_index, this->grid_ ) ){
            std::vector<UBlock<T> >&  Ucol = this->Ucol( LBj(snode_index, this->grid_) );
            Int&  UcolSize = this->UcolSize_[ LBj(snode_index, this->grid_) ];
            if(Ucol.size()!=UcolSize){
              doSend= true;
              Ucol.resize(UcolSize,UBlock<T>());
            }
          }
        }


        if(doSend){

          SendRecvSizesCD(superList[next_lidx], superList[next_lidx].size(),buffers);
          IRecvContentCD(superList[next_lidx], superList[next_lidx].size(),buffers);
          WaitContentLCD(superList[next_lidx], superList[next_lidx].size(),buffers);
          WaitContentUCD(superList[next_lidx], superList[next_lidx].size(),buffers);
          buffers.WaitAllSend();

        }
      }

      if(next_lidx < superList.size()){
        //Resize L and U first
        for (Int supidx=0; supidx<superList[next_lidx].size(); supidx++){ 
          Int snode_index = superList[next_lidx][supidx];
          if( MYROW( this->grid_ ) == PROW( snode_index, this->grid_ ) ){
            std::vector<LBlock<T> >&  Lrow = this->Lrow( LBi(snode_index, this->grid_) );
            Int&  LrowSize = this->LrowSize_[ LBi(snode_index, this->grid_) ];
            Lrow.resize(LrowSize,LBlock<T>());
          }

          if( MYCOL( this->grid_ ) == PCOL( snode_index, this->grid_ ) ){
            std::vector<UBlock<T> >&  Ucol = this->Ucol( LBj(snode_index, this->grid_) );
            Int&  UcolSize = this->UcolSize_[ LBj(snode_index, this->grid_) ];
            Ucol.resize(UcolSize,UBlock<T>());
          }
        }

        SendRecvSizesCD(superList[next_lidx],superList[next_lidx].size(),nextCDBuffers);
      }


#if ( _DEBUGlevel_ >= 1 )
      statusOFS << std::endl << "Communication to the Schur complement." << std::endl << std::endl; 
#endif

      // Senders
      for (Int supidx=0; supidx<stepSuper; supidx++){
        SuperNodeBufferTypeUnsym & snode = arrSuperNodes[supidx];
        std::vector<MPI_Request> & mpireqsSendLToBelow = arrMpireqsSendLToBelow[supidx];
        std::vector<MPI_Request> & mpireqsSendLToRight = arrMpireqsSendLToRight[supidx];
        std::vector<MPI_Request> & mpireqsSendUToBelow = arrMpireqsSendUToBelow[supidx];
        std::vector<MPI_Request> & mpireqsSendUToRight = arrMpireqsSendUToRight[supidx];

#if ( _DEBUGlevel_ >= 1 )
        statusOFS << std::endl <<  "["<<snode.Index<<"] "
          << "Communication for the Lrow part." << std::endl << std::endl; 
#endif
        // Communication for the Lrow part.
        if( MYROW( this->grid_ ) == PROW( snode.Index, this->grid_ ) ){
          // Pack the data in Lrow
          TIMER_START(Serialize_LrowL);
          std::stringstream sstm;
          std::vector<Int> mask( LBlockMask::TOTAL_NUMBER, 1 );
          std::vector<LBlock<T> >&  Lrow = this->Lrow( LBi(snode.Index, this->grid_) );
          std::vector<UBlock<T> >&  Urow = this->U( LBi(snode.Index, this->grid_) );


          // All blocks are to be sent down.
          serialize( (Int)Lrow.size(), sstm, NO_MASK );
          for( Int ib = 0; ib < Lrow.size(); ib++ ){
            assert( Lrow[ib].blockIdx > snode.Index );
            serialize( Lrow[ib], sstm, mask );
          }
          snode.SstrLrowSend.resize( Size( sstm ) );
          sstm.read( &snode.SstrLrowSend[0], snode.SstrLrowSend.size() );
          snode.SizeSstrLrowSend = snode.SstrLrowSend.size();
          TIMER_STOP(Serialize_LrowL);

          for( Int iProcRow = 0; iProcRow < this->grid_->numProcRow; iProcRow++ ){
            if( MYROW( this->grid_ ) != iProcRow &&
                this->isSendToBelow_( iProcRow,snode.Index ) == true ){

              // Use Isend to send to multiple targets
              MPI_Isend( &snode.SizeSstrLrowSend, 1, MPI_INT,  
                  iProcRow, IDX_TO_TAG2(snode.Index,supidx,SELINV_TAG_LROW_SIZE), this->grid_->colComm, &mpireqsSendLToBelow[2*iProcRow] );
              MPI_Isend( (void*)&snode.SstrLrowSend[0], snode.SizeSstrLrowSend, MPI_BYTE, 
                  iProcRow, IDX_TO_TAG2(snode.Index,supidx,SELINV_TAG_LROW_CONTENT), 
                  this->grid_->colComm, &mpireqsSendLToBelow[2*iProcRow+1] );
#if ( _DEBUGlevel_ >= 1 )
              statusOFS<< "["<<snode.Index<<"] "<<"Lrow SENT "<<std::endl;
              for(Int ib=0;ib<Lrow.size();++ib){statusOFS<<Lrow[ib].blockIdx<<" ";}
              statusOFS<<std::endl;
              statusOFS << "["<<snode.Index<<"] "<<  "Sending Lrow "
                << snode.SizeSstrLrowSend << " BYTES to P" << PNUM(iProcRow,MYCOL(this->grid_),this->grid_)
                << std::endl <<  std::endl; 
#endif
            } // Send 
          } // for (iProcRow)
        } // if I am the sender

#if ( _DEBUGlevel_ >= 1 )
        statusOFS << std::endl << "["<<snode.Index<<"] "<< "Communication for the L part." << std::endl << std::endl; 
#endif
        // Communication for the L (Lcol) part.
        if( MYCOL( this->grid_ ) == PCOL( snode.Index, this->grid_ ) ){

          TIMER_START(Serialize_LcolL);
          // Pack the data in L 
          std::stringstream sstm;
          std::vector<Int> mask( LBlockMask::TOTAL_NUMBER, 1 );
          mask[LBlockMask::NZVAL] = 0; // nzval is excluded 

          std::vector<LBlock<T> >&  Lcol = this->L( LBj(snode.Index, this->grid_) );
          // All blocks except for the diagonal block are to be sent right


          Int startIdx = ( MYROW( this->grid_ ) == PROW( snode.Index, this->grid_ ) )?1:0;
          serialize( (Int)Lcol.size() - startIdx, sstm, NO_MASK );
          for( Int ib = startIdx; ib < Lcol.size(); ib++ ){
            assert( Lcol[ib].blockIdx > snode.Index );

#if ( _DEBUGlevel_ >= 2 )
            //                statusOFS << std::endl << "["<<snode.Index<<"] "<<  "Serializing Block index " << Lcol[ib].blockIdx << std::endl;
#endif
            serialize( Lcol[ib], sstm, mask );
          }
          snode.SstrLcolSend.resize( Size( sstm ) );
          sstm.read( &snode.SstrLcolSend[0], snode.SstrLcolSend.size() );
          snode.SizeSstrLcolSend = snode.SstrLcolSend.size();
          TIMER_STOP(Serialize_LcolL);

          for( Int iProcCol = 0; iProcCol < this->grid_->numProcCol ; iProcCol++ ){
            if( MYCOL( this->grid_ ) != iProcCol &&
                this->isSendToRight_( iProcCol, snode.Index ) == true ){
              // Use Isend to send to multiple targets
              MPI_Isend( &snode.SizeSstrLcolSend, 1, MPI_INT,  
                  iProcCol, IDX_TO_TAG2(snode.Index,supidx,SELINV_TAG_L_SIZE), 
                  this->grid_->rowComm, &mpireqsSendLToRight[2*iProcCol] );
              MPI_Isend( (void*)&snode.SstrLcolSend[0], snode.SizeSstrLcolSend, MPI_BYTE, 
                  iProcCol, IDX_TO_TAG2(snode.Index,supidx,SELINV_TAG_L_CONTENT), 
                  this->grid_->rowComm, &mpireqsSendLToRight[2*iProcCol+1] );
#if ( _DEBUGlevel_ >= 1 )
              statusOFS<< "["<<snode.Index<<"] "<<"L SENT "<<std::endl;
              for(Int ib=startIdx;ib<Lcol.size();++ib){statusOFS<<Lcol[ib].blockIdx<<" ";}
              statusOFS<<std::endl;
              statusOFS << "["<<snode.Index<<"] "<<  "Sending L "
                << snode.SizeSstrLcolSend << " BYTES to P" << PNUM(MYROW(this->grid_),iProcCol,this->grid_)
                << std::endl <<  std::endl; 
#endif
            } // Send 
          } // for (iProcCol)
        } // if I am the sender

#if ( _DEBUGlevel_ >= 1 )
        statusOFS << std::endl <<  "["<<snode.Index<<"] "
          << "Communication for the U part." << std::endl << std::endl; 
#endif
        // Communication for the U (Urow) part.
        if( MYROW( this->grid_ ) == PROW( snode.Index, this->grid_ ) ){

          TIMER_START(Serialize_UcolU);
          // Pack the data in U 
          std::stringstream sstm;
          std::vector<Int> mask( UBlockMask::TOTAL_NUMBER, 1 );
          mask[UBlockMask::NZVAL] = 0; // nzval is excluded 

          std::vector<UBlock<T> >&  Urow = this->U( LBi(snode.Index, this->grid_) );
          // All blocks except for the diagonal block are to be sent right

          serialize( (Int)Urow.size(), sstm, NO_MASK );
          for( Int jb = 0; jb < Urow.size(); jb++ ){
            assert( Urow[jb].blockIdx > snode.Index );
#if ( _DEBUGlevel_ >= 2 )
            //                statusOFS << std::endl << "["<<snode.Index<<"] "<<  "Serializing Block index " << Urow[jb].blockIdx << std::endl;
#endif
            serialize( Urow[jb], sstm, mask );
          }
          snode.SstrUrowSend.resize( Size( sstm ) );
          sstm.read( &snode.SstrUrowSend[0], snode.SstrUrowSend.size() );
          snode.SizeSstrUrowSend = snode.SstrUrowSend.size();
          TIMER_STOP(Serialize_UcolU);

          for( Int iProcRow = 0; iProcRow < this->grid_->numProcRow; iProcRow++ ){
            if( MYROW( this->grid_ ) != iProcRow &&
                this->isSendToBelow_( iProcRow,snode.Index ) == true ){
              // Use Isend to send to multiple targets
              MPI_Isend( &snode.SizeSstrUrowSend, 1, MPI_INT,  
                  iProcRow, IDX_TO_TAG2(snode.Index,supidx,SELINV_TAG_U_SIZE), this->grid_->colComm, &mpireqsSendUToBelow[2*iProcRow] );
              MPI_Isend( (void*)&snode.SstrLrowSend[0], snode.SizeSstrUrowSend, MPI_BYTE, 
                  iProcRow, IDX_TO_TAG2(snode.Index,supidx,SELINV_TAG_U_CONTENT), 
                  this->grid_->colComm, &mpireqsSendUToBelow[2*iProcRow+1] );
#if ( _DEBUGlevel_ >= 1 )
              statusOFS << std::endl << "["<<snode.Index<<"] "<<  "Sending U " << snode.SizeSstrUrowSend << " BYTES"<< std::endl <<  std::endl; 
              statusOFS<< "["<<snode.Index<<"] "<<"U SENT "<<std::endl;
              for(Int ib=0;ib<Urow.size();++ib){statusOFS<<Urow[ib].blockIdx<<" ";}
              statusOFS<<std::endl;
              statusOFS << "["<<snode.Index<<"] "<<  "Sending U "
                << snode.SizeSstrUrowSend << " BYTES to P" << PNUM(iProcRow,MYCOL(this->grid_),this->grid_)
                << std::endl <<  std::endl; 
#endif
            } // Send 
          } // for (iProcRow)
        } // if I am the sender

#if ( _DEBUGlevel_ >= 1 )
        statusOFS << "["<<snode.Index<<"] "
          << "Communication for the Ucol part." << std::endl 
          << std::endl; 
#endif
        // Communication for the Ucol part.
        if( MYCOL( this->grid_ ) == PCOL( snode.Index, this->grid_ ) ){
          // Pack the data in Ucol
          TIMER_START(Serialize_UcolU);
          std::stringstream sstm;
          std::vector<Int> mask( UBlockMask::TOTAL_NUMBER, 1 );
          std::vector<UBlock<T> >&  Ucol = 
            this->Ucol( LBj(snode.Index, this->grid_) );
          // All blocks are to be sent down.
          serialize( (Int)Ucol.size(), sstm, NO_MASK );
          for( Int jb = 0; jb < Ucol.size(); jb++ ){
            UBlock<T> & UB = Ucol[jb];
            assert( UB.blockIdx > snode.Index );
            serialize( UB, sstm, mask );
          }
          snode.SstrUcolSend.resize( Size( sstm ) );
          sstm.read( &snode.SstrUcolSend[0], snode.SstrUcolSend.size() );
          snode.SizeSstrUcolSend = snode.SstrUcolSend.size();
          TIMER_STOP(Serialize_UcolU);

          for( Int iProcCol = 0; 
              iProcCol < this->grid_->numProcCol ; iProcCol++ ){
            if( MYCOL( this->grid_ ) != iProcCol &&
                this->isSendToRight_( iProcCol, snode.Index ) == true ){
              // Use Isend to send to multiple targets
              MPI_Isend( &snode.SizeSstrUcolSend, 1, MPI_INT,  
                  iProcCol, IDX_TO_TAG2(snode.Index,supidx,SELINV_TAG_UCOL_SIZE), 
                  this->grid_->rowComm, &mpireqsSendUToRight[2*iProcCol] );
              MPI_Isend( &snode.SstrUcolSend[0],snode.SizeSstrUcolSend,
                  MPI_BYTE, iProcCol, 
                  IDX_TO_TAG2(snode.Index,supidx,SELINV_TAG_UCOL_CONTENT), 
                  this->grid_->rowComm, 
                  &mpireqsSendUToRight[2*iProcCol+1] );
#if ( _DEBUGlevel_ >= 2 )
              statusOFS<< "["<<snode.Index<<"] "<<"Ucol SENT "<<std::endl;
              for(Int ib=0;ib<Ucol.size();++ib){
                statusOFS<<Ucol[ib].blockIdx<<" ";
              }
              statusOFS<<std::endl;
#endif
#if ( _DEBUGlevel_ >= 1 )
              statusOFS << "["<<snode.Index<<"] "<<  "Sending Ucol "
                << snode.SizeSstrUcolSend << " BYTES to P" 
                << PNUM(MYROW(this->grid_),iProcCol,this->grid_)
                << std::endl <<  std::endl; 
#endif
            } // Send 
          } // for (iProcCol)
        } // if I am the sender



      } //Senders

      //TODO Ideally, we should not receive data in sequence 
      // but in any order with ksup packed with the data

      TIMER_START(WaitContentLU);
      // Receivers (Size)
      for (Int supidx=0; supidx<stepSuper ; supidx++){
        SuperNodeBufferTypeUnsym & snode = arrSuperNodes[supidx];
        MPI_Request * mpireqsRecvLFromAbove = 
          &arrMpireqsRecvLSizeFromAny[supidx*2];
        MPI_Request * mpireqsRecvLFromLeft = 
          &arrMpireqsRecvLSizeFromAny[supidx*2+1];
        MPI_Request * mpireqsRecvUFromAbove = 
          &arrMpireqsRecvUSizeFromAny[supidx*2];
        MPI_Request * mpireqsRecvUFromLeft = 
          &arrMpireqsRecvUSizeFromAny[supidx*2+1];

        // Receive the size first
        if( this->isRecvFromAbove_( snode.Index ) && 
            MYROW( this->grid_ ) != PROW( snode.Index, this->grid_ ) ){
          Int sender = PROW( snode.Index, this->grid_ ); 
          MPI_Irecv( &snode.SizeSstrLrowRecv, 1, MPI_INT, sender, 
              IDX_TO_TAG2(snode.Index,supidx,SELINV_TAG_LROW_SIZE),
              this->grid_->colComm, mpireqsRecvLFromAbove );
#if ( _DEBUGlevel_ >= 1 )
          statusOFS << "["<<snode.Index<<"] "<<  "Receiving Lrow"
            <<" size on tag "<<IDX_TO_TAG2(snode.Index,supidx,SELINV_TAG_LROW_SIZE)
            <<" from P" << PNUM(sender,MYCOL(this->grid_),this->grid_)
            << std::endl << std::endl; 
#endif
          MPI_Irecv( &snode.SizeSstrUrowRecv, 1, MPI_INT, sender, 
              IDX_TO_TAG2(snode.Index,supidx,SELINV_TAG_U_SIZE),
              this->grid_->colComm, mpireqsRecvUFromAbove );
#if ( _DEBUGlevel_ >= 1 )
          statusOFS << "["<<snode.Index<<"] "<<  "Receiving U"
            <<" size on tag "<<IDX_TO_TAG2(snode.Index,supidx,SELINV_TAG_U_SIZE)
            <<" from P" << PNUM(sender,MYCOL(this->grid_),this->grid_)
            << std::endl << std::endl; 
#endif
        } // if I need to receive from up


        if( this->isRecvFromLeft_( snode.Index ) &&
            MYCOL( this->grid_ ) != PCOL( snode.Index, this->grid_ ) ){
          Int sender = PCOL( snode.Index, this->grid_ ); 
          MPI_Irecv( &snode.SizeSstrLcolRecv, 1, MPI_INT, sender, 
              IDX_TO_TAG2(snode.Index,supidx,SELINV_TAG_L_SIZE),
              this->grid_->rowComm, mpireqsRecvLFromLeft );
#if ( _DEBUGlevel_ >= 1 )
          statusOFS << "["<<snode.Index<<"] "<<  "Receiving L"
            <<" size on tag "<<IDX_TO_TAG2(snode.Index,supidx,SELINV_TAG_L_SIZE)
            <<" from P" << PNUM(MYROW(this->grid_),sender,this->grid_)
            << std::endl << std::endl; 
#endif
          MPI_Irecv( &snode.SizeSstrUcolRecv, 1, MPI_INT, sender, 
              IDX_TO_TAG2(snode.Index,supidx,SELINV_TAG_UCOL_SIZE),
              this->grid_->rowComm, mpireqsRecvUFromLeft );
#if ( _DEBUGlevel_ >= 1 )
          statusOFS << "["<<snode.Index<<"] "<<  "Receiving Ucol"
            <<" size on tag "<<IDX_TO_TAG2(snode.Index,supidx,SELINV_TAG_UCOL_SIZE)
            <<" from P" << PNUM(MYROW(this->grid_),sender,this->grid_)
            << std::endl << std::endl; 
#endif
        } // if I need to receive from left
      }
      TIMER_STOP(WaitContentLU);

      //Wait to receive all the sizes for L
      TIMER_START(WaitContentLU);
      TIMER_START(WaitSize_LrowL);
      mpi::Waitall(arrMpireqsRecvLSizeFromAny);
      TIMER_STOP(WaitSize_LrowL);

      // Receivers (Content)
      for (Int supidx=0; supidx<stepSuper ; supidx++){
        SuperNodeBufferTypeUnsym & snode = arrSuperNodes[supidx];

        MPI_Request * mpireqsRecvFromAbove =
          &arrMpireqsRecvLContentFromAny[supidx*2];
        MPI_Request * mpireqsRecvFromLeft = 
          &arrMpireqsRecvLContentFromAny[supidx*2+1];

        TIMER_START(Alloc_Buffer_Recv_LrowL);
        if( this->isRecvFromAbove_( snode.Index ) && 
            MYROW( this->grid_ ) != PROW( snode.Index, this->grid_ ) ){
          snode.SstrLrowRecv.resize( snode.SizeSstrLrowRecv );
          Int sender = PROW( snode.Index, this->grid_ ); 
          MPI_Irecv( &snode.SstrLrowRecv[0], snode.SizeSstrLrowRecv, MPI_BYTE, 
              sender, IDX_TO_TAG2(snode.Index,supidx,SELINV_TAG_LROW_CONTENT), 
              this->grid_->colComm, mpireqsRecvFromAbove );
#if ( _DEBUGlevel_ >= 1 )
          statusOFS << "["<<snode.Index<<"] "<<  "Receiving Lrow "
            << snode.SizeSstrLrowRecv << " BYTES from P" 
            << PNUM(sender,MYCOL(this->grid_),this->grid_)
            << std::endl <<  std::endl; 
#endif
        } // if I need to receive from up
        TIMER_STOP(Alloc_Buffer_Recv_LrowL);

        TIMER_START(Alloc_Buffer_Recv_LcolL);
        if( this->isRecvFromLeft_( snode.Index ) &&
            MYCOL( this->grid_ ) != PCOL( snode.Index, this->grid_ ) ){
          snode.SstrLcolRecv.resize( snode.SizeSstrLcolRecv );
          Int sender = PCOL( snode.Index, this->grid_ ); 
          MPI_Irecv( &snode.SstrLcolRecv[0], snode.SizeSstrLcolRecv, MPI_BYTE, 
              sender, IDX_TO_TAG2(snode.Index,supidx,SELINV_TAG_L_CONTENT), 
              this->grid_->rowComm,
              mpireqsRecvFromLeft );
#if ( _DEBUGlevel_ >= 1 )
          statusOFS << "["<<snode.Index<<"] "<<  "Receiving L "
            << snode.SizeSstrLcolRecv << " BYTES from P" 
            << PNUM(MYROW(this->grid_),sender,this->grid_)
            << std::endl <<  std::endl; 
#endif
        } // if I need to receive from left
        TIMER_STOP(Alloc_Buffer_Recv_LcolL);
      }
      TIMER_STOP(WaitContentLU);


      //Wait to receive all the sizes for U
      TIMER_START(WaitContentLU);
      TIMER_START(WaitSize_UcolU);
      mpi::Waitall(arrMpireqsRecvUSizeFromAny);
      TIMER_STOP(WaitSize_UcolU);

      // Receivers (Content)
      for (Int supidx=0; supidx<stepSuper ; supidx++){
        SuperNodeBufferTypeUnsym & snode = arrSuperNodes[supidx];

        MPI_Request * mpireqsRecvFromAbove = 
          &arrMpireqsRecvUContentFromAny[supidx*2];
        MPI_Request * mpireqsRecvFromLeft = 
          &arrMpireqsRecvUContentFromAny[supidx*2+1];

        TIMER_START(Alloc_Buffer_Recv_UrowL);
        if( this->isRecvFromAbove_( snode.Index ) && 
            MYROW( this->grid_ ) != PROW( snode.Index, this->grid_ ) ){
          snode.SstrUrowRecv.resize( snode.SizeSstrUrowRecv );
          Int sender = PROW( snode.Index, this->grid_ ); 
          MPI_Irecv( &snode.SstrUrowRecv[0], snode.SizeSstrUrowRecv, MPI_BYTE,
              sender, IDX_TO_TAG2(snode.Index,supidx,SELINV_TAG_U_CONTENT), 
              this->grid_->colComm, mpireqsRecvFromAbove );
#if ( _DEBUGlevel_ >= 1 )
          statusOFS << "["<<snode.Index<<"] "<<  "Receiving U "
            << snode.SizeSstrUrowRecv << " BYTES from P" 
            << PNUM(sender,MYCOL(this->grid_),this->grid_)
            << std::endl <<  std::endl; 
#endif
        } // if I need to receive from up
        TIMER_STOP(Alloc_Buffer_Recv_UrowL);

        TIMER_START(Alloc_Buffer_Recv_UcolL);
        if( this->isRecvFromLeft_( snode.Index ) &&
            MYCOL( this->grid_ ) != PCOL( snode.Index, this->grid_ ) ){
          snode.SstrUcolRecv.resize( snode.SizeSstrUcolRecv );
          Int sender = PCOL( snode.Index, this->grid_ ); 
          MPI_Irecv( &snode.SstrUcolRecv[0], snode.SizeSstrUcolRecv, MPI_BYTE, 
              sender, IDX_TO_TAG2(snode.Index,supidx,SELINV_TAG_UCOL_CONTENT), 
              this->grid_->rowComm,
              mpireqsRecvFromLeft );
#if ( _DEBUGlevel_ >= 1 )
          statusOFS << "["<<snode.Index<<"] "<<  "Receiving Ucol "
            << snode.SizeSstrUcolRecv << " BYTES from P" 
            << PNUM(MYROW(this->grid_),sender,this->grid_)
            << std::endl <<  std::endl; 
#endif
        } // if I need to receive from left
        TIMER_STOP(Alloc_Buffer_Recv_UcolL);
      }

      TIMER_STOP(WaitContentLU);


      TIMER_START(Compute_Sinv_LU);
      {
        Int gemmProcessed = 0;
        Int gemmToDo = 0;
        //      Int toRecvGemm = 0;
        //copy the list of supernodes we need to process
        std::vector<Int> readySupidx;
        //find local things to do
        for(Int supidx = 0;supidx<stepSuper;supidx++){
          SuperNodeBufferTypeUnsym & snode = arrSuperNodes[supidx];
          if( this->isRecvFromAbove_( snode.Index ) 
              && this->isRecvFromLeft_( snode.Index )){
            gemmToDo+=2;
            if( MYCOL( this->grid_ ) == PCOL( snode.Index, this->grid_ ) ){
              snode.isReady+=2;
            }

            if(  MYROW( this->grid_ ) == PROW( snode.Index, this->grid_ ) ){
              snode.isReady+=2;
            }

            if(snode.isReady==4){
              readySupidx.push_back(supidx);
#if ( _DEBUGlevel_ >= 1 )
              statusOFS<<"Locally processing ["<<snode.Index<<"]"<<std::endl;
#endif
            }
          }
          else{
            TIMER_START(Reduce_Sinv_L_Send);
            if( this->isRecvFromLeft_( snode.Index )  
                && MYCOL( this->grid_ ) != PCOL( snode.Index, this->grid_ ) ){
              MPI_Request & mpireqsSendToLeft = arrMpireqsSendLToLeft[supidx];
              // Dummy 0-b send If I was a receiver, I need to send my data to
              // proc in column of snode.Index
              MPI_Isend( NULL, 0, MPI_BYTE, PCOL(snode.Index,this->grid_),
                  IDX_TO_TAG2(snode.Index,supidx,SELINV_TAG_L_REDUCE), 
                  this->grid_->rowComm, &mpireqsSendToLeft );

#if ( _DEBUGlevel_ >= 1 )
              Int dst = PNUM(MYROW(this->grid_),
                  PCOL(snode.Index,this->grid_),this->grid_);
              statusOFS << "["<<snode.Index<<"] "<< " LReduce P"
                << MYPROC(this->grid_) << " has sent "
                << 0 << " bytes to "
                << dst << std::endl;
#endif
            }// if( isRecvFromLeft_( snode.Index ))
            TIMER_STOP(Reduce_Sinv_L_Send);



          }
        }

#if ( _DEBUGlevel_ >= 1 )
        statusOFS<<std::endl<<"gemmToDo ="<<gemmToDo<<std::endl;
#endif


#if defined (PROFILE) 
        end_SendULWaitContentFirst=0;
        begin_SendULWaitContentFirst=0;
#endif

        while(gemmProcessed<gemmToDo)
        {
          Int reqidx = MPI_UNDEFINED;
          Int supidx = -1;

          TIMER_START(WaitContentLU);
          //while I don't have anything to do, wait for data to arrive 
          do{
            int reqIndicesL[arrMpireqsRecvLContentFromAny.size()];
            int numRecv = 0; 

            //then process with the remote ones

#if defined(PROFILE)
            if(begin_SendULWaitContentFirst==0){
              begin_SendULWaitContentFirst=1;
              TIMER_START(WaitContent_LrowL_First);
            }
#endif

            TIMER_START(WaitContent_LrowL);
            MPI_Waitsome(2*stepSuper, &arrMpireqsRecvLContentFromAny[0],
                &numRecv, reqIndicesL, MPI_STATUSES_IGNORE);

            for(int i =0;i<numRecv;i++){
              reqidx = reqIndicesL[i];
              //I've received something
              if(reqidx!=MPI_UNDEFINED){ 

                supidx = reqidx/2;
                SuperNodeBufferTypeUnsym & snode = arrSuperNodes[supidx];
                snode.isReady++;

#if ( _DEBUGlevel_ >= 1 )
                statusOFS << "["<<snode.Index<<"] "<<"Received data for L"
                  << " reqidx%2=" << reqidx%2 
                  << " is ready ?"<<snode.isReady<<std::endl;
#endif
                //if we received both L and U, the supernode is ready
                if(snode.isReady==4){
                  readySupidx.push_back(supidx);
#if defined(PROFILE)
                  if(end_SendULWaitContentFirst==0){
                    TIMER_STOP(WaitContent_LrowL_First);
                    end_SendULWaitContentFirst=1;
                  }
#endif
                }
              }

            }//end for waitsome
            TIMER_STOP(WaitContent_LrowL);


            int reqIndicesU[arrMpireqsRecvUContentFromAny.size()];
            numRecv = 0; 


            TIMER_START(WaitContent_UcolU);
            MPI_Waitsome(2*stepSuper, &arrMpireqsRecvUContentFromAny[0],
                &numRecv, reqIndicesU, MPI_STATUSES_IGNORE);

            for(int i =0;i<numRecv;i++){
              reqidx = reqIndicesU[i];
              //I've received something
              if(reqidx!=MPI_UNDEFINED){ 
                supidx = reqidx/2;
                SuperNodeBufferTypeUnsym & snode = arrSuperNodes[supidx];
                snode.isReady++;

#if ( _DEBUGlevel_ >= 1 )
                statusOFS << "["<<snode.Index<<"] "<<"Received data for U"
                  << " reqidx%2=" << reqidx%2 
                  << " is ready ?"<<snode.isReady<<std::endl;
#endif
                //if we received both L and U, the supernode is ready
                if(snode.isReady==4){
                  readySupidx.push_back(supidx);
#if defined(PROFILE)
                  if(end_SendULWaitContentFirst==0){
                    TIMER_STOP(WaitContent_LrowL_First);
                    end_SendULWaitContentFirst=1;
                  }
#endif
                }
              }

            }//end for waitsome
            TIMER_STOP(WaitContent_UcolU);

          } while(readySupidx.size()==0);
          TIMER_STOP(WaitContentLU);

          //If I have some work to do 
          if(readySupidx.size()>0)
          {
            supidx = readySupidx.back();
            readySupidx.pop_back();
            SuperNodeBufferTypeUnsym & snode = arrSuperNodes[supidx];


            // Only the processors received information participate in the Gemm 
            if( this->isRecvFromAbove_( snode.Index )
                && this->isRecvFromLeft_( snode.Index ) ){
              std::vector<LBlock<T> > LcolRecv;
              std::vector<LBlock<T> > LrowRecv;
              std::vector<UBlock<T> > UcolRecv;
              std::vector<UBlock<T> > UrowRecv;

              //TODO REMOVE THIS THIS IS ONLY FOR DEBUGING PURPOSE
              NumMat<T> * pAinvBuf = new NumMat<T>();
              NumMat<T> * pUBuf = new NumMat<T>();
              NumMat<T> * pLBuf = new NumMat<T>();
              NumMat<T> & AinvBuf = *pAinvBuf;
              NumMat<T> & UBuf = *pUBuf;
              NumMat<T> & LBuf = *pLBuf;

              UnpackData(snode, LcolRecv, LrowRecv, UcolRecv, UrowRecv);

              SelInv_lookup_indexes(snode, LcolRecv, LrowRecv, 
                  UcolRecv, UrowRecv,AinvBuf,LBuf, UBuf);


#if ( _DEBUGlevel_ >= 2 )
              statusOFS << "["<<snode.Index<<"] " << "LBuf: ";
              statusOFS << LBuf << std::endl;
              statusOFS << "["<<snode.Index<<"] " << "UBuf: ";
              statusOFS << UBuf << std::endl;
#endif

              NumMat<T> LUpdateBuf;
              LUpdateBuf.Resize( AinvBuf.m(), 
                  SuperSize( snode.Index, this->super_ ) );



              TIMER_START(Compute_Sinv_L_Resize);
              snode.LUpdateBuf.Resize( AinvBuf.m(), 
                  SuperSize( snode.Index, this->super_ ) );
              TIMER_STOP(Compute_Sinv_L_Resize);

              TIMER_START(Compute_Sinv_LT_GEMM);
              blas::Gemm('N', 'T', AinvBuf.m(), LBuf.m(), AinvBuf.n(),
                  MINUS_ONE<T>(), AinvBuf.Data(), AinvBuf.m(), 
                  LBuf.Data(), LBuf.m(), ZERO<T>(), 
                  snode.LUpdateBuf.Data(), snode.LUpdateBuf.m() );
              TIMER_STOP(Compute_Sinv_LT_GEMM);

              TIMER_START(Compute_Sinv_U_Resize);
              snode.UUpdateBuf.Resize( SuperSize( snode.Index, this->super_ ),
                  AinvBuf.n() );
              TIMER_STOP(Compute_Sinv_U_Resize);

              TIMER_START(Compute_Sinv_U_GEMM);
              blas::Gemm('N', 'N', UBuf.m(), AinvBuf.n(), AinvBuf.m(),
                  MINUS_ONE<T>(), UBuf.Data(), UBuf.m(),
                  AinvBuf.Data(), AinvBuf.m(), ZERO<T>(),
                  snode.UUpdateBuf.Data(), snode.UUpdateBuf.m() );
              TIMER_STOP(Compute_Sinv_U_GEMM);

#if ( _DEBUGlevel_ >= 2 )
              statusOFS << "["<<snode.Index<<"] " << "snode.LUpdateBuf: ";
              statusOFS << snode.LUpdateBuf << std::endl;
              statusOFS << "["<<snode.Index<<"] " << "snode.UUpdateBuf: ";
              statusOFS << snode.UUpdateBuf << std::endl;
#endif

              //TODO REMOVE THIS THIS IS ONLY FOR DEBUGING PURPOSE
              delete pLBuf;
              delete pUBuf;
              delete pAinvBuf;


            } // if Gemm is to be done locally


            TIMER_START(Reduce_Sinv_L_Send);
            // If I was a receiver, I need to send my data to proc in column
            // of snode.Index
            if( this->isRecvFromAbove_( snode.Index )  ){
              if( this->isRecvFromLeft_( snode.Index ) 
                  && MYCOL( this->grid_ ) != PCOL( snode.Index, this->grid_ ) ){
                MPI_Request & reqsSendToLeft = arrMpireqsSendLToLeft[supidx];
                MPI_Isend( snode.LUpdateBuf.Data(), snode.LUpdateBuf.ByteSize(),
                    MPI_BYTE, PCOL(snode.Index,this->grid_),
                    IDX_TO_TAG2(snode.Index,supidx,SELINV_TAG_L_REDUCE), 
                    this->grid_->rowComm, &reqsSendToLeft );

#if ( _DEBUGlevel_ >= 1 )
                Int dst = PNUM(MYROW(this->grid_),
                    PCOL(snode.Index,this->grid_),this->grid_);
                statusOFS << "["<<snode.Index<<"] "<< " LReduce P"
                  << MYPROC(this->grid_) << " has sent "
                  << snode.LUpdateBuf.ByteSize() << " bytes to "
                  << dst << std::endl;
#endif
              }//Sender
            }
            TIMER_STOP(Reduce_Sinv_L_Send);

            gemmProcessed+=2;
#if ( _DEBUGlevel_ >= 1 )
            statusOFS << "["<<snode.Index<<"] "<<"gemmProcessed = "
              << gemmProcessed<<"/"<<gemmToDo<<std::endl;
#endif
          }
        }

      }
      TIMER_STOP(Compute_Sinv_LU);







      //Reduce Sinv L to the processors in PCOL(ksup,this->grid_)
      TIMER_START(Reduce_Sinv_L);

      for (Int supidx=0; supidx<stepSuper; supidx++){
        SuperNodeBufferTypeUnsym & snode = arrSuperNodes[supidx];
        if( MYCOL( this->grid_ ) == PCOL( snode.Index, this->grid_ ) ){
          //determine the number of rows in LUpdateBufReduced
          Int numRowLUpdateBuf;
          std::vector<LBlock<T> >&  Lcol = 
            this->L( LBj( snode.Index, this->grid_ ) );

          // If I own the diagonal block, skip the diagonal block
          Int offset =
            (MYROW( this->grid_ )==PROW(snode.Index,this->grid_))?1:0;
          snode.RowLocalPtr.resize( Lcol.size() + 1 - offset );
          snode.BlockIdxLocal.resize( Lcol.size() - offset );
          snode.RowLocalPtr[0] = 0;
          for( Int ib = 0; ib < snode.BlockIdxLocal.size(); ib++ ){
            snode.RowLocalPtr[ib+1] = snode.RowLocalPtr[ib]
              + Lcol[ib+offset].numRow;
            snode.BlockIdxLocal[ib] = Lcol[ib+offset].blockIdx;
          }
          numRowLUpdateBuf = *snode.RowLocalPtr.rbegin();


          if( numRowLUpdateBuf > 0 ){
            if( snode.LUpdateBuf.m() == 0 && snode.LUpdateBuf.n() == 0 ){
              snode.LUpdateBuf.Resize( numRowLUpdateBuf,
                  SuperSize( snode.Index, this->super_ ) );
              // Fill zero is important
              SetValue( snode.LUpdateBuf, ZERO<T>() );
            }
          }

#if ( _DEBUGlevel_ >= 2 )
          statusOFS << "["<<snode.Index<<"] "<<"LUpdateBuf Before Reduction: "
            <<  snode.LUpdateBuf << std::endl << std::endl; 
#endif

          Int totCountRecv = 0;
          Int numRecv = this->CountSendToRight(snode.Index);
          NumMat<T>  LUpdateBufRecv(numRowLUpdateBuf,
              SuperSize( snode.Index, this->super_ ) );
          for( Int countRecv = 0; countRecv < numRecv ; ++countRecv ){
            //Do the blocking recv
            MPI_Status stat;
            Int size = 0;
            TIMER_START(L_RECV);
            MPI_Recv(LUpdateBufRecv.Data(), LUpdateBufRecv.ByteSize(), 
                MPI_BYTE, MPI_ANY_SOURCE,
                IDX_TO_TAG2(snode.Index,supidx,SELINV_TAG_L_REDUCE),
                this->grid_->rowComm,&stat);
            TIMER_STOP(L_RECV);
            MPI_Get_count(&stat, MPI_BYTE, &size);
            //if the processor contributes
            if(size>0){
#if ( _DEBUGlevel_ >= 1 )
              Int src = PNUM(MYROW(this->grid_),stat.MPI_SOURCE,this->grid_);
              statusOFS << "["<<snode.Index<<"] "<< " LReduce P"
                << MYPROC(this->grid_)<<" has received "<< size 
                << " bytes from " << src << std::endl;
#endif
#if ( _DEBUGlevel_ >= 2 )
              statusOFS << "["<<snode.Index<<"] "<<   "LUpdateBufRecv: "
                << LUpdateBufRecv << std::endl << std::endl; 
#endif
              //do the sum
              blas::Axpy(snode.LUpdateBuf.Size(), ONE<T>(), 
                  LUpdateBufRecv.Data(), 1, 
                  snode.LUpdateBuf.Data(), 1 );
            }
          } // for (iProcCol)
#if ( _DEBUGlevel_ >= 2 ) 
          statusOFS << "["<<snode.Index<<"] "<<   "LUpdateBuf After Reduction: "
            <<  snode.LUpdateBuf << std::endl << std::endl; 
#endif
        } // Receiver
      }

      TIMER_STOP(Reduce_Sinv_L);
      //--------------------- End of reduce of LUpdateBuf-------------------------


      mpi::Waitall( arrMpireqsSendLToLeft );


      if(next_lidx < superList.size()){
        IRecvContentCD(superList[next_lidx], superList[next_lidx].size(),nextCDBuffers);
      }





      TIMER_START(Update_Diagonal);
      for (Int supidx=0; supidx<stepSuper; supidx++){
        SuperNodeBufferTypeUnsym & snode = arrSuperNodes[supidx];

        ComputeDiagUpdate(snode);

        if( MYCOL( this->grid_ ) == PCOL( snode.Index, this->grid_ ) ){
          if( MYROW( this->grid_ ) != PROW( snode.Index, this->grid_ ) ){
            if(this->isSendToDiagonal_(snode.Index)){
              //send to above
              MPI_Request & mpireqsSendToAbove = arrMpireqsSendToAbove[supidx];
              MPI_Isend( snode.DiagBuf.Data(), snode.DiagBuf.ByteSize(), 
                  MPI_BYTE, PROW(snode.Index,this->grid_) ,
                  IDX_TO_TAG2(snode.Index,supidx,SELINV_TAG_D_REDUCE), 
                  this->grid_->colComm, &mpireqsSendToAbove );

#if ( _DEBUGlevel_ >= 1 )
              statusOFS << "["<<snode.Index<<"] "<< " P"<<MYROW(this->grid_)
                <<" has sent "<< snode.DiagBuf.ByteSize() 
                << " bytes of DiagBuf to " 
                << PROW(snode.Index,this->grid_) << std::endl;
#endif
            }
          }
        }
      }

      TIMER_STOP(Update_Diagonal);



      TIMER_START(Reduce_Diagonal);

      for (Int supidx=0; supidx<stepSuper; supidx++){
        SuperNodeBufferTypeUnsym & snode = arrSuperNodes[supidx];
        if( MYCOL( this->grid_ ) == PCOL( snode.Index, this->grid_ ) ){
          if( MYROW( this->grid_ ) == PROW( snode.Index, this->grid_ ) ){
            if(snode.DiagBuf.Size()==0){
              snode.DiagBuf.Resize( SuperSize( snode.Index, this->super_ ),
                  SuperSize( snode.Index, this->super_ ));
              SetValue(snode.DiagBuf, ZERO<T>());
            }
            //receive from below
            Int totCountRecv = 0;
            Int numRecv = this->CountRecvFromBelow(snode.Index);
            NumMat<T>  DiagBufRecv(snode.DiagBuf.m(),snode.DiagBuf.n());

            for( Int countRecv = 0; countRecv < numRecv ; ++countRecv ){
              //Do the blocking recv
              MPI_Status stat;
              Int size = 0;
              TIMER_START(D_RECV);
              MPI_Recv(DiagBufRecv.Data(), DiagBufRecv.ByteSize(), MPI_BYTE, 
                  MPI_ANY_SOURCE,IDX_TO_TAG2(snode.Index,supidx,SELINV_TAG_D_REDUCE),
                  this->grid_->colComm,&stat);
              TIMER_STOP(D_RECV);
              MPI_Get_count(&stat, MPI_BYTE, &size);
              //if the processor contributes
              if(size>0){
                // Add DiagBufRecv to diagonal block.
                blas::Axpy(snode.DiagBuf.Size(), ONE<T>(),
                    DiagBufRecv.Data(), 1,
                    snode.DiagBuf.Data(), 1 );
              }
            }
            LBlock<T> &  LB = this->L( LBj( snode.Index, this->grid_ ) )[0];
            blas::Axpy( LB.numRow * LB.numCol, ONE<T>(), 
                snode.DiagBuf.Data(), 1,
                LB.nzval.Data(), 1 );
          }

        } 
      }


      TIMER_STOP(Reduce_Diagonal);



      //Reduce U Sinv  to the processors in PROW(ksup,this->grid_)
      TIMER_START(Reduce_Sinv_U);

      for (Int supidx=0; supidx<stepSuper; supidx++){
        SuperNodeBufferTypeUnsym & snode = arrSuperNodes[supidx];


        if( MYROW( this->grid_ ) == PROW( snode.Index, this->grid_ ) ){
          //determine the number of rows in UUpdateBufReduced
          Int numColUUpdateBuf;
          //FIXME U must be revised to store the same structure as L ?
          std::vector<UBlock<T> >&  Urow = 
            this->U( LBi( snode.Index, this->grid_ ) );

          snode.ColLocalPtr.resize( Urow.size() + 1 );
          snode.BlockIdxLocalU.resize( Urow.size() );
          snode.ColLocalPtr[0] = 0;

          //            std::vector<LBlock<T> >&  Lrow = this->Lrow( LBi( snode.Index, this->grid_ ) );
          //            if(Lrow.size()>=Urow.size()){
          //              statusOFS<<"UReduce first case"<<std::endl;
          //              std::vector<Int> colPtrL(Lrow.size()+1);
          //              colPtrL[0] = 0;
          //              for( Int ib = 0; ib < Lrow.size(); ib++ ){
          //                colPtrL[ib+1] = colPtrL[ib] + Lrow[ib].numCol;
          //              }
          //
          //
          //
          //              for( Int jb = 0; jb < Urow.size(); jb++ ){
          //                Int indexL =0;
          //                for( Int ib = 0; ib < Lrow.size(); ib++ ){
          //                  if(Lrow[ib].blockIdx == Urow[jb].blockIdx){
          //                    indexL = ib;
          //                    break;
          //                  }
          //                }
          //                statusOFS<<jb<<" vs "<<indexL<<std::endl;
          //
          //                snode.ColLocalPtr[jb] = colPtrL[indexL];
          //                //snode.ColLocalPtr[jb+1] = snode.ColLocalPtr[jb] + Urow[jb].numCol;
          //                snode.BlockIdxLocalU[jb] = Lrow[indexL].blockIdx;
          //                //snode.BlockIdxLocalU[jb] = Urow[jb].blockIdx;
          //              }
          //              snode.ColLocalPtr.back()=colPtrL.back();
          //
          //              statusOFS<<colPtrL<<std::endl;
          //            }
          //            else{
          //              statusOFS<<"UReduce second case"<<std::endl;
          Int urowsize = Urow.size();
          snode.ColLocalPtr[0] = 0;
          for( Int jb = 0; jb < Urow.size(); jb++ ){
            UBlock<T> & UB = Urow[jb];
            snode.ColLocalPtr[jb+1] = snode.ColLocalPtr[jb]+UB.numCol;
            snode.BlockIdxLocalU[jb] = UB.blockIdx;
          }
          //            }

#if ( _DEBUGlevel_ >= 2 )
          statusOFS << "["<<snode.Index<<"] "<<"UReduce collocalptr "
            << snode.ColLocalPtr << std::endl;
          statusOFS << "["<<snode.Index<<"] "<<"UReduce blockidxlocalU "
            << snode.BlockIdxLocalU << std::endl;
#endif

          numColUUpdateBuf = *snode.ColLocalPtr.rbegin();

          if( numColUUpdateBuf > 0 ){
            if( snode.UUpdateBuf.m() == 0 && snode.UUpdateBuf.n() == 0 ){
              snode.UUpdateBuf.Resize( SuperSize( snode.Index, this->super_ ),
                  numColUUpdateBuf );
              // Fill zero is important
              SetValue( snode.UUpdateBuf, ZERO<T>() );
            }
          }

#if ( _DEBUGlevel_ >= 2 )
          statusOFS << "["<<snode.Index<<"] "<<"UUpdateBuf Before Reduction: "
            << snode.UUpdateBuf << std::endl << std::endl;
#endif

          Int totCountRecv = 0;

          Int numRecv = this->CountSendToBelow(snode.Index);

          NumMat<T>  UUpdateBufRecv(SuperSize( snode.Index, this->super_ ),
              numColUUpdateBuf );

#if ( _DEBUGlevel_ >= 1 )
          statusOFS << "["<<snode.Index<<"] "<< " UReduce P"
            << MYPROC(this->grid_) << " can receive at most "
            << UUpdateBufRecv.ByteSize() << " bytes" << std::endl;
#endif
          for( Int countRecv = 0; countRecv < numRecv ; ++countRecv ){
            //Do the blocking recv
            MPI_Status stat;
            Int size = 0;
            TIMER_START(U_RECV);
            MPI_Recv(UUpdateBufRecv.Data(), UUpdateBufRecv.ByteSize(),
                MPI_BYTE, MPI_ANY_SOURCE,
                IDX_TO_TAG2(snode.Index,supidx,SELINV_TAG_U_REDUCE),
                this->grid_->colComm,&stat);
            TIMER_STOP(U_RECV);
            MPI_Get_count(&stat, MPI_BYTE, &size);
            //if the processor contributes
            if(size>0){

#if ( _DEBUGlevel_ >= 1 )
              Int src = PNUM(stat.MPI_SOURCE,MYCOL(this->grid_),this->grid_);
              statusOFS << "["<<snode.Index<<"] "<< " UReduce P"
                << MYPROC(this->grid_)<<" has received "
                << size << " bytes from P" << src << std::endl;
#endif

#if ( _DEBUGlevel_ >= 2 )
              statusOFS << "["<<snode.Index<<"] "<< "UUpdateBufRecv: "
                <<  UUpdateBufRecv << std::endl << std::endl; 
#endif
              //do the sum
              blas::Axpy(snode.UUpdateBuf.Size(), ONE<T>(), 
                  UUpdateBufRecv.Data(), 1,
                  snode.UUpdateBuf.Data(), 1);
            }
          } // for (iProcRow)
#if ( _DEBUGlevel_ >= 2 ) 
          statusOFS << "["<<snode.Index<<"] "<<"UUpdateBuf After Reduction: "
            <<  snode.UUpdateBuf << std::endl << std::endl; 
#endif
        } // Receiver
        else{


          TIMER_START(Reduce_Sinv_U_Send);
          if(!this->isRecvFromLeft_( snode.Index ) && this->isRecvFromAbove_( snode.Index )   ){
            MPI_Request & mpireqsSendToAbove = arrMpireqsSendUToAbove[supidx];
            // Dummy 0-b send If I was a receiver, I need to send my data to 
            // proc in row of snode.Index
            MPI_Isend( NULL, 0, MPI_BYTE, PROW(snode.Index,this->grid_),
                IDX_TO_TAG2(snode.Index,supidx,SELINV_TAG_U_REDUCE), 
                this->grid_->colComm, &mpireqsSendToAbove );

#if ( _DEBUGlevel_ >= 1 )
            Int dst = PNUM(PROW(snode.Index,this->grid_),
                MYCOL(this->grid_),this->grid_);
            statusOFS << "["<<snode.Index<<"] "<< " UReduce P"
              << MYPROC(this->grid_) << " has sent "
              << 0 << " bytes to "
              << dst << std::endl;
#endif
          }// if( isRecvFromAbove_( snode.Index ) )
          TIMER_STOP(Reduce_Sinv_U_Send);


          TIMER_START(Reduce_Sinv_U_Send);
          // If I was a receiver, I need to send my data to proc in row 
          // of snode.Index
          if( this->isRecvFromAbove_( snode.Index )  && this->isRecvFromLeft_( snode.Index )  ){
            MPI_Request & reqsSendToAbove = arrMpireqsSendUToAbove[supidx];
            TIMER_START(Reduce_Sinv_U_Send_Isend);
            MPI_Isend( snode.UUpdateBuf.Data(), snode.UUpdateBuf.ByteSize(),
                MPI_BYTE, PROW(snode.Index,this->grid_),
                IDX_TO_TAG2(snode.Index,supidx,SELINV_TAG_U_REDUCE), 
                this->grid_->colComm, &reqsSendToAbove );
            TIMER_STOP(Reduce_Sinv_U_Send_Isend);
#if ( _DEBUGlevel_ >= 1 )
            Int dst = PNUM(PROW(snode.Index,this->grid_),
                MYCOL(this->grid_),this->grid_);
            statusOFS << "["<<snode.Index<<"] "<< " UReduce P"
              << MYPROC(this->grid_) << " has sent "
              << snode.UUpdateBuf.ByteSize() << " bytes to "
              << dst << std::endl;
#endif
          }//Sender
          TIMER_STOP(Reduce_Sinv_U_Send);


        }
      }

      TIMER_STOP(Reduce_Sinv_U);
      //--------------------- End of reduce of UUpdateBuf-------------------------

      mpi::Waitall( arrMpireqsSendUToAbove );






      //Deallocate Lrow and Ucol
      for (Int supidx=0; supidx<stepSuper; supidx++){ 
        SuperNodeBufferTypeUnsym & snode = arrSuperNodes[supidx];
        if( MYROW( this->grid_ ) == PROW( snode.Index, this->grid_ )){ 
          std::vector<LBlock<T> >&  Lrow = this->Lrow( LBi(snode.Index, this->grid_) );
          Lrow.clear();
        }
        if( MYCOL( this->grid_ ) == PCOL( snode.Index, this->grid_ )){ 
          std::vector<UBlock<T> >&  Ucol = this->Ucol( LBj(snode.Index, this->grid_) );
          Ucol.clear();
        }
      }


      TIMER_START(Update_L);

      for (Int supidx=0; supidx<stepSuper; supidx++){
        SuperNodeBufferTypeUnsym & snode = arrSuperNodes[supidx];

#if ( _DEBUGlevel_ >= 1 )
        statusOFS << std::endl << "["<<snode.Index<<"] "
          << "Finish updating the L part by filling LUpdateBufReduced"
          << " back to L" << std::endl << std::endl; 
#endif

        if( MYCOL( this->grid_ ) == PCOL( snode.Index, this->grid_ ) 
            && snode.LUpdateBuf.m() > 0 ){
          std::vector<LBlock<T> >& Lcol = this->L(LBj(snode.Index,this->grid_));
          //Need to skip the diagonal block if present
          Int startBlock = 
            (MYROW( this->grid_ ) == PROW( snode.Index, this->grid_ ))?1:0;
          for( Int ib = startBlock; ib < Lcol.size(); ib++ ){
            LBlock<T> & LB = Lcol[ib];
            if(1|| (LB.numRow>0 && LB.numCol>0)){
              lapack::Lacpy( 'A', LB.numRow, LB.numCol, 
                  &snode.LUpdateBuf(snode.RowLocalPtr[ib-startBlock], 0),
                  snode.LUpdateBuf.m(), LB.nzval.Data(), LB.numRow );
            }
          }
        } // Finish updating L	
      } // for (snode.Index) : Main loop


      TIMER_STOP(Update_L);




      TIMER_START(Update_U);

      for (Int supidx=0; supidx<stepSuper; supidx++){
        SuperNodeBufferTypeUnsym & snode = arrSuperNodes[supidx];

#if ( _DEBUGlevel_ >= 1 )
        statusOFS << "["<<snode.Index<<"] "
          << "Finish updating the U part by filling UUpdateBufReduced"
          << " back to U" << std::endl << std::endl; 
#endif

        if( MYROW( this->grid_ ) == PROW( snode.Index, this->grid_ ) 
            && snode.UUpdateBuf.m() > 0 ){
          std::vector<UBlock<T> >& Urow = this->U(LBi(snode.Index,this->grid_));

          for( Int jb = 0; jb < Urow.size(); jb++ ){
            UBlock<T> & UB = Urow[jb];
            if(1|| (UB.numRow>0 && UB.numCol>0)){
#if ( _DEBUGlevel_ >= 1 )
              statusOFS << "["<<snode.Index<<"] "<<"Putting colptr "
                << snode.ColLocalPtr[jb] << " of UUpdateBuf "
                << snode.ColLocalPtr << " "
                << "in UB "<<UB.blockIdx<<std::endl;
#endif
              //Indices follow L order... look at Lrow
              lapack::Lacpy( 'A', UB.numRow, UB.numCol, 
                  &snode.UUpdateBuf( 0, snode.ColLocalPtr[jb] ),
                  snode.UUpdateBuf.m(), UB.nzval.Data(), UB.numRow );
#if ( _DEBUGlevel_ >= 2 )
              statusOFS<< "["<<snode.Index<<"] "<<"UB: "<<UB.nzval<<std::endl;
#endif
            }
          }
        } // Finish updating U
      } // for (snode.Index) : Main loop


      TIMER_STOP(Update_U);




      if(next_lidx < superList.size()){
        WaitContentLCD(superList[next_lidx], superList[next_lidx].size(),nextCDBuffers);
        WaitContentUCD(superList[next_lidx], superList[next_lidx].size(),nextCDBuffers);
      }



      TIMER_START(Barrier);
      if(next_lidx < superList.size()){
        nextCDBuffers.WaitAllSend();
      }
      mpi::Waitall(arrMpireqsRecvLContentFromAny);
      mpi::Waitall(arrMpireqsRecvUContentFromAny);
      //Sync for reduce L
      //      mpi::Waitall( arrMpireqsSendToLeft );
      //Sync for reduce D
      mpi::Waitall(arrMpireqsSendToAbove);

      for (Int supidx=0; supidx<stepSuper; supidx++){
        Int ksup = superList[lidx][supidx];
        std::vector<MPI_Request> & mpireqsSendLToRight = arrMpireqsSendLToRight[supidx];
        std::vector<MPI_Request> & mpireqsSendLToBelow = arrMpireqsSendLToBelow[supidx];
        std::vector<MPI_Request> & mpireqsSendUToRight = arrMpireqsSendUToRight[supidx];
        std::vector<MPI_Request> & mpireqsSendUToBelow = arrMpireqsSendUToBelow[supidx];

        mpi::Waitall( mpireqsSendLToRight );
        mpi::Waitall( mpireqsSendLToBelow );
        mpi::Waitall( mpireqsSendUToRight );
        mpi::Waitall( mpireqsSendUToBelow );

      }
      TIMER_STOP(Barrier);



#ifdef LIST_BARRIER
#ifndef ALL_BARRIER
      if (this->options_->maxPipelineDepth!=-1)
#endif
      {
        MPI_Barrier(this->grid_->comm);
      }
#endif

    }

  template<typename T> 
    void PMatrixUnsym<T>::SelInv	(  )
    {
      this->SelInv_P2p	(  );
    } 		// -----  end of method PMatrixUnsym::SelInv  ----- 

  template<typename T> 
    void PMatrixUnsym<T>::SelInv_P2p	(  )
    {
      TIMER_START(SelInv_P2p);



      Int numSuper = this->NumSuper(); 

      // Main loop
      std::vector<std::vector<Int> > & superList = this->WorkingSet();
      Int numSteps = superList.size();

      for (Int lidx=0; lidx<numSteps ; lidx++){
        Int stepSuper = superList[lidx].size(); 
        this->SelInvIntra_P2p(lidx);
      }


      TIMER_STOP(SelInv_P2p);

      return ;
    } 		// -----  end of method PMatrixUnsym::SelInv_P2p  ----- 













  template<typename T> 
    void PMatrixUnsym<T>::PreSelInv	(  )
    {

      Int numSuper = this->NumSuper(); 

#if ( _DEBUGlevel_ >= 1 )
      statusOFS << std::endl << "L(i,k) <- L(i,k) * L(k,k)^{-1}"
        << std::endl << std::endl; 
#endif

      for( Int ksup = 0; ksup < numSuper; ksup++ ){
        if( MYCOL( this->grid_ ) == PCOL( ksup, this->grid_ ) ){
          // Broadcast the diagonal L block
          NumMat<T> nzvalLDiag;
          std::vector<LBlock<T> >& Lcol = this->L( LBj( ksup, this->grid_ ) );
          if( MYROW( this->grid_ ) == PROW( ksup, this->grid_ ) ){
            nzvalLDiag = Lcol[0].nzval;
            if( nzvalLDiag.m() != SuperSize(ksup, this->super_) ||
                nzvalLDiag.n() != SuperSize(ksup, this->super_) ){
              ErrorHandling( 
                  "The size of the diagonal block of L is wrong." );
            }
          } // Owns the diagonal block
          else {
            nzvalLDiag.Resize(SuperSize(ksup, this->super_), SuperSize(ksup, this->super_));
          }
          MPI_Bcast( (void*)nzvalLDiag.Data(), nzvalLDiag.ByteSize(),
              MPI_BYTE, PROW( ksup, this->grid_ ), this->grid_->colComm );

          // Triangular solve
          for( Int ib = 0; ib < Lcol.size(); ib++ ){
            LBlock<T> & LB = Lcol[ib];
            if( LB.blockIdx > ksup  && (1 || (LB.numCol>0 && LB.numRow>0))){
#if ( _DEBUGlevel_ >= 2 )
              // Check the correctness of the triangular solve 
              //for the first local column
              //              if( LBj( ksup, this->grid_ ) == 0 ){
              //                statusOFS << "Diag   L(" << ksup << ", " << ksup << "): "
              //                  << nzvalLDiag << std::endl;
              //                statusOFS << "Before solve L(" << LB.blockIdx << ", " << ksup 
              //                  << "): " << LB.nzval << std::endl;
              //              }
#endif
#if defined( CHECK_NORM )
              NumMat<T> Ljk = LB.nzval;
#endif
              blas::Trsm( 'R', 'L', 'N', 'U', LB.numRow, LB.numCol,
                  ONE<T>(), nzvalLDiag.Data(), LB.numCol, 
                  LB.nzval.Data(), LB.numRow );
#if defined( CHECK_NORM )
              NumMat<T> res = LB.nzval;
              //Compute L'(jk)U(kk) which should be Ljk
              blas::Trmm('R','L','N','U',res.m(),res.n(),ONE<T>(),nzvalLDiag.Data(),nzvalLDiag.m(),res.Data(),res.m());
              blas::Axpy(res.Size(), MINUS_ONE<T>(), Ljk.Data(), 1, res.Data(), 1 );
              double norm = lapack::Lange('F',res.m(),res.n(),res.Data(),res.m(),Ljk.Data());

              statusOFS << "After solve norm of residual of L(" << LB.blockIdx << ", " << ksup
                << "): " << norm << std::endl;






              // Check the correctness of the triangular solve
              // for the first local column
              //              if( LBj( ksup, this->grid_ ) == 0 ){
              //                statusOFS << "After solve  L(" << LB.blockIdx << ", " << ksup 
              //                  << "): " << LB.nzval << std::endl;
              //              }
#endif
            }
          }
        } // if( MYCOL( this->grid_ ) == PCOL( ksup, this->grid_ ) )
      } // for (ksup)




#if ( _DEBUGlevel_ >= 1 )
      statusOFS << std::endl << "U(k,j) <- U(k,k)^{-1} * U(k,j)" 
        << std::endl << std::endl; 
#endif
      for( Int ksup = 0; ksup < numSuper; ksup++ ){
        if( MYROW( this->grid_ ) == PROW( ksup, this->grid_ ) ){
          // Broadcast the diagonal L block
          NumMat<T> nzvalUDiag;
          std::vector<UBlock<T> >& Urow = this->U( LBi( ksup, this->grid_ ) );
          if( MYCOL( this->grid_ ) == PCOL( ksup, this->grid_ ) ){
            std::vector<LBlock<T> >& Lcol = this->L( LBj( ksup, this->grid_ ) );
            nzvalUDiag = Lcol[0].nzval;
            if( nzvalUDiag.m() != SuperSize(ksup, this->super_) ||
                nzvalUDiag.n() != SuperSize(ksup, this->super_) ){
              ErrorHandling( 
                  "The size of the diagonal block of U is wrong." );
            }
          } // Owns the diagonal block
          else {
            nzvalUDiag.Resize(SuperSize(ksup, this->super_), SuperSize(ksup, this->super_));
          }
          MPI_Bcast( (void*)nzvalUDiag.Data(), nzvalUDiag.ByteSize(),
              MPI_BYTE, PCOL( ksup, this->grid_ ), this->grid_->rowComm );

          // Triangular solve
          for( Int jb = 0; jb < Urow.size(); jb++ ){
            UBlock<T> & UB = Urow[jb];
            if( UB.blockIdx > ksup && (1||(UB.numCol>0 && UB.numRow>0))){
#if ( _DEBUGlevel_ >= 2 )
              // Check the correctness of the triangular solve for the first local column
              //              if( LBi( ksup, this->grid_ ) == 0 ){
              //                statusOFS << "Diag U(" << ksup << ", " << ksup << "): " 
              //                  << nzvalUDiag << std::endl;
              //                statusOFS << "Before solve U(" << ksup << ", " << UB.blockIdx
              //                  << "): " << UB.nzval << std::endl;
              //              }
#endif
#if defined( CHECK_NORM )
              NumMat<T> Ukj = UB.nzval;
#endif
              blas::Trsm( 'L', 'U', 'N', 'N', UB.numRow, UB.numCol, 
                  ONE<T>(), nzvalUDiag.Data(), UB.numRow,
                  UB.nzval.Data(), UB.numRow );
#if defined( CHECK_NORM )
              NumMat<T> res = UB.nzval;
              //Compute U(kk) * U^(kj) which should be Ukj
              blas::Trmm('L','U','N','N',res.m(),res.n(),ONE<T>(),nzvalUDiag.Data(),nzvalUDiag.m(),res.Data(),res.m());
              blas::Axpy(res.Size(), MINUS_ONE<T>(), Ukj.Data(), 1, res.Data(), 1 );
              double norm = lapack::Lange('F',res.m(),res.n(),res.Data(),res.m(),Ukj.Data());
              statusOFS << "After solve, norm of residual of U(" << ksup << ", " << UB.blockIdx
                << "): " << norm << std::endl;


              // Check the correctness of the triangular solve for the first local column
              //              if( LBi( ksup, this->grid_ ) == 0 ){
              //                statusOFS << "After solve  U(" << ksup << ", " << UB.blockIdx
              //                  << "): " << UB.nzval << std::endl;
              //              }
#endif
            }
          }
        } // if( MYROW( this->grid_ ) == PROW( ksup, this->grid_ ) )
      } // for (ksup)






#if ( _DEBUGlevel_ >= 1 )
      statusOFS << std::endl << "L(i,i) <- [L(k,k) * U(k,k)]^{-1}" << std::endl 
        << std::endl; 
#endif

      for( Int ksup = 0; ksup < numSuper; ksup++ ){
        if( MYROW( this->grid_ ) == PROW( ksup, this->grid_ ) &&
            MYCOL( this->grid_ ) == PCOL( ksup, this->grid_ )	){
          IntNumVec ipiv( SuperSize( ksup, this->super_ ) );
          // Note that the pivoting vector ipiv should follow the FORTRAN
          // notation by adding the +1
          for(Int i = 0; i < SuperSize( ksup, this->super_ ); i++){
            ipiv[i] = i + 1;
          }
          LBlock<T> & LB = (this->L( LBj( ksup, this->grid_ ) ))[0];
#if ( _DEBUGlevel_ >= 2 )
          // Check the correctness of the matrix inversion 
          // for the first local column
          //          statusOFS << "Factorized A (" << ksup << ", " << ksup << "): "
          //            << LB.nzval << std::endl;
#endif

#if defined( CHECK_NORM )
          NumMat<T> Lkk = LB.nzval;
#endif
          lapack::Getri( SuperSize( ksup, this->super_ ), LB.nzval.Data(), 
              SuperSize( ksup, this->super_ ), ipiv.Data() );

#if defined( CHECK_NORM )
          NumMat<T> res (LB.nzval.m(), LB.nzval.n());
          NumMat<T> Akk (LB.nzval.m(), LB.nzval.n());
          //rebuild A
          SetValue(Akk,ZERO<T>());
          //copy u into A
          for(Int i = 0; i<Akk.m();++i){
            for(Int j = i; j<Akk.n();++j){
              Akk(i,j)= Lkk(i,j);
            }
          }
          blas::Trmm('L','L','N','U',Akk.m(),Akk.n(),ONE<T>(),Lkk.Data(),Lkk.m(),Akk.Data(),Akk.m());

          //              statusOFS << "After inversion, original A(" << ksup << ", " << ksup
          //                  << "): " << Akk << std::endl;

          //Compute U(kk) * U'(kj) which should be Ukj
          blas::Gemm('N','N',res.m(),res.n(),res.m(),ONE<T>(),LB.nzval.Data(),res.m(),Akk.Data(),res.m(),ZERO<T>(),res.Data(),res.m());
          for(Int i = 0; i<res.m();++i){
            res(i,i)-=ONE<T>();
          }

          double norm = lapack::Lange('F',res.m(),res.n(),res.Data(),res.m(),Lkk.Data());
          statusOFS << "After inversion, norm of residual of A(" << ksup << ", " << ksup
            << "): " << norm << std::endl;




          // Check the correctness of the matrix inversion 
          // for the first local column
          //          statusOFS << "Inverted   A (" << ksup << ", " << ksup << "): " 
          //            << LB.nzval << std::endl;
#endif
        } // if I need to invert the diagonal block
      } // for (ksup)




      return ;
    } 		// -----  end of method PMatrixUnsym::PreSelInv  ----- 


  template<typename T>
    void PMatrixUnsym<T>::ConstructCommunicationPattern	(  )
    {
      ConstructCommunicationPattern_P2p();
    } 		// -----  end of method PMatrixUnsym::ConstructCommunicationPattern  ----- 



  template<typename T>
    void PMatrixUnsym<T>::ConstructCommunicationPattern_P2p	(  )
    {

      TIMER_START(ConstructCommunicationPattern);

      Int numSuper = this->NumSuper();

      TIMER_START(Allocate);

      this->isSendToBelow_.Resize(this->grid_->numProcRow, numSuper);
      this->isSendToRight_.Resize(this->grid_->numProcCol, numSuper);
      this->isSendToDiagonal_.Resize( numSuper );
      SetValue( this->isSendToBelow_, false );
      SetValue( this->isSendToRight_, false );
      SetValue( this->isSendToDiagonal_, false );

      this->isSendToCrossDiagonal_.Resize(this->grid_->numProcCol+1, numSuper );
      SetValue( this->isSendToCrossDiagonal_, false );
      this->isRecvFromCrossDiagonal_.Resize(this->grid_->numProcRow+1, numSuper );
      SetValue( this->isRecvFromCrossDiagonal_, false );

      this->isRecvFromAbove_.Resize( numSuper );
      this->isRecvFromLeft_.Resize( numSuper );
      this->isRecvFromBelow_.Resize( this->grid_->numProcRow, numSuper );
      SetValue( this->isRecvFromAbove_, false );
      SetValue( this->isRecvFromBelow_, false );
      SetValue( this->isRecvFromLeft_, false );

      TIMER_STOP(Allocate);


      TIMER_START(GetEtree);
      std::vector<Int> snodeEtree(this->NumSuper());
      this->GetEtree(snodeEtree);
      TIMER_STOP(GetEtree);

      // localColBlockRowIdx stores the nonzero block indices for each local block column.
      // The nonzero block indices including contribution from both L and U.
      // Dimension: numLocalBlockCol x numNonzeroBlock
      std::vector<std::vector<LBlock<T> > * > colBlockRowBuf; 
      colBlockRowBuf.resize( numSuper );
      std::vector<std::vector<UBlock<T> > * > rowBlockColBuf; 
      rowBlockColBuf.resize( numSuper );
      for( Int ksup = 0; ksup < numSuper; ksup++ ){
        if( MYCOL( this->grid_ ) == PCOL( ksup, this->grid_ ) 
            ||  MYROW( this->grid_ ) == PROW( ksup, this->grid_ ) ){
          colBlockRowBuf[ksup] = new std::vector<LBlock<T> >();
          if( MYROW( this->grid_ ) == PROW( ksup, this->grid_ ) ){
            rowBlockColBuf[ksup] = new std::vector<UBlock<T> >();
          }
        }
      }




      TIMER_START(Column_communication);
      for( Int ksup = 0; ksup < numSuper; ksup++ ){
        // All block columns perform independently
        if( MYCOL( this->grid_ ) == PCOL( ksup, this->grid_ ) ){
          // Communication
          std::vector< LBlock<T> > & Lcol = this->L( LBj(ksup, this->grid_ ) );
          std::vector< LBlock<T> > & LcolRecv = *colBlockRowBuf[ksup];
          std::vector<char> sendBuf;
          std::vector<char> recvBuf;
          Int localSize = 0;

          std::stringstream sstms,sstmr;
          std::vector<Int> mask( LBlockMask::TOTAL_NUMBER, 1 );
          mask[LBlockMask::NZVAL] = 0; // nzval is excluded 

          // All blocks are to be allgathered within column.
          Int numLBlock, totalLBlock;
          numLBlock = Lcol.size();
          mpi::Allreduce ( &numLBlock, &totalLBlock, 1, MPI_SUM, this->grid_->colComm);

          //serialize( (Int)Lcol.size(), sstms, NO_MASK );
          for( Int ib = 0; ib < Lcol.size(); ib++ ){
            serialize( Lcol[ib], sstms, mask );
          }
          sendBuf.resize( Size( sstms ) );
          sstms.read( &sendBuf[0], sendBuf.size() );
          localSize = sendBuf.size();

          //allgather to get the localsizes and displacements
          TIMER_START(Allgatherv_Column_communication);
          if( this->grid_ -> mpisize != 1 )
          {
            mpi::Allgatherv( sendBuf, recvBuf, this->grid_->colComm );
          }
          else{
            recvBuf = sendBuf;
          }
          TIMER_STOP(Allgatherv_Column_communication);


          sstmr.write( &recvBuf[0], recvBuf.size() );

          LcolRecv.resize(totalLBlock);

          //deserialize( numLBlock, sstmr, NO_MASK );

          for( Int ib = 0; ib < totalLBlock; ib++ ){
            deserialize( LcolRecv[ib], sstmr, mask );
          }


          //Sort and make it unique
          std::sort(LcolRecv.begin(),LcolRecv.end(),LBlockComparator<T>);
          //auto last = std::unique(LcolRecv.begin(),LcolRecv.end(),LBlockComparator<T>);
          //LcolRecv.resize(last - LcolRecv.begin());
#if ( _DEBUGlevel_ >= 1  )
          statusOFS<<"["<<ksup<<"] "<<"Lcol1: ";for(auto it = Lcol.begin();it != Lcol.end();++it){statusOFS<<it->blockIdx<<" ";}statusOFS<<endl;  
          statusOFS<<"["<<ksup<<"] "<<"LcolRecv1: ";for(auto it = LcolRecv.begin();it != LcolRecv.end();++it){statusOFS<<it->blockIdx<<" ";}statusOFS<<endl;  
#endif
        } // if( MYCOL( this->grid_ ) == PCOL( ksup, this->grid_ ) )
      }
      TIMER_STOP(Column_communication);


      TIMER_START(Row_communication);
      for( Int ksup = 0; ksup < numSuper; ksup++ ){
        // All block columns perform independently
        if( MYROW( this->grid_ ) == PROW( ksup, this->grid_ ) ){
          // Communication
          std::vector< UBlock<T> > & Urow = this->U( LBi(ksup, this->grid_ ) );
          std::vector< UBlock<T> > & UrowRecv = *rowBlockColBuf[ksup];
          std::vector<char> sendBuf;
          std::vector<char> recvBuf;
          Int localSize = 0;

          std::stringstream sstms,sstmr;
          std::vector<Int> mask( UBlockMask::TOTAL_NUMBER, 1 );
          mask[UBlockMask::NZVAL] = 0; // nzval is excluded 

          Int numUBlock, totalUBlock;
          numUBlock = Urow.size();
          mpi::Allreduce ( &numUBlock, &totalUBlock, 1, MPI_SUM, this->grid_->rowComm);
          // All blocks are to be allgathered within column.
          //serialize( (Int)Urow.size(), sstms, NO_MASK );
          for( Int jb = 0; jb < Urow.size(); jb++ ){
            serialize( Urow[jb], sstms, mask );
          }
          sendBuf.resize( Size( sstms ) );
          sstms.read( &sendBuf[0], sendBuf.size() );
          localSize = sendBuf.size();

          //allgather to get the localsizes and displacements
          TIMER_START(Allgatherv_Row_communication);
          if( this->grid_ -> mpisize != 1 )
            mpi::Allgatherv( sendBuf, recvBuf, this->grid_->rowComm );
          else
            recvBuf = sendBuf;
          TIMER_STOP(Allgatherv_Row_communication);

          sstmr.write( &recvBuf[0], recvBuf.size() );

          //deserialize( numUBlock, sstmr, NO_MASK );

          UrowRecv.resize(totalUBlock);
          for( Int jb = 0; jb < totalUBlock; jb++ ){
            deserialize( UrowRecv[jb], sstmr, mask );
          }

          //Sort and make it unique
          std::sort(UrowRecv.begin(),UrowRecv.end(),UBlockComparator<T>);
          //          auto last = std::unique(UrowRecv.begin(),UrowRecv.end(),UBlockComparator<T>);
          //          UrowRecv.resize(last - UrowRecv.begin());

#if ( _DEBUGlevel_ >= 1  )
          statusOFS<<"["<<ksup<<"] "<<"Urow1: ";for(auto it = Urow.begin();it != Urow.end();++it){statusOFS<<it->blockIdx<<" ";}statusOFS<<endl;  
          statusOFS<<"["<<ksup<<"] "<<"UrowRecv1: ";for(auto it = UrowRecv.begin();it != UrowRecv.end();++it){statusOFS<<it->blockIdx<<" ";}statusOFS<<endl;  
#endif

        } // if( MYROW( this->grid_ ) == PROW( ksup, this->grid_ ) )
      }
      TIMER_STOP(Row_communication);

      //Broadcast from diagonal processor and merge
      TIMER_START(Row_communication);
      for( Int ksup = 0; ksup < numSuper; ksup++ ){
        // All block columns perform independently
        if( MYROW( this->grid_ ) == PROW( ksup, this->grid_ ) ){
          // Communication
          std::vector< LBlock<T> > * pLcolRecv;
          std::vector< LBlock<T> > & LcolSend = *colBlockRowBuf[ksup];
          std::vector< UBlock<T> > & UrowRecv = *rowBlockColBuf[ksup];
          std::vector< LBlock<T> > Union;
          std::vector<char> sendBuf;

          std::vector<Int> mask( LBlockMask::TOTAL_NUMBER, 1 );
          mask[LBlockMask::NZVAL] = 0; // nzval is excluded 
          if( this->grid_ -> mpisize != 1 ){
            pLcolRecv = new std::vector<LBlock<T> >();
            //diagonal processor
            if( MYCOL( this->grid_ ) == PCOL( ksup, this->grid_ ) ){
              Int localSize = 0;
              std::stringstream sstms;

              serialize( (Int)LcolSend.size(), sstms, NO_MASK );
              for( Int ib = 0; ib < LcolSend.size(); ib++ ){
                serialize( LcolSend[ib], sstms, mask );
              }
              sendBuf.resize( Size( sstms ) );
              sstms.read( &sendBuf[0], sendBuf.size() );
              localSize = sendBuf.size();
            }
            TIMER_START(Bcast_Row_communication);
            mpi::Bcast( sendBuf, PCOL(ksup,this->grid_), this->grid_->rowComm );
            TIMER_STOP(Bcast_Row_communication);

            std::stringstream sstmr;
            sstmr.write( &sendBuf[0], sendBuf.size() );

            Int numLBlock;
            deserialize( numLBlock, sstmr, NO_MASK );
            pLcolRecv->resize(numLBlock);
            for( Int ib = 0; ib < numLBlock; ib++ ){
              deserialize( (*pLcolRecv)[ib], sstmr, mask );
            }
          }
          else{
            pLcolRecv = &LcolSend;
          }

          //LcolRecv is sorted and with unique elements
          //UrowRecv is sorted and with unique elements

          auto result = back_inserter(Union);
          auto firstUrow = UrowRecv.begin();
          auto lastUrow = UrowRecv.end();
          auto firstLcol = pLcolRecv->begin();
          auto lastLcol = pLcolRecv->end();
          while (true)
          {
            if (firstUrow==lastUrow){
              result = std::copy(firstLcol,lastLcol,result);
              break;
            }
            if (firstLcol==lastLcol){
              for(auto it = firstUrow; it != lastUrow; ++it){
                LBlock<T> LB;
                LB.blockIdx = it->blockIdx;
                LB.numRow = it->numCol;
                LB.numCol = it->numRow;
                LB.rows = it->cols;
                *result = LB;
                ++result;
              }
              break;
            }

            if (firstUrow->blockIdx<firstLcol->blockIdx) {
              LBlock<T> LB;
              LB.blockIdx = firstUrow->blockIdx;
              LB.numRow = firstUrow->numCol;
              LB.numCol = firstUrow->numRow;
              LB.rows = firstUrow->cols;
              *result = LB;
              ++firstUrow;
            }
            else if (firstLcol->blockIdx<firstUrow->blockIdx) { 
              *result = *firstLcol; 
              ++firstLcol; 
            }
            else {
              //compute set union between rows and cols
              std::set<Int> uRowCol;
              uRowCol.insert(&firstLcol->rows[0],
                  &firstLcol->rows[0]+firstLcol->numRow);
              uRowCol.insert(&firstUrow->cols[0],
                  &firstUrow->cols[0]+firstUrow->numCol);

              LBlock<T> LB;
              LB.blockIdx = firstUrow->blockIdx;
              LB.numRow = uRowCol.size();
              LB.numCol = firstUrow->numRow;
              LB.rows.Resize(uRowCol.size());
              std::copy(uRowCol.begin(),uRowCol.end(),&LB.rows[0]);
              *result = LB;
              ++firstUrow; 
              ++firstLcol; 
            }
            ++result;
          }

#if ( _DEBUGlevel_ >= 1  )
          statusOFS<<"["<<ksup<<"] "<<"Lcol: ";for(auto it = pLcolRecv->begin();it != pLcolRecv->end();++it){statusOFS<<it->blockIdx<<" ";}statusOFS<<endl;  
          statusOFS<<"["<<ksup<<"] "<<"Urow: ";for(auto it = UrowRecv.begin();it != UrowRecv.end();++it){statusOFS<<it->blockIdx<<" ";}statusOFS<<endl;  
          statusOFS<<"["<<ksup<<"] "<<"Union: ";for(auto it = Union.begin();it != Union.end();++it){statusOFS<<it->blockIdx<<" ";}statusOFS<<endl;  
#endif

          //replace LcolSend by the union
          LcolSend.swap(Union);

          if( this->grid_ -> mpisize != 1 ){
            delete pLcolRecv;
          }
        } // if( MYROW( this->grid_ ) == PROW( ksup, this->grid_ ) )
      }
      TIMER_STOP(Row_communication);


      TIMER_START(Col_communication);
      for( Int ksup = 0; ksup < numSuper; ksup++ ){
        // All block columns perform independently
        if( MYCOL( this->grid_ ) == PCOL( ksup, this->grid_ ) ){
          // Communication
          std::vector< LBlock<T> > * pUnionRecv;
          std::vector< LBlock<T> > & UnionSend = *colBlockRowBuf[ksup];

          std::vector<char> sendBuf;

          std::vector<Int> mask( LBlockMask::TOTAL_NUMBER, 1 );
          mask[LBlockMask::NZVAL] = 0; // nzval is excluded 
          if( this->grid_ -> mpisize != 1 ){
            pUnionRecv = new std::vector<LBlock<T> >();
            //diagonal processor
            if( MYROW( this->grid_ ) == PROW( ksup, this->grid_ ) ){
              Int localSize = 0;
              std::stringstream sstms;

              serialize( (Int)UnionSend.size(), sstms, NO_MASK );
              for( Int ib = 0; ib < UnionSend.size(); ib++ ){
                serialize( UnionSend[ib], sstms, mask );
              }
              sendBuf.resize( Size( sstms ) );
              sstms.read( &sendBuf[0], sendBuf.size() );
              localSize = sendBuf.size();
            }
            TIMER_START(Bcast_Col_communication);
            mpi::Bcast( sendBuf, PROW(ksup,this->grid_), this->grid_->colComm );
            TIMER_STOP(Bcast_Col_communication);

            std::stringstream sstmr;
            sstmr.write( &sendBuf[0], sendBuf.size() );

            Int numLBlock;
            deserialize( numLBlock, sstmr, NO_MASK );
            pUnionRecv->resize(numLBlock);
            for( Int ib = 0; ib < numLBlock; ib++ ){
              deserialize( (*pUnionRecv)[ib], sstmr, mask );
            }


            UnionSend.resize(pUnionRecv->size());
            std::copy(pUnionRecv->begin(),pUnionRecv->end(),UnionSend.begin());
          }
          else{
            pUnionRecv = &UnionSend;
          }

#if ( _DEBUGlevel_ >= 1  )
          statusOFS<<"["<<ksup<<"] "<<"Union: ";for(auto it = UnionSend.begin();it != UnionSend.end();++it){statusOFS<<it->blockIdx<<" ";}statusOFS<<endl;  
#endif

          if( this->grid_ -> mpisize != 1 ){
            delete pUnionRecv;
          }
        } // if( MYCOL( this->grid_ ) == PCOL( ksup, this->grid_ ) )
      }
      TIMER_STOP(Col_communication);


      for( Int ksup = 0; ksup < numSuper; ksup++ ){
        // All block columns perform independently
        if( MYROW( this->grid_ ) == PROW( ksup, this->grid_ ) ){

          std::vector< LBlock<T> > & Union = *colBlockRowBuf[ksup];
          std::vector< UBlock<T> > & Urow = this->U( LBi(ksup, this->grid_ ) );
          Int & LrowSize = this->LrowSize_[ LBi(ksup, this->grid_ ) ];
          TIMER_START(Computing_Lrow_size);
          //Allocate Lrow and extend Urow
          LrowSize = 0;
          for(auto it = Union.begin(); it!=Union.end();++it){
            Int Idx = it->blockIdx;
            if(Idx > ksup && (Idx % this->grid_->numProcCol) == MYCOL(this->grid_)  ){ 
              ++LrowSize;
            }
          }
          TIMER_STOP(Computing_Lrow_size);

          TIMER_START(Extending_Urow);
          Urow.reserve(LrowSize);

          for(auto it = Union.begin(); it!=Union.end();++it){
            Int Idx = it->blockIdx;
            if(Idx > ksup && (Idx % this->grid_->numProcCol) == MYCOL(this->grid_)  ){ 
              bool isFound = false;
              Int nextIdx = -1;
              for(Int jb = 0; jb < Urow.size(); ++jb){
                UBlock<T> & UB = Urow[jb];
                if(UB.blockIdx == Idx){
                  isFound = true;
                  nextIdx = jb;
                  break;
                }
                if(UB.blockIdx > Idx){
                  nextIdx = jb;
                  break;
                }
              }
              if(!isFound){
                //push_back
                UBlock<T> UB;
                UB.blockIdx = Idx;
                UB.numRow = it->numCol;
                UB.numCol = it->numRow;
                UB.cols = it->rows;
                Urow.push_back(UB);
                UBlock<T> & UBl = Urow.back();
                UBl.nzval.Resize(UB.numRow,UB.numCol);
                SetValue(UBl.nzval,ZERO<T>());
              }
              else{
                //make sure blocks are the same size
                UBlock<T> & UB = Urow[nextIdx];
                assert(UB.numRow == it->numCol);
                if( UB.numCol != it->numRow ){
                  NumMat<T> tmpNzval = UB.nzval;
                  IntNumVec tmpCols = UB.cols;

                  UB.numRow = it->numCol;
                  UB.numCol = it->numRow;
                  UB.cols = it->rows;
                  UB.nzval.Resize(UB.numRow,UB.numCol);
                  SetValue(UB.nzval,ZERO<T>());

                  //now put nzvals back in place
                  Int jOldCols = 0;
                  for(Int j = 0; j<UB.numCol; ++j){
                    Int newCol = UB.cols[j];
                    if(jOldCols<tmpCols.m()){
                      Int oldCol = tmpCols[jOldCols];
                      if(newCol == oldCol){
                        T * nzcolPtr = tmpNzval.VecData(jOldCols);
                        std::copy(nzcolPtr,nzcolPtr+UB.numRow,UB.nzval.VecData(j));
                        jOldCols++;
                      }
                    }
                    else{
                      break;
                    }
                  }
                  assert(jOldCols>=tmpCols.m());

                }

              }
            }
          }
#if ( _DEBUGlevel_ >= 1  )
          statusOFS<<"["<<ksup<<"] "<<"LrowSize = "<<LrowSize<<std::endl;
#endif
          std::sort(Urow.begin(),Urow.end(),UBlockComparator<T>);
#if ( _DEBUGlevel_ >= 1  )
          statusOFS<<"["<<ksup<<"] "<<"Urow extended: ";for(auto it = Urow.begin();it != Urow.end();++it){statusOFS<<it->blockIdx<<" ";}statusOFS<<endl;  
#endif
          TIMER_STOP(Extending_Urow);
        } // if( MYROW( this->grid_ ) == PROW( ksup, this->grid_ ) )
      } // for(ksup)

      for( Int ksup = 0; ksup < numSuper; ksup++ ){
        // All block columns perform independently
        if( MYCOL( this->grid_ ) == PCOL( ksup, this->grid_ ) ){

          std::vector< LBlock<T> > & Union = *colBlockRowBuf[ksup];
          std::vector< LBlock<T> > & Lcol = this->L( LBj(ksup, this->grid_ ) );
          Int & UcolSize = this->UcolSize_[ LBj(ksup, this->grid_ ) ];
          TIMER_START(Computing_Ucol_size);
          //Allocate Ucol and extend Lcol
          UcolSize = 0;
          for(auto it = Union.begin(); it!=Union.end();++it){
            Int Idx = it->blockIdx;
            if(Idx > ksup && (Idx % this->grid_->numProcRow) == MYROW(this->grid_)  ){ 
              ++UcolSize;
            }
          }
          TIMER_STOP(Computing_Ucol_size);

          TIMER_START(Extending_Lcol);
          Lcol.reserve(UcolSize);

          for(auto it = Union.begin(); it!=Union.end();++it){
            Int Idx = it->blockIdx;
            if(Idx > ksup && (Idx % this->grid_->numProcRow) == MYROW(this->grid_)  ){ 
              bool isFound = false;
              Int nextIdx = -1;
              for(Int ib = 0; ib < Lcol.size(); ++ib){
                LBlock<T> & LB = Lcol[ib];
                if(LB.blockIdx == Idx){
                  isFound = true;
                  nextIdx = ib;
                  break;
                }
                if(LB.blockIdx > Idx){
                  nextIdx = ib;
                  break;
                }
              }
              if(!isFound){
                //push_back
                Lcol.push_back(*it);
                LBlock<T> & LB = Lcol.back();
                LB.nzval.Resize(LB.numRow,LB.numCol);
                SetValue(LB.nzval,ZERO<T>());
              }
              else{
                //make sure blocks are the same size
                LBlock<T> & LB = Lcol[nextIdx];
                assert(LB.numCol == it->numCol);
                if( LB.numRow != it->numRow ){
                  NumMat<T> tmpNzval = LB.nzval;
                  IntNumVec tmpRows = LB.rows;

                  LB.numRow = it->numRow;
                  LB.numCol = it->numCol;
                  LB.rows = it->rows;
                  LB.nzval.Resize(LB.numRow,LB.numCol);
                  SetValue(LB.nzval,ZERO<T>());

                  //now put nzvals back in place
                  Int iOldRows = 0;
                  for(Int i = 0; i<LB.numRow; ++i){
                    Int newRow = LB.rows[i];
                    if(iOldRows<tmpRows.m()){
                      Int oldRow = tmpRows[iOldRows];
                      if(newRow == oldRow){
                        for(Int j = 0; j<LB.numCol; ++j){
                          LB.nzval(i,j) = tmpNzval(iOldRows,j);
                        }
                        iOldRows++;
                      }
                    }
                    else{
                      break;
                    }
                  }
                  assert(iOldRows>=tmpRows.m());

                }

              }
            }
          }
#if ( _DEBUGlevel_ >= 1  )
          statusOFS<<"["<<ksup<<"] "<<"UcolSize = "<<UcolSize<<std::endl;
#endif
          std::sort(Lcol.begin(),Lcol.end(),LBlockComparator<T>);
#if ( _DEBUGlevel_ >= 1  )
          statusOFS<<"["<<ksup<<"] "<<"Lcol extended: ";for(auto it = Lcol.begin();it != Lcol.end();++it){statusOFS<<it->blockIdx<<" ";}statusOFS<<endl;  
#endif
          TIMER_STOP(Extending_Lcol);
        } // if( MYCOL( this->grid_ ) == PCOL( ksup, this->grid_ ) )
      } // for(ksup)





      //pointers to next non zero block
      std::vector<Int> nextNZColBlockRow;
      nextNZColBlockRow.resize( this->NumLocalBlockCol(), 0 );

      std::vector<Int> nextNZRowBlockCol;
      nextNZRowBlockCol.resize( this->NumLocalBlockRow() , 0);

      if(1){
        for( Int ksup = 0; ksup < numSuper; ksup++ ){
          // All block columns perform independently
          std::vector<Int> tAllBlockRowIdx;
          if( MYCOL( this->grid_ ) == PCOL( ksup, this->grid_ ) ){
            // Communication
            std::vector< LBlock<T> > & Union = *colBlockRowBuf[ksup];

            std::vector<Int> tlocalBlockRowIdx;
            //look at previous supernodes (ideally look at children in the etree only)

            for( Int jsup = 0; jsup < ksup; jsup++ ){
              //look into U
              Int pjsup = snodeEtree[jsup];
              if(pjsup<=ksup && MYROW( this->grid_ ) == PROW( jsup, this->grid_ ) ){
                //            if( MYCOL( this->grid_ ) == PCOL( jsup, this->grid_ ) ){
                std::vector< UBlock<T> > & Urow = this->U( LBi(jsup, this->grid_ ) );



#if ( _DEBUGlevel_ >= 1  )
                statusOFS<<"["<<ksup<<"] "<<"Urow["<<jsup<<"]: "; for(auto it = Urow.begin(); it != Urow.end(); ++it){ statusOFS<<it->blockIdx<<" "; } statusOFS<<endl;
#endif

                Int & nextIdx = nextNZRowBlockCol[ LBi(jsup, this->grid_) ]; 
                while(nextIdx<Urow.size()){ 
                  if(Urow[nextIdx].blockIdx < ksup){ 
                    ++nextIdx; 
                  }
                  else{
                    break;
                  }
                }
                if(nextIdx<Urow.size()){
                  if(Urow[nextIdx].blockIdx == ksup){
                    tlocalBlockRowIdx.push_back(jsup);
                  }
                }
                //            }
              }
            }





            TIMER_START(Allgatherv_Column_communication);
            if( this->grid_ -> mpisize != 1 )
              mpi::Allgatherv( tlocalBlockRowIdx, tAllBlockRowIdx, this->grid_->colComm );
            else
              tAllBlockRowIdx = tlocalBlockRowIdx;

            TIMER_STOP(Allgatherv_Column_communication);

            //add dummy LBlocks to Union


            for(auto it = tAllBlockRowIdx.begin(); it != tAllBlockRowIdx.end(); ++it ){
              LBlock<T> LB;
              LB.numCol = 0;
              LB.numRow = 0;
              LB.blockIdx = *it;
              Union.push_back(LB);
            }

            std::sort(Union.begin(),Union.end(),LBlockComparator<T>);

#if ( _DEBUGlevel_ >= 1  )
            statusOFS<<"["<<ksup<<"] "<<"NEW Union: ";
            for(auto it = Union.begin(); it != Union.end(); ++it){ statusOFS<<it->blockIdx<<" "; }
            statusOFS<<endl;
#endif

          } // if( MYCOL( this->grid_ ) == PCOL( ksup, this->grid_ ) )
        }

      }

      for( Int ksup = 0; ksup < numSuper; ksup++ ){
        // All block columns perform independently
        std::vector<Int> tAllBlockColIdx;
        if( MYROW( this->grid_ ) == PROW( ksup, this->grid_ ) ){
          // Communication
          std::vector< LBlock<T> > & Union = *colBlockRowBuf[ksup];
          std::vector<Int> tlocalBlockColIdx;
          //look at previous supernodes (ideally look at children in the etree only)
          for( Int jsup = 0; jsup < ksup; jsup++ ){
            //look into U
            Int pjsup = snodeEtree[jsup];
            if(pjsup<=ksup && MYCOL( this->grid_ ) == PCOL( jsup, this->grid_ ) ){
              std::vector< LBlock<T> > & Lcol = this->L( LBj(jsup, this->grid_ ) );



#if ( _DEBUGlevel_ >= 1  )
              statusOFS<<"["<<ksup<<"] "<<"Lcol["<<jsup<<"]: "; for(auto it = Lcol.begin(); it != Lcol.end(); ++it){ statusOFS<<it->blockIdx<<" "; } statusOFS<<endl;
#endif

              Int & nextIdx = nextNZColBlockRow[ LBj(jsup, this->grid_) ]; 
              while(nextIdx<Lcol.size()){ 
                if(Lcol[nextIdx].blockIdx < ksup){ 
                  ++nextIdx; 
                }
                else{
                  break;
                }
              }
              if(nextIdx<Lcol.size()){
                if(Lcol[nextIdx].blockIdx == ksup){
                  tlocalBlockColIdx.push_back(jsup);
                }
              }
            }
          }










          TIMER_START(Allgatherv_Row_communication);
          if( this->grid_ -> mpisize != 1 )
            mpi::Allgatherv( tlocalBlockColIdx, tAllBlockColIdx, this->grid_->rowComm );
          else
            tAllBlockColIdx = tlocalBlockColIdx;

          TIMER_STOP(Allgatherv_Row_communication);

          //add dummy LBlocks to Union


          for(auto it = tAllBlockColIdx.begin(); it != tAllBlockColIdx.end(); ++it ){
            LBlock<T> LB;
            LB.numCol = 0;
            LB.numRow = 0;
            LB.blockIdx = *it;
            Union.push_back(LB);
          }

          std::sort(Union.begin(),Union.end(),LBlockComparator<T>);
          auto last = std::unique(Union.begin(),Union.end(),LBlockEqualComparator<T>);
          Union.erase(last,Union.end());

#if ( _DEBUGlevel_ >= 1  )
          statusOFS<<"["<<ksup<<"] "<<"NEW Union: ";
          for(auto it = Union.begin(); it != Union.end(); ++it){ statusOFS<<it->blockIdx<<" "; }
          statusOFS<<endl;
#endif

        } // if( MYROW( this->grid_ ) == PROW( ksup, this->grid_ ) )
      }





      TIMER_START(STB_RFA);
      for( Int ksup = 0; ksup < numSuper - 1; ksup++ ){
        // Loop over all the supernodes to the right of ksup


        Int jsup = snodeEtree[ksup];
        while(jsup<numSuper){
          Int jsupLocalBlockCol = LBj( jsup, this->grid_ );
          Int jsupProcCol = PCOL( jsup, this->grid_ );
          if( MYCOL( this->grid_ ) == jsupProcCol ){
            std::vector< LBlock<T> > & Union = *colBlockRowBuf[jsup];

            // SendToBelow / RecvFromAbove only if (ksup, jsup) is nonzero.
            Int isFound = 0;
            for(auto it = Union.begin();it!=Union.end();++it){
              if(it->blockIdx == ksup){
                isFound = 1;
                break;
              }
              else if(it->blockIdx > ksup){
                break;
              }
            }
            if( isFound ) {
              for( auto si = Union.begin(); si != Union.end(); si++	 ){
                Int isup = si->blockIdx;
                Int isupProcRow = PROW( isup, this->grid_ );
                if( isup > ksup ){
                  if( MYROW( this->grid_ ) == isupProcRow ){
                    this->isRecvFromAbove_(ksup) = true;
                  }
                  if( MYROW( this->grid_ ) == PROW( ksup, this->grid_ ) ){
                    this->isSendToBelow_( isupProcRow, ksup ) = true;
                  }
                } // if( isup > ksup )
              } // for (si)
            } // if( isFound )

          } // if( MYCOL( this->grid_ ) == PCOL( jsup, this->grid_ ) )
          jsup = snodeEtree[jsup];
        } // for(jsup)
      } // for(ksup)

#if ( _DEBUGlevel_ >= 1 )
      statusOFS << std::endl << "this->isSendToBelow:" << std::endl;
      for(int j = 0;j< this->isSendToBelow_.n();j++){
        statusOFS << "["<<j<<"] ";
        for(int i =0; i < this->isSendToBelow_.m();i++){
          statusOFS<< this->isSendToBelow_(i,j) << " ";
        }
        statusOFS<<std::endl;
      }

      statusOFS << std::endl << "this->isRecvFromAbove:" << std::endl;
      for(int j = 0;j< this->isRecvFromAbove_.m();j++){
        statusOFS << "["<<j<<"] "<< this->isRecvFromAbove_(j)<<std::endl;
      }
#endif



      TIMER_STOP(STB_RFA);








      TIMER_START(STR_RFL_RFB);


      for( Int ksup = 0; ksup < numSuper - 1; ksup++ ){
        // Loop over all the supernodes below ksup

        Int isup = snodeEtree[ksup];
        while(isup<numSuper){
          Int isupLocalBlockRow = LBi( isup, this->grid_ );
          Int isupProcRow       = PROW( isup, this->grid_ );
          if( MYROW( this->grid_ ) == isupProcRow ){
            std::vector< LBlock<T> > & Union = *colBlockRowBuf[isup];

            // SendToRight / RecvFromLeft only if (isup, ksup) is nonzero.
            Int isFound = 0;
            for(auto it = Union.begin();it!=Union.end();++it){
              if(it->blockIdx == ksup){
                isFound = 1;
                break;
              }
              else if(it->blockIdx > ksup){
                break;
              }
            }
            if( isFound ) {
              for( auto si = Union.begin(); si != Union.end(); si++ ){
                Int jsup = si->blockIdx;
                Int jsupProcCol = PCOL( jsup, this->grid_ );
                if( jsup > ksup ){
                  if( MYCOL( this->grid_ ) == jsupProcCol ){
                    this->isRecvFromLeft_(ksup) = true;
                  }
                  if( MYCOL( this->grid_ ) == PCOL( ksup, this->grid_ ) ){
                    this->isSendToRight_( jsupProcCol, ksup ) = true;
                  }
                }
              } // for (si)
            } // if( isFound )
          } // if( MYROW( this->grid_ ) == isupProcRow )

          if( MYCOL( this->grid_ ) == PCOL(ksup, this->grid_) ){
            if( MYROW( this->grid_ ) == PROW( ksup, this->grid_ ) ){ 
              this->isRecvFromBelow_(isupProcRow,ksup) = true;
            }    
            else if (MYROW(this->grid_) == isupProcRow){
              this->isSendToDiagonal_(ksup)=true;
            }    
          } // if( MYCOL( this->grid_ ) == PCOL(ksup, this->grid_) )
          isup = snodeEtree[isup];
        } // for (isup)
      }	 // for (ksup)

#if ( _DEBUGlevel_ >= 1 )
      statusOFS << std::endl << "this->isSendToRight:" << std::endl;
      for(int j = 0;j< this->isSendToRight_.n();j++){
        statusOFS << "["<<j<<"] ";
        for(int i =0; i < this->isSendToRight_.m();i++){
          statusOFS<< this->isSendToRight_(i,j) << " ";
        }
        statusOFS<<std::endl;
      }

      statusOFS << std::endl << "this->isRecvFromLeft:" << std::endl;
      for(int j = 0;j< this->isRecvFromLeft_.m();j++){
        statusOFS << "["<<j<<"] "<< this->isRecvFromLeft_(j)<<std::endl;
      }

      statusOFS << std::endl << "this->isRecvFromBelow:" << std::endl;
      for(int j = 0;j< this->isRecvFromBelow_.n();j++){
        statusOFS << "["<<j<<"] ";
        for(int i =0; i < this->isRecvFromBelow_.m();i++){
          statusOFS<< this->isRecvFromBelow_(i,j) << " ";
        }
        statusOFS<<std::endl;
      }
#endif


      TIMER_STOP(STR_RFL_RFB);




      TIMER_START(STCD_RFCD);


      for( Int ksup = 0; ksup < numSuper - 1; ksup++ ){
        if( MYCOL( this->grid_ ) == PCOL( ksup, this->grid_ ) ){
          std::vector< LBlock<T> > & Union = *colBlockRowBuf[ksup];
          for( auto si = Union.begin(); si != Union.end(); si++ ){
            Int isup = si->blockIdx;
            Int isupProcRow = PROW( isup, this->grid_ );
            Int isupProcCol = PCOL( isup, this->grid_ );
            if( isup > ksup && MYROW( this->grid_ ) == isupProcRow ){
              this->isSendToCrossDiagonal_(this->grid_->numProcCol, ksup ) = true;
              this->isSendToCrossDiagonal_(isupProcCol, ksup ) = true;
            }
          } // for (si)
        } // if( MYCOL( this->grid_ ) == PCOL( ksup, this->grid_ ) )
      } // for (ksup)

      for( Int ksup = 0; ksup < numSuper - 1; ksup++ ){
        if( MYROW( this->grid_ ) == PROW( ksup, this->grid_ ) ){
          std::vector< LBlock<T> > & Union = *colBlockRowBuf[ksup];
          for( auto si = Union.begin(); si != Union.end(); si++ ){
            Int jsup = si->blockIdx;
            Int jsupProcCol = PCOL( jsup, this->grid_ );
            Int jsupProcRow = PROW( jsup, this->grid_ );
            if( jsup > ksup && MYCOL(this->grid_) == jsupProcCol ){
              this->isRecvFromCrossDiagonal_(this->grid_->numProcRow, ksup ) = true;
              this->isRecvFromCrossDiagonal_(jsupProcRow, ksup ) = true;
            }
          } // for (si)
        } // if( MYROW( this->grid_ ) == PROW( ksup, this->grid_ ) )
      } // for (ksup)
#if ( _DEBUGlevel_ >= 1 )
      statusOFS << std::endl << "this->isSendToCrossDiagonal:" << std::endl;
      for(int j =0; j < this->isSendToCrossDiagonal_.n();j++){
        if(this->isSendToCrossDiagonal_(this->grid_->numProcCol,j)){
          statusOFS << "["<<j<<"] ";
          for(int i =0; i < this->isSendToCrossDiagonal_.m()-1;i++){
            if(this->isSendToCrossDiagonal_(i,j))
            {
              statusOFS<< PNUM(PROW(j,this->grid_),i,this->grid_)<<" ";
            }
          }
          statusOFS<<std::endl;
        }
      }

      statusOFS << std::endl << "this->isRecvFromCrossDiagonal:" << std::endl;
      for(int j =0; j < this->isRecvFromCrossDiagonal_.n();j++){
        if(this->isRecvFromCrossDiagonal_(this->grid_->numProcRow,j)){
          statusOFS << "["<<j<<"] ";
          for(int i =0; i < this->isRecvFromCrossDiagonal_.m()-1;i++){
            if(this->isRecvFromCrossDiagonal_(i,j))
            {
              statusOFS<< PNUM(i,PCOL(j,this->grid_),this->grid_)<<" ";
            }
          }
          statusOFS<<std::endl;
        }
      }


#endif


      TIMER_STOP(STCD_RFCD);

      for( Int ksup = 0; ksup < numSuper; ksup++ ){
        if( MYCOL( this->grid_ ) == PCOL( ksup, this->grid_ ) 
            ||  MYROW( this->grid_ ) == PROW( ksup, this->grid_ ) ){
          delete colBlockRowBuf[ksup];
          if( MYROW( this->grid_ ) == PROW( ksup, this->grid_ ) ){
            delete rowBlockColBuf[ksup];
          }
        }
      }

      //Build the list of supernodes based on the elimination tree from SuperLU
      this->GetWorkSet(snodeEtree,this->WorkingSet());

      TIMER_STOP(ConstructCommunicationPattern);


      return ;
    } 		// -----  end of method PMatrixUnsym::ConstructCommunicationPattern_P2p  ----- 


  } // namespace PEXSI

#endif //_PEXSI_PSELINV_UNSYM_IMPL_HPP_
