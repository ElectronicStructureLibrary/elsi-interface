/*
   Copyright (c) 2012 The Regents of the University of California,
   through Lawrence Berkeley National Laboratory.  

   Author: Mathias Jacquelin and Lin Lin

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
/// @file ngchol_interf_impl.hpp
/// @brief Implementation of interface with NGCHOL.
/// @date 2014-07-08 Original version
#ifndef _PEXSI_NGCHOL_INTERF_IMPL_HPP_
#define _PEXSI_NGCHOL_INTERF_IMPL_HPP_

// Interface with NGCHOL
#include "ngchol.hpp"
#include "ngchol/sp_blas.hpp"

// Interface with PSelInv
#include "pexsi/pselinv.hpp"

// Interface with sparse matrix (CSC format)
#include "pexsi/sparse_matrix.hpp"
#include "pexsi/environment.hpp"
#include "pexsi/sparse_matrix.hpp"
#include "pexsi/NumMat.hpp"
#include "pexsi/NumVec.hpp"

// Interface with LAPACK
#include "pexsi/lapack.hpp"


#include <list>

namespace PEXSI{

template<typename T> void NGCHOLMatrixToSuperNode( 
    LIBCHOLESKY::SupernodalMatrix<T>& SMat,
    SuperNodeType& super ){
#ifndef _RELEASE_
	PushCallStack("NGCHOLMatrixToSuperNode");
#endif
  Int n = SMat.Size();

  // perm
  const LIBCHOLESKY::IntNumVec& SPerm = SMat.GetOrdering().perm;
  super.perm.Resize( SPerm.m() );
  for( Int i = 0; i < SPerm.m(); i++ ){
    super.perm[i] = SPerm[i];
  }
  
  // permInv
  super.permInv.Resize( n );
  for( Int i = 0; i < n; i++ ){
    super.permInv[i] = i;
  }
  std::sort( super.permInv.Data(), super.permInv.Data() + n,
      IndexComp<IntNumVec&>(super.perm) );

  LIBCHOLESKY::IntNumVec& XSuper = SMat.GetSupernodalPartition();
  Int numSuper = XSuper.m() - 1;

  // superPtr
  IntNumVec& superPtr = super.superPtr;
  superPtr.Resize( numSuper + 1 );
  for( Int i = 0; i < numSuper + 1; i++ ){
    superPtr[i] = XSuper[i] - 1;
  }

  // superIdx
  IntNumVec& superIdx = super.superIdx;
  superIdx.Resize( n );
  const LIBCHOLESKY::IntNumVec& superMember = SMat.GetSupMembership();
  for( Int i = 0; i < n; i++ ){
    superIdx(i) = superMember[i] - 1;
  }

  // etree
  LIBCHOLESKY::ETree& etree = SMat.GetETree();
  super.etree.Resize(n);
  for( Int i = 0; i < n; i++ ){
    super.etree[i] = etree.PostParent(i);
  }

#ifndef _RELEASE_
	PopCallStack();
#endif
}  // -----  end of NGCHOLMatrixToSuperNode ----- 



template<typename T> void NGCHOLMatrixToPMatrix( 
    LIBCHOLESKY::SupernodalMatrix<T>& SMat,
    PMatrix<T>& PMat ){
#ifndef _RELEASE_
	PushCallStack("NGCHOLMatrixToPMatrix");
#endif
  // This routine assumes that the g, supernode and options of PMatrix
  // has been set outside this routine.
  
  Int mpirank, mpisize;
  const GridType *g = PMat.Grid();

  // FIXME Check PMatrix and SupernodalMatrix has the same communicator
  MPI_Comm comm = g->comm;
  MPI_Comm colComm = g->colComm;
  MPI_Comm rowComm = g->rowComm;
  MPI_Comm_size(comm, &mpisize); 
  MPI_Comm_rank(comm, &mpirank);

  Int nprow = g->numProcRow, npcol = g->numProcCol;

  Int mpirow = mpirank / npcol;  
  Int mpicol = mpirank % npcol;

  Int n = SMat.Size();
  PMat.ColBlockIdx().clear();
  PMat.RowBlockIdx().clear();
  PMat.ColBlockIdx().resize( PMat.NumLocalBlockCol() );
  PMat.RowBlockIdx().resize( PMat.NumLocalBlockRow() );

  const IntNumVec& superIdx = PMat.SuperNode()->superIdx;
#if ( _DEBUGlevel_ >= 1 )
      statusOFS << "superIdx = " << superIdx << std::endl;
#endif


  // for loop over all supernodes
  //
  //   if the current processor owns the supernode (the ownership is block
  //   cyclic)
  //     serialize the information
  //     broadcast the information to all processors
  //   elseif 
  //     receive the information from the processor owning the
  //     supernode, and save the information in a deserialized buffer
  //   endif
  //
  //   (Local operation from now)
  //   if the current processor owns the correct processor column
  //     for loop over the local blocks
  //       for loop over each row
  //         if the current row number belong to the current processor
  //           add the row index to a vector for LBuf
  //         endif
  //       endfor
  //     endfor
  //
  //     Allocate the sign of LBuf for saving the nonzero values
  //     Convert the format of rowind
  //
  //     for loop over the local blocks
  //       for loop over each row
  //         if the current row number belong to the current processor
  //           Append the nonzero values to nzval
  //         endif
  //       endfor
  //     endfor
  //   endif
  //
  //   discard the temporary information in the buffer for all
  //   processors
  //
  // endfor
  //
  // Perform SendToCrossDiagonal to fill the U part.
  //
  // Remark:
  //   1.  The communiation is performed in a supernode-by-supernode
  //   way.  This may not be very fast.  A first improvement is to
  //   broadcast all local supernodes at once and then do local
  //   post-processing.

  Int numSuper = PMat.NumSuper();
  LIBCHOLESKY::Icomm snodeIcomm;
  std::vector<char> buffer;
  for( Int iSuper = 0; iSuper < numSuper; iSuper++ ){
    LIBCHOLESKY::SuperNode<T> snode;

    if( mpirank == ( iSuper % mpisize ) ){
      // Get the local supernode
      Int iSuperLocal = iSuper / mpisize;
      LIBCHOLESKY::SuperNode<T>& snodetmp = 
        SMat.GetLocalSupernode(iSuperLocal);
#if ( _DEBUGlevel_ >= 1 )
      statusOFS << "iSuper = " << iSuper << ", iSuperLocal = " <<
        iSuperLocal << ", id = " << snodetmp.Id() << ", size = " << 
        snodetmp.Size() << ", #Block = " << snodetmp.NZBlockCnt() <<
        std::endl;
    statusOFS << "snode (before bcast) = " << snodetmp << std::endl;
#endif
      // Serialize the information in the current supernode
//      std::stringstream sstm;
//      serialize( snodetmp.Id(), sstm, NO_MASK );
//      serialize( snodetmp.Size(), sstm, NO_MASK );
//      serialize( snodetmp.NZBlockCnt(), sstm, NO_MASK );
//      for( Int blkidx = 0; blkidx < snodetmp.NZBlockCnt(); blkidx++){
//        LIBCHOLESKY::NZBlockDesc & nzblk_desc =
//          snodetmp.GetNZBlockDesc( blkidx );
//      } // for (blkidx)

      LIBCHOLESKY::Serialize( snodeIcomm, snodetmp );
      Int msgSize = snodeIcomm.size();
#if ( _DEBUGlevel_ >= 1 )
      statusOFS << "msgSize = " << msgSize << std::endl;
#endif

      // Communicate the supernode
      MPI_Bcast( &msgSize, 1, MPI_INT, mpirank, comm );
      MPI_Bcast( snodeIcomm.front(), msgSize, MPI_CHAR, mpirank, comm );
      // Copy the data from buffer to snode
      LIBCHOLESKY::Deserialize( snodeIcomm.front(), snode );
    } // if owning the supernode
    else{
      // Receive the supernode
      Int rootRank = ( iSuper % mpisize );
      Int msgSize;
      MPI_Bcast( &msgSize, 1, MPI_INT, rootRank, comm );
      buffer.resize(msgSize);
      MPI_Bcast( &buffer[0], msgSize, MPI_CHAR, rootRank, comm );
      LIBCHOLESKY::Deserialize( &buffer[0], snode );
    } // if not owning the supernode but is in the receiving column

#if ( _DEBUGlevel_ >= 1 )
    statusOFS << "All communication is finished." << std::endl;
#endif

    // Local operation from now
    if( mpicol == ( iSuper % npcol ) ){
      Int jb = iSuper / npcol;
      std::vector<LBlock<T> >& Lcol = PMat.L(jb);
      std::set<Int> blkSet;
      Int superSize = snode.Size();

      // Count the number of blocks in the supernode belonging to this
      // processor
      for( Int blkidx = 0; blkidx < snode.NZBlockCnt(); blkidx++ ){
        LIBCHOLESKY::NZBlockDesc desc = snode.GetNZBlockDesc( blkidx );
        Int nrows = snode.NRows(blkidx);
        Int firstRow = desc.GIndex - 1;
        Int lastRow = firstRow + nrows - 1;
#if ( _DEBUGlevel_ >= 1 )
    statusOFS << "firstRow = " << firstRow << ", lastRow = " << lastRow << std::endl;
#endif
        for( Int i = superIdx(firstRow); i <= superIdx(lastRow); i++ ){
          if( mpirow == ( i % nprow ) ){
            blkSet.insert( i );
          } // if the current processor is in the right processor row
        }
      } // for ( blkidx )
      
      Int numBlkLocal = blkSet.size();
      std::vector<Int> blkVec;
      blkVec.insert( blkVec.end(), blkSet.begin(), blkSet.end() );
      Lcol.resize( blkVec.size() );

      // Allocate the nonzero rows and nzvals 
#if ( _DEBUGlevel_ >= 1 )
      statusOFS << "Lcol.size = " << Lcol.size() << std::endl;
      statusOFS << "blkSet.size = " << blkSet.size() << std::endl;
      statusOFS << "blkVec = " << blkVec << std::endl;
#endif

      std::vector<std::vector<Int> > rowsBlk( Lcol.size() );
      std::vector<std::vector<T> > nzvalBlk( Lcol.size() );

      for( Int blkidx = 0; blkidx < snode.NZBlockCnt(); blkidx++ ){
        LIBCHOLESKY::NZBlockDesc desc = snode.GetNZBlockDesc( blkidx );
        Int nrows = snode.NRows(blkidx);
        Int firstRow = desc.GIndex - 1;
        Int lastRow = firstRow + nrows - 1;
        T* nzval = snode.GetNZval( desc.Offset );
        std::vector<Int>::iterator vi;
        Int pos;
        for( Int i = firstRow; i <= lastRow; i++ ){
          vi = find( blkVec.begin(), blkVec.end(), superIdx(i) );
          if( vi != blkVec.end() ){
            pos = vi - blkVec.begin();
            rowsBlk[pos].push_back(i);
            nzvalBlk[pos].insert(
                nzvalBlk[pos].end(),
                nzval, nzval + superSize ); 
          }
          nzval += superSize;
        }
      } // for ( blkidx )
      
      // Save the information to Lcol
      for ( Int iblk = 0; iblk < Lcol.size(); iblk++ ){
        std::vector<Int>& rows = rowsBlk[iblk];
        std::vector<T>& nzval = nzvalBlk[iblk];
        LBlock<T>& LB = Lcol[iblk];
        LB.blockIdx = blkVec[iblk];
        LB.numRow = rows.size();
        LB.numCol = superSize;
        LB.rows = IntNumVec( LB.numRow, true, &rows[0] );
        if( LB.numRow * LB.numCol != nzval.size() ){
          std::ostringstream msg;
          msg << "message size does not match for the blockIdx " << LB.blockIdx << std::endl
            << "LB.numRow * LB.numCol = " << LB.numRow * LB.numCol << std::endl
            << "nzval.size            = " << nzval.size() << std::endl;
          throw std::runtime_error( msg.str().c_str() );
        }
        // Convert the row major format to column major format
        Transpose( NumMat<T>( LB.numCol, LB.numRow, true,
            &nzval[0] ), LB.nzval );
      }
    } // if the current processor is in the right processor column

    // Set the MPI Barrier
    MPI_Barrier( comm );
  }


  PMatrixLtoU( PMat );

#ifndef _RELEASE_
	PopCallStack();
#endif
}  // -----  end of NGCHOLMatrixToPMatrix ----- 

template<typename T> void PMatrixLtoU( PMatrix<T>& PMat )
{

#ifndef _RELEASE_
	PushCallStack("PMatrixLtoU");
#endif
  //Send L to U
  Int mpirank, mpisize;
  const GridType *g = PMat.Grid();

  // FIXME Check PMatrix and SupernodalMatrix has the same communicator
  MPI_Comm comm = g->comm;
  MPI_Comm colComm = g->colComm;
  MPI_Comm rowComm = g->rowComm;
  MPI_Comm_size(comm, &mpisize); 
  MPI_Comm_rank(comm, &mpirank);

  Int nprow = g->numProcRow, npcol = g->numProcCol;


  Int mpirow = mpirank / npcol;  
  Int mpicol = mpirank % npcol;

  Int numSuper = PMat.NumSuper();
  for( Int ksup = 0; ksup < numSuper; ksup++ ){
#if ( _DEBUGlevel_ >= 1 )
statusOFS<<"----------------------- "<< ksup<<std::endl;
#endif
     //If I'm in the supernodal column
     std::vector<Int> all_proc_list;
     std::vector<Int> all_blocks_cnt;
     std::vector<Int> sizes;
     std::vector<Int> displs;

     std::vector<Int> sender_proc;
     std::vector<std::list<LBlock<T> * > > blocks_to_receiver;

     std::vector<Int> receiver_proc;
     if( mpicol == ( ksup % npcol ) ){
        Int jb = ksup / npcol;
        std::vector<LBlock<T> >& Lcol = PMat.L(jb);
        Int root = (ksup % nprow);
        //compute the list of receiving processors based on my local structure
        std::set<Int> proc_set;
        blocks_to_receiver.resize(npcol);
        Int startBlk = mpirow==root?1:0;
        for ( Int iblk = startBlk; iblk < Lcol.size(); iblk++ ){
          LBlock<T>& LB = Lcol[iblk];
          //column of the target processor
          Int snode_idx = LB.blockIdx;
          Int tgt_pcol =  snode_idx % npcol;
          blocks_to_receiver[tgt_pcol].push_back(&LB);
          proc_set.insert(tgt_pcol);
        }
        //Insert the set into the vector to be able to send it
        receiver_proc.insert(receiver_proc.begin(),proc_set.begin(),proc_set.end());
        
        //Now do a gatherv on the root
        mpi::Gatherv(receiver_proc,all_proc_list,sizes,displs,root, colComm);

        //Do a gatherv of the local blocks to each processors
        std::vector<Int> blocks_cnt(receiver_proc.size());
        for(Int j = 0; j< receiver_proc.size();++j){
          Int pcol = receiver_proc[j];
          std::list<LBlock<T> *> & blocks = blocks_to_receiver[pcol];
          blocks_cnt[j] = blocks.size();
        }
        mpi::Gatherv(blocks_cnt,all_blocks_cnt,root, colComm);
  

        //On the root, convert from a sender array to a receiver array
        if(mpirow == root){
          //sender_proc[i] contains the set of sender to processor column i
          std::vector<std::set<Int> > sender_procs(npcol);
          std::vector<Int> urow_sizes(npcol,0);      
          for(Int prow = 0; prow < nprow; ++prow){
            Int * recv_list = &all_proc_list[displs[prow]];
            Int * recv_blocks = &all_blocks_cnt[displs[prow]];
            Int size = sizes[prow];
            for(Int i = 0; i<size;++i){
              Int pcol = recv_list[i];
              sender_procs[pcol].insert(prow);
              Int ucol_contrib = recv_blocks[i];
              urow_sizes[pcol]+=ucol_contrib;
            }
          }
          //now prepare the data structures for a scatterv along the rows
          all_blocks_cnt = urow_sizes;
          all_proc_list.clear();
          sizes.resize(npcol);
          displs.resize(npcol);
          Int totalsize = 0;
          for(Int pcol = 0; pcol < npcol; ++pcol){
            sizes[pcol] = sender_procs[pcol].size()*sizeof(Int);
            displs[pcol] = totalsize;
            totalsize += sizes[pcol];
          }
          //put the senders in the all_proc_list_array
          all_proc_list.reserve(totalsize / sizeof(Int) );
          for(Int pcol = 0; pcol < npcol; ++pcol){
            all_proc_list.insert(all_proc_list.end(),sender_procs[pcol].begin(),sender_procs[pcol].end());
          }
        }
     }

     //If I'm in the supernodal row
     if( mpirow == ( ksup % nprow ) ){
        Int root = (ksup % npcol);
        //scatter the sizes
        Int localSize = 0;
        MPI_Scatter(mpicol==root?&sizes[0]:NULL,sizeof(Int),MPI_BYTE,
                      &localSize,sizeof(Int),MPI_BYTE, root, rowComm);
        sender_proc.resize(localSize / sizeof(Int));
        //Now do the scatterv; 
        if(mpicol==root){
          MPI_Scatterv(&all_proc_list[0],&sizes[0],&displs[0],MPI_BYTE,
              &sender_proc[0],localSize,MPI_BYTE, root, rowComm);
        }
        else{
          MPI_Scatterv(NULL,NULL,NULL,MPI_BYTE,
              &sender_proc[0],localSize,MPI_BYTE, root, rowComm);
        }
        
        Int urowSize = 0;
        MPI_Scatter(mpicol==root?&all_blocks_cnt[0]:NULL,sizeof(Int),MPI_BYTE,
                      &urowSize,sizeof(Int),MPI_BYTE, root, rowComm);
        //Resize Urow
        Int ib = ksup / nprow;
        std::vector<UBlock<T> >& Urow = PMat.U(ib);
        Urow.resize(urowSize);


        //At this point we have both a sender AND a receiver list
        //and Urows are resized
     }

     //Communicate this supernode
     //If I'm a sender
     if( mpicol == ( ksup % npcol ) ){
        std::vector<Int> mask( LBlockMask::TOTAL_NUMBER, 1 );
        //for each target col
        for(Int pcol = 0; pcol < npcol; pcol++){
          Int pnum = (ksup % nprow)*npcol + pcol;
          if(pnum!=mpirank){
            //Serialize everything in the blocks_to_receiver list
            std::list<LBlock<T> *> & blocks = blocks_to_receiver[pcol];
            if(blocks.size()>0){
              std::stringstream sstm;
              Int numLBlocks = blocks.size();
#if ( _DEBUGlevel_ >= 1 )
statusOFS<<"Sending "<<numLBlocks<<" LBlocks"<<std::endl;
#endif
              serialize( numLBlocks , sstm, NO_MASK);
              typename std::list<LBlock<T> *>::iterator it;
              for(it = blocks.begin(); 
                  it!=blocks.end();++it){
                LBlock<T> & LB = *(*it);
#if ( _DEBUGlevel_ >= 1 )
statusOFS<<"Sent LB: "<<LB<<std::endl;
#endif
                serialize(LB, sstm,mask);
              }
              mpi::Send(sstm, pnum, PMat.IdxToTag(ksup,PMatrix<T>::SELINV_TAG_L_SIZE),
                  PMat.IdxToTag(ksup,PMatrix<T>::SELINV_TAG_L_CONTENT), comm);
            }
          }
        }
     }
     //If I'm a receiver 
     if( mpirow == ( ksup % nprow ) ){
      std::vector<Int> mask( LBlockMask::TOTAL_NUMBER, 1 );
        Int ib = ksup / nprow;
        std::vector<UBlock<T> >& Urow = PMat.U(ib);

        //for each target row
        Int idx = 0;
        for(Int i = 0; i < sender_proc.size(); ++i){

          Int prow = sender_proc[i];
          Int pnum = (prow)*npcol + (ksup % npcol);
          if(pnum != mpirank){
            std::stringstream sstm;
            mpi::Recv(sstm,pnum,PMat.IdxToTag(ksup,PMatrix<T>::SELINV_TAG_L_SIZE),
                PMat.IdxToTag(ksup,PMatrix<T>::SELINV_TAG_L_CONTENT), comm);

            //now deserialize and put everything in U
            Int numLBlocks = 0;
            deserialize(numLBlocks, sstm, NO_MASK);
#if ( _DEBUGlevel_ >= 1 )
statusOFS<<"Receiving "<<numLBlocks<<" LBlocks"<<std::endl;
#endif
            for(Int i = 0; i<numLBlocks;++i){
              LBlock<T> LB;
              deserialize(LB, sstm, mask);

#if ( _DEBUGlevel_ >= 1 )
statusOFS<<"Received LB: "<<LB<<std::endl;
#endif
              //put this LBlock in the appropriate UBlock
                UBlock<T> & UB = Urow[idx];

                UB.blockIdx = LB.blockIdx;
                UB.numCol = LB.numRow;
                UB.numRow = LB.numCol;
                UB.cols = LB.rows;
                Transpose(LB.nzval,UB.nzval);
                ++idx; 
            }
          }
          else{
            Int pcol = mpicol ;
            //do the copy locally
            Int ib = ksup / nprow;
            std::vector<UBlock<T> >& Urow = PMat.U(ib);
            //Serialize everything in the blocks_to_receiver list
            std::list<LBlock<T> *> & blocks = blocks_to_receiver[pcol];
            if(blocks.size()>0){
              typename std::list<LBlock<T> *>::iterator it;
              for(it = blocks.begin(); 
                  it!=blocks.end();++it){
                LBlock<T> & LB = *(*it);
                UBlock<T> & UB = Urow[idx];

                UB.blockIdx = LB.blockIdx;
                UB.numCol = LB.numRow;
                UB.numRow = LB.numCol;
                UB.cols = LB.rows;
                Transpose(LB.nzval,UB.nzval);
                ++idx; 
              }
            }
          }

          
        }
     }
  }
#ifndef _RELEASE_
	PopCallStack();
#endif

}  // -----  end of PMatrixLToU ----- 


}


#endif //_PEXSI_NGCHOL_INTERF_IMPL_HPP_

