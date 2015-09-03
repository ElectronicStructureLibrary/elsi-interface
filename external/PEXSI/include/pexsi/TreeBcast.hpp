#ifndef _PEXSI_TREE_HPP_
#define _PEXSI_TREE_HPP_

#include "pexsi/environment.hpp"
#include "pexsi/timer.h"

#include <vector>
#include <algorithm>
#include <string>
#include <random>

// options to switch from a flat bcast/reduce tree to a binary tree

#ifndef FTREE_LIMIT
#define FTREE_LIMIT 24
#endif



namespace PEXSI{

class TreeBcast{
  protected:
    Int myRoot_;
    MPI_Comm comm_;
    vector<Int> myDests_;
    Int myRank_;
    Int msgSize_;
    bool isReady_;
    Int mainRoot_;
    Int tag_;
    Int numRecv_;


#ifdef COMM_PROFILE
protected:
    Int myGRank_;
    vector<int> Granks_;
public:
    void SetGlobalComm(const MPI_Comm & pGComm){
      MPI_Comm_rank(pGComm,&myGRank_);
      MPI_Group group2 = MPI_GROUP_NULL;
      MPI_Comm_group(pGComm, &group2);
      MPI_Group group1 = MPI_GROUP_NULL;
      MPI_Comm_group(comm_, &group1);

      Int size;
      MPI_Comm_size(comm_,&size);
      Granks_.resize(size);
      vector<int> Lranks(size);
      for(int i = 0; i<size;++i){Lranks[i]=i;}
      MPI_Group_translate_ranks(group1, size, &Lranks[0],group2, &Granks_[0]);
    }
#endif


  
    virtual void buildTree(Int * ranks, Int rank_cnt)=0;
  public:
    TreeBcast(){
      comm_ = MPI_COMM_WORLD;
      myRank_=-1;
      myRoot_ = -1; 
      msgSize_ = -1;
      numRecv_ = -1;
      tag_=-1;
      mainRoot_=-1;
      isReady_ = false;
    }


    TreeBcast(const MPI_Comm & pComm, Int * ranks, Int rank_cnt,Int msgSize){
      comm_ = pComm;
      MPI_Comm_rank(comm_,&myRank_);
      myRoot_ = -1; 
      msgSize_ = msgSize;

      numRecv_ = 0;
      tag_=-1;
      mainRoot_=ranks[0];
      isReady_ = false;
    }

    TreeBcast(const TreeBcast & Tree){
      Copy(Tree); 
    }
 
    virtual void Copy(const TreeBcast & Tree){
      comm_ = Tree.comm_;
      myRank_ = Tree.myRank_;
      myRoot_ = Tree.myRoot_; 
      msgSize_ = Tree.msgSize_;

      numRecv_ = Tree.numRecv_;
      tag_= Tree.tag_;
      mainRoot_= Tree.mainRoot_;
      isReady_ = Tree.isReady_;
      myDests_ = Tree.myDests_;
    }

    virtual TreeBcast * clone() const = 0; 


    
    static TreeBcast * Create(const MPI_Comm & pComm, Int * ranks, Int rank_cnt, Int msgSize,double rseed);

    virtual inline Int GetNumRecvMsg(){return numRecv_;}
    virtual inline Int GetNumMsgToRecv(){return 1;}
    inline void SetDataReady(bool rdy){ isReady_=rdy;}
    inline void SetTag(Int tag){ tag_ = tag;}


    Int * GetDests(){ return &myDests_[0];}
    Int GetDest(Int i){ return myDests_[i];}
    Int GetDestCount(){ return myDests_.size();}
    Int GetRoot(){ return myRoot_;}
    Int GetMsgSize(){ return msgSize_;}

    void ForwardMessage( char * data, size_t size, int tag, MPI_Request * requests ){
                  for( Int idxRecv = 0; idxRecv < myDests_.size(); ++idxRecv ){
                    Int iProc = myDests_[idxRecv];
                    // Use Isend to send to multiple targets
                    MPI_Isend( data, size, MPI_BYTE, 
                        iProc, tag,comm_, &requests[2*iProc+1] );

#ifdef COMM_PROFILE
      PROFILE_COMM(myGRank_,Granks_[iProc],tag,msgSize_);
#endif
                  } // for (iProc)
    }


};

class FTreeBcast: public TreeBcast{
  protected:
    virtual void buildTree(Int * ranks, Int rank_cnt){

      Int idxStart = 0;
      Int idxEnd = rank_cnt;



      myRoot_ = ranks[0];

      if(myRank_==myRoot_){
        myDests_.insert(myDests_.begin(),&ranks[1],&ranks[0]+rank_cnt);
      }

#if ( _DEBUGlevel_ >= 1 )
      statusOFS<<"My root is "<<myRoot_<<std::endl;
      statusOFS<<"My dests are ";
      for(int i =0;i<myDests_.size();++i){statusOFS<<myDests_[i]<<" ";}
      statusOFS<<std::endl;
#endif
    }



  public:
    FTreeBcast(const MPI_Comm & pComm, Int * ranks, Int rank_cnt, Int msgSize):TreeBcast(pComm,ranks,rank_cnt,msgSize){
      //build the binary tree;
      buildTree(ranks,rank_cnt);
    }


    virtual FTreeBcast * clone() const{
      FTreeBcast * out = new FTreeBcast(*this);
      return out;
    } 
};



class BTreeBcast: public TreeBcast{
  protected:
////  virtual void buildTree(Int * ranks, Int rank_cnt){
////          Int numLevel = floor(log2(rank_cnt));
////          Int numRoots = 0;
////          for(Int level=0;level<numLevel;++level){
////            numRoots = std::min( rank_cnt, numRoots + (Int)pow(2,level));
////            Int numNextRoots = std::min(rank_cnt,numRoots + (Int)pow(2,(level+1)));
////            Int numReceivers = numNextRoots - numRoots;
////            for(Int ip = 0; ip<numRoots;++ip){
////              Int p = ranks[ip];
////              for(Int ir = ip; ir<numReceivers;ir+=numRoots){
////                Int r = ranks[numRoots+ir];
////                if(r==myRank_){
////                  myRoot_ = p;
////                }
////
////                if(p==myRank_){
////                  myDests_.push_back(r);
////                }
////              }
////            }
////          }
////  }
////
    virtual void buildTree(Int * ranks, Int rank_cnt){

      Int idxStart = 0;
      Int idxEnd = rank_cnt;



      Int prevRoot = ranks[0];
      while(idxStart<idxEnd){
        Int curRoot = ranks[idxStart];
        Int listSize = idxEnd - idxStart;

        if(listSize == 1){
          if(curRoot == myRank_){
            myRoot_ = prevRoot;
            break;
          }
        }
        else{
          Int halfList = floor(ceil(double(listSize) / 2.0));
          Int idxStartL = idxStart+1;
          Int idxStartH = idxStart+halfList;

          if(curRoot == myRank_){
            if ((idxEnd - idxStartH) > 0 && (idxStartH - idxStartL)>0){
              Int childL = ranks[idxStartL];
              Int childR = ranks[idxStartH];

              myDests_.push_back(childL);
              myDests_.push_back(childR);
            }
            else if ((idxEnd - idxStartH) > 0){
              Int childR = ranks[idxStartH];
              myDests_.push_back(childR);
            }
            else{
              Int childL = ranks[idxStartL];
              myDests_.push_back(childL);
            }
            myRoot_ = prevRoot;
            break;
          } 

          if( myRank_ < ranks[idxStartH]){
            idxStart = idxStartL;
            idxEnd = idxStartH;
          }
          else{
            idxStart = idxStartH;
          }
          prevRoot = curRoot;
        }

      }

#if ( _DEBUGlevel_ >= 1 )
      statusOFS<<"My root is "<<myRoot_<<std::endl;
      statusOFS<<"My dests are ";
      for(int i =0;i<myDests_.size();++i){statusOFS<<myDests_[i]<<" ";}
      statusOFS<<std::endl;
#endif
    }



  public:
    BTreeBcast(const MPI_Comm & pComm, Int * ranks, Int rank_cnt, Int msgSize):TreeBcast(pComm,ranks,rank_cnt,msgSize){
      //build the binary tree;
      buildTree(ranks,rank_cnt);
    }

    virtual BTreeBcast * clone() const{
      BTreeBcast * out = new BTreeBcast(*this);
      return out;
    }

};



class ModBTreeBcast: public TreeBcast{
  protected:
    double rseed_;

    virtual void buildTree(Int * ranks, Int rank_cnt){

      Int idxStart = 0;
      Int idxEnd = rank_cnt;

      //sort the ranks with the modulo like operation
      if(rank_cnt>1){
        //Int new_idx = (int)((rand()+1.0) * (double)rank_cnt / ((double)RAND_MAX+1.0));

//        srand(ranks[0]+rank_cnt);
        Int new_idx = (int)((rank_cnt - 0) * ( (double)this->rseed_ / (double)RAND_MAX ) + 0);// (this->rseed_)%(rank_cnt-1)+1;
//        statusOFS<<new_idx<<endl;



        Int * new_start = &ranks[new_idx];

        //for(int i =0;i<rank_cnt;++i){statusOFS<<ranks[i]<<" ";} statusOFS<<std::endl;

//        Int * new_start = std::lower_bound(&ranks[1],&ranks[0]+rank_cnt,ranks[0]);
        //just swap the two chunks   r[0] | r[1] --- r[new_start-1] | r[new_start] --- r[end]
        // becomes                   r[0] | r[new_start] --- r[end] | r[1] --- r[new_start-1] 
        std::rotate(&ranks[1], new_start, &ranks[0]+rank_cnt);

        //for(int i =0;i<rank_cnt;++i){statusOFS<<ranks[i]<<" ";} statusOFS<<std::endl;
      }

      Int prevRoot = ranks[0];
      while(idxStart<idxEnd){
        Int curRoot = ranks[idxStart];
        Int listSize = idxEnd - idxStart;

        if(listSize == 1){
          if(curRoot == myRank_){
            myRoot_ = prevRoot;
            break;
          }
        }
        else{
          Int halfList = floor(ceil(double(listSize) / 2.0));
          Int idxStartL = idxStart+1;
          Int idxStartH = idxStart+halfList;

          if(curRoot == myRank_){
            if ((idxEnd - idxStartH) > 0 && (idxStartH - idxStartL)>0){
              Int childL = ranks[idxStartL];
              Int childR = ranks[idxStartH];

              myDests_.push_back(childL);
              myDests_.push_back(childR);
            }
            else if ((idxEnd - idxStartH) > 0){
              Int childR = ranks[idxStartH];
              myDests_.push_back(childR);
            }
            else{
              Int childL = ranks[idxStartL];
              myDests_.push_back(childL);
            }
            myRoot_ = prevRoot;
            break;
          } 

          //not true anymore ?
          //first half to 
TIMER_START(FIND_RANK);
          Int * pos = std::find(&ranks[idxStartL], &ranks[idxStartH], myRank_);
TIMER_STOP(FIND_RANK);
          if( pos != &ranks[idxStartH]){
            idxStart = idxStartL;
            idxEnd = idxStartH;
          }
          else{
            idxStart = idxStartH;
          }
          prevRoot = curRoot;
        }

      }

#if ( _DEBUGlevel_ >= 1 ) || defined(REDUCE_VERBOSE)
      statusOFS<<"My root is "<<myRoot_<<std::endl;
      statusOFS<<"My dests are ";
      for(int i =0;i<myDests_.size();++i){statusOFS<<myDests_[i]<<" ";}
      statusOFS<<std::endl;
#endif
    }



  public:
    ModBTreeBcast(const MPI_Comm & pComm, Int * ranks, Int rank_cnt, Int msgSize, double rseed):TreeBcast(pComm,ranks,rank_cnt,msgSize){
      //build the binary tree;
      rseed_ = rseed;
      buildTree(ranks,rank_cnt);
    }

    virtual void Copy(const ModBTreeBcast & Tree){
      comm_ = Tree.comm_;
      myRank_ = Tree.myRank_;
      myRoot_ = Tree.myRoot_; 
      msgSize_ = Tree.msgSize_;

      numRecv_ = Tree.numRecv_;
      tag_= Tree.tag_;
      mainRoot_= Tree.mainRoot_;
      isReady_ = Tree.isReady_;
      myDests_ = Tree.myDests_;

      rseed_ = Tree.rseed_;
      myRank_ = Tree.myRank_;
      myRoot_ = Tree.myRoot_; 
      msgSize_ = Tree.msgSize_;

      numRecv_ = Tree.numRecv_;
      tag_= Tree.tag_;
      mainRoot_= Tree.mainRoot_;
      isReady_ = Tree.isReady_;
      myDests_ = Tree.myDests_;
    }
 
    virtual ModBTreeBcast * clone() const{
      ModBTreeBcast * out = new ModBTreeBcast(*this);
      return out;
    }

};


class RandBTreeBcast: public TreeBcast{
  protected:
    virtual void buildTree(Int * ranks, Int rank_cnt){

      Int idxStart = 0;
      Int idxEnd = rank_cnt;

      //random permute ranks
      if(rank_cnt>1){
        for(int i =0;i<rank_cnt;++i){statusOFS<<ranks[i]<<" ";} statusOFS<<std::endl;
        srand(ranks[0]);
        std::random_shuffle(&ranks[1],&ranks[0]+rank_cnt);
        for(int i =0;i<rank_cnt;++i){statusOFS<<ranks[i]<<" ";} statusOFS<<std::endl;

      }

      Int prevRoot = ranks[0];
      while(idxStart<idxEnd){
        Int curRoot = ranks[idxStart];
        Int listSize = idxEnd - idxStart;

        if(listSize == 1){
          if(curRoot == myRank_){
            myRoot_ = prevRoot;
            break;
          }
        }
        else{
          Int halfList = floor(ceil(double(listSize) / 2.0));
          Int idxStartL = idxStart+1;
          Int idxStartH = idxStart+halfList;

          if(curRoot == myRank_){
            if ((idxEnd - idxStartH) > 0 && (idxStartH - idxStartL)>0){
              Int childL = ranks[idxStartL];
              Int childR = ranks[idxStartH];

              myDests_.push_back(childL);
              myDests_.push_back(childR);
            }
            else if ((idxEnd - idxStartH) > 0){
              Int childR = ranks[idxStartH];
              myDests_.push_back(childR);
            }
            else{
              Int childL = ranks[idxStartL];
              myDests_.push_back(childL);
            }
            myRoot_ = prevRoot;
            break;
          } 

          //not true anymore ?
          //first half to 
          Int * pos = std::find(&ranks[idxStartL], &ranks[idxStartH], myRank_);
          if( pos != &ranks[idxStartH]){
            idxStart = idxStartL;
            idxEnd = idxStartH;
          }
          else{
            idxStart = idxStartH;
          }
          prevRoot = curRoot;
        }

      }

#if ( _DEBUGlevel_ >= 1 ) || defined(REDUCE_VERBOSE)
      statusOFS<<"My root is "<<myRoot_<<std::endl;
      statusOFS<<"My dests are ";
      for(int i =0;i<myDests_.size();++i){statusOFS<<myDests_[i]<<" ";}
      statusOFS<<std::endl;
#endif
    }



  public:
    RandBTreeBcast(const MPI_Comm & pComm, Int * ranks, Int rank_cnt, Int msgSize):TreeBcast(pComm,ranks,rank_cnt,msgSize){
      //build the binary tree;
      buildTree(ranks,rank_cnt);
    }

    virtual RandBTreeBcast * clone() const{
      RandBTreeBcast * out = new RandBTreeBcast(*this);
      return out;
    }

};




class PalmTreeBcast: public TreeBcast{
  protected:
  virtual void buildTree(Int * ranks, Int rank_cnt){
          Int numLevel = floor(log2(rank_cnt));
          Int numRoots = 0;
          for(Int level=0;level<numLevel;++level){
            numRoots = std::min( rank_cnt, numRoots + (Int)pow(2,level));
            Int numNextRoots = std::min(rank_cnt,numRoots + (Int)pow(2,(level+1)));
            Int numReceivers = numNextRoots - numRoots;
            for(Int ip = 0; ip<numRoots;++ip){
              Int p = ranks[ip];
              for(Int ir = ip; ir<numReceivers;ir+=numRoots){
                Int r = ranks[numRoots+ir];
                if(r==myRank_){
                  myRoot_ = p;
                }

                if(p==myRank_){
                  myDests_.push_back(r);
                }
              }
            }
          }

#if ( _DEBUGlevel_ >= 1 )
      statusOFS<<"My root is "<<myRoot_<<std::endl;
      statusOFS<<"My dests are ";
      for(int i =0;i<myDests_.size();++i){statusOFS<<myDests_[i]<<" ";}
      statusOFS<<std::endl;
#endif
    }



  public:
    PalmTreeBcast(const MPI_Comm & pComm, Int * ranks, Int rank_cnt, Int msgSize):TreeBcast(pComm,ranks,rank_cnt,msgSize){
      //build the binary tree;
      buildTree(ranks,rank_cnt);
    }

    virtual PalmTreeBcast * clone() const{
      PalmTreeBcast * out = new PalmTreeBcast(*this);
      return out;
    }

};


template< typename T>
class TreeReduce: public TreeBcast{
  protected:

    T * myData_;
    MPI_Request sendRequest_;

    //char * myLocalBuffer_;
    NumVec<char> myLocalBuffer_;
    //char * myRecvBuffers_;
    NumVec<char> myRecvBuffers_;
    NumVec<T *> remoteData_;
    NumVec<MPI_Request> myRequests_;
    NumVec<MPI_Status> myStatuses_;
    NumVec<int> recvIdx_;

    bool fwded_;
    bool isAllocated_;
    Int numRecvPosted_;

  public:
    TreeReduce(const MPI_Comm & pComm, Int * ranks, Int rank_cnt, Int msgSize):TreeBcast(pComm,ranks,rank_cnt,msgSize){
      myData_ = NULL;
      sendRequest_ = MPI_REQUEST_NULL;
      fwded_=false;
      isAllocated_=false;
      numRecvPosted_= 0;
    }


    virtual TreeReduce * clone() const = 0; 

    TreeReduce(const TreeReduce & Tree){
      this->Copy(Tree);
    }

    virtual void Copy(const TreeReduce & Tree){
      this->comm_ = Tree.comm_;
      this->myRank_ = Tree.myRank_;
      this->myRoot_ = Tree.myRoot_; 
      this->msgSize_ = Tree.msgSize_;

      this->numRecv_ = Tree.numRecv_;
      this->tag_= Tree.tag_;
      this->mainRoot_= Tree.mainRoot_;
      this->isReady_ = Tree.isReady_;
      this->myDests_ = Tree.myDests_;


      this->myData_ = Tree.myData_;
      this->sendRequest_ = Tree.sendRequest_;
      this->fwded_= Tree.fwded_;
      this->isAllocated_= Tree.isAllocated_;
      this->numRecvPosted_= Tree.numRecvPosted_;

      this->myLocalBuffer_ = Tree.myLocalBuffer_;
      this->myRecvBuffers_ = Tree.myRecvBuffers_;
      this->remoteData_ = Tree.remoteData_;
      this->myRequests_ = Tree.myRequests_;
      this->myStatuses_ = Tree.myStatuses_;
      this->recvIdx_ = Tree.recvIdx_;
    }
 


    bool IsAllocated(){return isAllocated_;}

    virtual ~TreeReduce(){
      CleanupBuffers();
    }


    static TreeReduce<T> * Create(const MPI_Comm & pComm, Int * ranks, Int rank_cnt, Int msgSize,double rseed);

    virtual inline Int GetNumMsgToRecv(){return GetDestCount();}

    virtual void AllocRecvBuffers(){
      remoteData_.Resize(GetDestCount());
      //SetValue(remoteData_,(T*)NULL);

      //assert(myRecvBuffers_==NULL);
      //myRecvBuffers_ = new char[GetDestCount()*msgSize_];


      myRecvBuffers_.Resize(GetDestCount()*msgSize_);
      //SetValue(myRecvBuffers_,(char)0);

      for( Int idxRecv = 0; idxRecv < GetDestCount(); ++idxRecv ){
        remoteData_[idxRecv] = (T*)&(myRecvBuffers_[idxRecv*msgSize_]);
        //Int nelem = msgSize_ / sizeof(T);
        //std::fill(remoteData_[idxRecv],remoteData_[idxRecv]+nelem,ZERO<T>());
      }

      myRequests_.Resize(GetDestCount());
      SetValue(myRequests_,MPI_REQUEST_NULL);
      myStatuses_.Resize(GetDestCount());
      recvIdx_.Resize(GetDestCount());

      sendRequest_ = MPI_REQUEST_NULL;

      isAllocated_ = true;
    }

    void CleanupBuffers(){
      myLocalBuffer_.Clear();
//      if(myLocalBuffer_!=NULL){
//        delete []myLocalBuffer_;
//        myLocalBuffer_=NULL;
//      }


      remoteData_.Clear();
//      myRecvBuffers_.Clear();
//      if(myRecvBuffers_!=NULL){
//        delete []myRecvBuffers_;
//        myRecvBuffers_=NULL;
//      }


      myRequests_.Clear();
      myStatuses_.Clear();
      recvIdx_.Clear();


//              if(myLocalBuffer_!=NULL){
//                delete [] myLocalBuffer_;
//              }
//              myLocalBuffer_=NULL;


      myData_ = NULL;
      sendRequest_ = MPI_REQUEST_NULL;
      fwded_=false;
      isAllocated_=false;
      isReady_=false;
      numRecv_ = 0;
      numRecvPosted_= 0;
    }


    void SetLocalBuffer(T * locBuffer){
      if(myData_!=NULL && myData_!=locBuffer){
        blas::Axpy(msgSize_/sizeof(T), ONE<T>(), myData_, 1, locBuffer, 1 );
        myLocalBuffer_.Clear(); 
      }

      myData_ = locBuffer;
    }

    inline bool AccumulationDone(){
      if(myRank_==myRoot_ && isAllocated_){
        isReady_=true;
      }
      return isReady_ && (numRecv_ == GetDestCount());
    }


    inline bool IsDone(){
      if(myRank_==myRoot_ && isAllocated_){
        isReady_=true;
      }

      bool retVal = AccumulationDone();
      if(myRoot_ != myRank_ && !fwded_){
        retVal = false;
      }

      if (retVal && myRoot_ != myRank_ && fwded_){
        //test the send request
        int flag = 0;
        MPI_Test(&sendRequest_,&flag,MPI_STATUS_IGNORE);
        retVal = flag==1;
      }
      
      return retVal;
    }

    //async wait and forward
    virtual bool Progress(){
      if(!isAllocated_){
        return true;
      }

      if(myRank_==myRoot_ && isAllocated_){
        isReady_=true;
      }

//      if(this->numRecvPosted_==0){ 
//        this->PostFirstRecv();
//      }

      bool retVal = AccumulationDone();
      if(isReady_ && !retVal){

        //assert(isAllocated_);

        //mpi_test_some on my requests
        int recvCount = -1;
        int reqCnt = GetDestCount();
        assert(reqCnt == myRequests_.m());
        MPI_Testsome(reqCnt,&myRequests_[0],&recvCount,&recvIdx_[0],&myStatuses_[0]);
        //if something has been received, accumulate and potentially forward it
        for(Int i = 0;i<recvCount;++i ){
          Int idx = recvIdx_[i];

         if(idx!=MPI_UNDEFINED){

          Int size = 0;
          MPI_Get_count(&myStatuses_[i], MPI_BYTE, &size);


#if ( _DEBUGlevel_ >= 1 ) || defined(REDUCE_VERBOSE)
        statusOFS<<myRank_<<" RECVD from "<<myStatuses_[i].MPI_SOURCE<<" on tag "<<tag_<<std::endl;
#endif
          if(size>0){
            //If myData is 0, allocate to the size of what has been received
            if(myData_==NULL){
              //assert(size==msgSize_);
              myLocalBuffer_.Resize(msgSize_);
    
              myData_ = (T*)&myLocalBuffer_[0];
              Int nelem = +msgSize_/sizeof(T);
              std::fill(myData_,myData_+nelem,ZERO<T>());
            }

            Reduce(idx,i);

          }

          numRecv_++;
          //MPI_Request_free(&myRequests_[idx]);
          }
        }

      }
      else if (isReady_ && sendRequest_ == MPI_REQUEST_NULL && myRoot_ != myRank_ && !fwded_){
        //free the unnecessary arrays
        myRecvBuffers_.Clear();
        myRequests_.Clear();
        myStatuses_.Clear();
        recvIdx_.Clear();

        //assert(isAllocated_);

        //Forward
        Forward();
        retVal = false;
      }
      else{
        retVal = IsDone();
        if(retVal){
          //free the unnecessary arrays
          myRecvBuffers_.Clear();
          myRequests_.Clear();
          myStatuses_.Clear();
          recvIdx_.Clear();
        }
      }

      return retVal;
    }

    //blocking wait
    void Wait(){
      if(isAllocated_){
        while(!Progress());
      }
    }

    T * GetLocalBuffer(){
       return myData_;
    }



    void CopyLocalBuffer(T* destBuffer){
       std::copy((char*)myData_,(char*)myData_+GetMsgSize(),(char*)destBuffer);
    }


    virtual void PostFirstRecv()
    {
      if(this->GetDestCount()>this->numRecvPosted_){
        for( Int idxRecv = 0; idxRecv < myDests_.size(); ++idxRecv ){
          Int iProc = myDests_[idxRecv];
          //assert(msgSize_>=0);
          MPI_Irecv( (char*)remoteData_[idxRecv], msgSize_, MPI_BYTE, 
              iProc, tag_,comm_, &myRequests_[idxRecv] );
          this->numRecvPosted_++;
        } // for (iProc)
      }
    }




    protected:
    void Reduce( Int idxRecv, Int idReq){
      //add thing to my data
      blas::Axpy(msgSize_/sizeof(T), ONE<T>(), remoteData_[idxRecv], 1, myData_, 1 );
    }

    void Forward(){ 
      //forward to my root if I have reseived everything
      Int iProc = myRoot_;
      // Use Isend to send to multiple targets
      if(myData_==NULL){
        MPI_Isend( NULL, 0, MPI_BYTE, 
            iProc, tag_,comm_, &sendRequest_ );
#ifdef COMM_PROFILE
      PROFILE_COMM(myGRank_,Granks_[iProc],tag_,0);
#endif
      }
      else{
        MPI_Isend( (char*)myData_, msgSize_, MPI_BYTE, 
            iProc, tag_,comm_, &sendRequest_ );
#ifdef COMM_PROFILE
      PROFILE_COMM(myGRank_,Granks_[iProc],tag_,msgSize_);
#endif
      }

#if ( _DEBUGlevel_ >= 1 ) || defined(REDUCE_VERBOSE)
        statusOFS<<myRank_<<" FWD to "<<iProc<<" on tag "<<tag_<<std::endl;
#endif

      fwded_ = true;
      
    }

};


template< typename T>
class FTreeReduce: public TreeReduce<T>{
    protected:
    virtual void buildTree(Int * ranks, Int rank_cnt){

      Int idxStart = 0;
      Int idxEnd = rank_cnt;



      this->myRoot_ = ranks[0];

      if(this->myRank_==this->myRoot_){
        this->myDests_.insert(this->myDests_.begin(),&ranks[1],&ranks[0]+rank_cnt);
      }

#if ( _DEBUGlevel_ >= 1 ) || defined(REDUCE_VERBOSE)
      statusOFS<<"My root is "<<this->myRoot_<<std::endl;
      statusOFS<<"My dests are ";
      for(int i =0;i<this->myDests_.size();++i){statusOFS<<this->myDests_[i]<<" ";}
      statusOFS<<std::endl;
#endif
    }

    void Reduce( ){
      //add thing to my data
      blas::Axpy(this->msgSize_/sizeof(T), ONE<T>(), this->remoteData_[0], 1, this->myData_, 1 );
    }



    public:
    FTreeReduce(const MPI_Comm & pComm, Int * ranks, Int rank_cnt, Int msgSize):TreeReduce<T>(pComm, ranks, rank_cnt, msgSize){
      buildTree(ranks,rank_cnt);
    }

    virtual void PostFirstRecv()
    {
//        if(!this->isAllocated_){
//          this->AllocRecvBuffers();
//        }
        if(this->isAllocated_ && this->GetDestCount()>this->numRecvPosted_){
          MPI_Irecv( (char*)this->remoteData_[0], this->msgSize_, MPI_BYTE, 
            MPI_ANY_SOURCE, this->tag_,this->comm_, &this->myRequests_[0] );
          this->numRecvPosted_++;
        }
    }

    virtual void AllocRecvBuffers(){
      if(this->GetDestCount()>0){
        this->remoteData_.Resize(1);

        this->myRecvBuffers_.Resize(this->msgSize_);

        this->remoteData_[0] = (T*)&(this->myRecvBuffers_[0]);

        this->myRequests_.Resize(1);
        SetValue(this->myRequests_,MPI_REQUEST_NULL);
        this->myStatuses_.Resize(1);
        this->recvIdx_.Resize(1);
      }

      this->sendRequest_ = MPI_REQUEST_NULL;

      this->isAllocated_ = true;
    }

    virtual bool Progress(){

      if(!this->isAllocated_){
        return true;
      }


      if(this->myRank_==this->myRoot_ && this->isAllocated_){
        this->isReady_=true;
      }
     
//      if(this->numRecvPosted_==0){ 
//        this->PostFirstRecv();
//      }

      bool retVal = this->AccumulationDone();
      if(this->isReady_ && !retVal){

        //assert(this->isAllocated_);

        //mpi_test_some on my requests
        int recvCount = -1;
        int reqCnt = 1;

        MPI_Testsome(reqCnt,&this->myRequests_[0],&recvCount,&this->recvIdx_[0],&this->myStatuses_[0]);
        //MPI_Waitsome(reqCnt,&myRequests_[0],&recvCount,&recvIdx_[0],&myStatuses_[0]);
        //if something has been received, accumulate and potentially forward it
        for(Int i = 0;i<recvCount;++i ){
          Int idx = this->recvIdx_[i];

         if(idx!=MPI_UNDEFINED){

          Int size = 0;
          MPI_Get_count(&this->myStatuses_[i], MPI_BYTE, &size);


#if ( _DEBUGlevel_ >= 1 ) || defined(REDUCE_VERBOSE)

        statusOFS<<this->myRank_<<" RECVD from "<<this->myStatuses_[i].MPI_SOURCE<<" on tag "<<this->tag_<<std::endl;
#endif
          if(size>0){
            //If myData is 0, allocate to the size of what has been received
            if(this->myData_==NULL){
              //assert(size==this->msgSize_);
              this->myLocalBuffer_.Resize(this->msgSize_);
    
              this->myData_ = (T*)&this->myLocalBuffer_[0];
              Int nelem = this->msgSize_/sizeof(T);
              std::fill(this->myData_,this->myData_+nelem,ZERO<T>());
            }

            this->Reduce();
          }

          this->numRecv_++;
          }
        }

        if(recvCount>0){
          this->PostFirstRecv();
        }
      }
      else if (this->isReady_ && this->sendRequest_ == MPI_REQUEST_NULL && this->myRoot_ != this->myRank_ && !this->fwded_){
        //free the unnecessary arrays
        this->myRecvBuffers_.Clear();
        this->myRequests_.Clear();
        this->myStatuses_.Clear();
        this->recvIdx_.Clear();

        //Forward
        this->Forward();
        retVal = false;
      }
      else{
        retVal = this->IsDone();
        if(retVal){
        //free the unnecessary arrays
        this->myRecvBuffers_.Clear();
        this->myRequests_.Clear();
        this->myStatuses_.Clear();
        this->recvIdx_.Clear();
        }
      }

      return retVal;
    }


    virtual FTreeReduce * clone() const{
      FTreeReduce * out = new FTreeReduce(*this);
      return out;
    }



};



template< typename T>
class BTreeReduce: public TreeReduce<T>{
    protected:
    virtual void buildTree(Int * ranks, Int rank_cnt){
      Int idxStart = 0;
      Int idxEnd = rank_cnt;



      Int prevRoot = ranks[0];
      while(idxStart<idxEnd){
        Int curRoot = ranks[idxStart];
        Int listSize = idxEnd - idxStart;

        if(listSize == 1){
          if(curRoot == this->myRank_){
            this->myRoot_ = prevRoot;
            break;
          }
        }
        else{
          Int halfList = floor(ceil(double(listSize) / 2.0));
          Int idxStartL = idxStart+1;
          Int idxStartH = idxStart+halfList;

          if(curRoot == this->myRank_){
            if ((idxEnd - idxStartH) > 0 && (idxStartH - idxStartL)>0){
              Int childL = ranks[idxStartL];
              Int childR = ranks[idxStartH];

              this->myDests_.push_back(childL);
              this->myDests_.push_back(childR);
            }
            else if ((idxEnd - idxStartH) > 0){
              Int childR = ranks[idxStartH];
              this->myDests_.push_back(childR);
            }
            else{
              Int childL = ranks[idxStartL];
              this->myDests_.push_back(childL);
            }
            this->myRoot_ = prevRoot;
            break;
          } 

          if( this->myRank_ < ranks[idxStartH]){
            idxStart = idxStartL;
            idxEnd = idxStartH;
          }
          else{
            idxStart = idxStartH;
          }
          prevRoot = curRoot;
        }

      }

#if ( _DEBUGlevel_ >= 1 ) || defined(REDUCE_VERBOSE)
      statusOFS<<"My root is "<<this->myRoot_<<std::endl;
      statusOFS<<"My dests are ";
      for(int i =0;i<this->myDests_.size();++i){statusOFS<<this->myDests_[i]<<" ";}
      statusOFS<<std::endl;
#endif
    }
    public:
    BTreeReduce(const MPI_Comm & pComm, Int * ranks, Int rank_cnt, Int msgSize):TreeReduce<T>(pComm, ranks, rank_cnt, msgSize){
      buildTree(ranks,rank_cnt);
    }

    virtual BTreeReduce * clone() const{
      BTreeReduce * out = new BTreeReduce(*this);
      return out;
    }
};


template< typename T>
class ModBTreeReduce: public TreeReduce<T>{
    protected:
    double rseed_;
    virtual void buildTree(Int * ranks, Int rank_cnt){

      Int idxStart = 0;
      Int idxEnd = rank_cnt;

      //sort the ranks with the modulo like operation
      if(rank_cnt>1){
        //generate a random position in [1 .. rand_cnt]
        //Int new_idx = (int)((rand()+1.0) * (double)rank_cnt / ((double)RAND_MAX+1.0));
        //srand(ranks[0]+rank_cnt);
        //Int new_idx = rseed_%(rank_cnt-1)+1;
        Int new_idx = (int)((rank_cnt - 0) * ( (double)this->rseed_ / (double)RAND_MAX ) + 0);// (this->rseed_)%(rank_cnt-1)+1;


        Int * new_start = &ranks[new_idx];
//        for(int i =0;i<rank_cnt;++i){statusOFS<<ranks[i]<<" ";} statusOFS<<std::endl;
        
//        Int * new_start = std::lower_bound(&ranks[1],&ranks[0]+rank_cnt,ranks[0]);
        //just swap the two chunks   r[0] | r[1] --- r[new_start-1] | r[new_start] --- r[end]
        // becomes                   r[0] | r[new_start] --- r[end] | r[1] --- r[new_start-1] 
        std::rotate(&ranks[1], new_start, &ranks[0]+rank_cnt);
//        for(int i =0;i<rank_cnt;++i){statusOFS<<ranks[i]<<" ";} statusOFS<<std::endl;
      }

      Int prevRoot = ranks[0];
      while(idxStart<idxEnd){
        Int curRoot = ranks[idxStart];
        Int listSize = idxEnd - idxStart;

        if(listSize == 1){
          if(curRoot == this->myRank_){
            this->myRoot_ = prevRoot;
            break;
          }
        }
        else{
          Int halfList = floor(ceil(double(listSize) / 2.0));
          Int idxStartL = idxStart+1;
          Int idxStartH = idxStart+halfList;

          if(curRoot == this->myRank_){
            if ((idxEnd - idxStartH) > 0 && (idxStartH - idxStartL)>0){
              Int childL = ranks[idxStartL];
              Int childR = ranks[idxStartH];

              this->myDests_.push_back(childL);
              this->myDests_.push_back(childR);
            }
            else if ((idxEnd - idxStartH) > 0){
              Int childR = ranks[idxStartH];
              this->myDests_.push_back(childR);
            }
            else{
              Int childL = ranks[idxStartL];
              this->myDests_.push_back(childL);
            }
            this->myRoot_ = prevRoot;
            break;
          } 

          //not true anymore ?
          //first half to 
TIMER_START(FIND_RANK);
          Int * pos = std::find(&ranks[idxStartL], &ranks[idxStartH], this->myRank_);
TIMER_STOP(FIND_RANK);
          if( pos != &ranks[idxStartH]){
            idxStart = idxStartL;
            idxEnd = idxStartH;
          }
          else{
            idxStart = idxStartH;
          }
          prevRoot = curRoot;
        }

      }

#if ( _DEBUGlevel_ >= 1 ) || defined(REDUCE_VERBOSE)
      statusOFS<<"My root is "<<this->myRoot_<<std::endl;
      statusOFS<<"My dests are ";
      for(int i =0;i<this->myDests_.size();++i){statusOFS<<this->myDests_[i]<<" ";}
      statusOFS<<std::endl;
#endif
    }
    public:
    ModBTreeReduce(const MPI_Comm & pComm, Int * ranks, Int rank_cnt, Int msgSize, double rseed):TreeReduce<T>(pComm, ranks, rank_cnt, msgSize){
      this->rseed_ = rseed;
      buildTree(ranks,rank_cnt);
    }

    virtual void Copy(const ModBTreeReduce & Tree){
      this->comm_ = Tree.comm_;
      this->myRank_ = Tree.myRank_;
      this->myRoot_ = Tree.myRoot_; 
      this->msgSize_ = Tree.msgSize_;

      this->numRecv_ = Tree.numRecv_;
      this->tag_= Tree.tag_;
      this->mainRoot_= Tree.mainRoot_;
      this->isReady_ = Tree.isReady_;
      this->myDests_ = Tree.myDests_;


      this->myData_ = Tree.myData_;
      this->sendRequest_ = Tree.sendRequest_;
      this->fwded_= Tree.fwded_;
      this->isAllocated_= Tree.isAllocated_;
      this->numRecvPosted_= Tree.numRecvPosted_;

      this->myLocalBuffer_ = Tree.myLocalBuffer_;
      this->myRecvBuffers_ = Tree.myRecvBuffers_;
      this->remoteData_ = Tree.remoteData_;
      this->myRequests_ = Tree.myRequests_;
      this->myStatuses_ = Tree.myStatuses_;
      this->recvIdx_ = Tree.recvIdx_;
      this->rseed_ = Tree.rseed_;
    }
 



    virtual ModBTreeReduce * clone() const{
      ModBTreeReduce * out = new ModBTreeReduce(*this);
      return out;
    }

};



    inline TreeBcast * TreeBcast::Create(const MPI_Comm & pComm, Int * ranks, Int rank_cnt, Int msgSize, double rseed){
      //get communicator size
      Int nprocs = 0;
      MPI_Comm_size(pComm, &nprocs);

//      return new PalmTreeBcast(pComm,ranks,rank_cnt,msgSize);
//      return new ModBTreeBcast(pComm,ranks,rank_cnt,msgSize, rseed);
//      return new RandBTreeBcast(pComm,ranks,rank_cnt,msgSize);

      if(nprocs<=FTREE_LIMIT){
        return new FTreeBcast(pComm,ranks,rank_cnt,msgSize);
      }
      else{
        return new ModBTreeBcast(pComm,ranks,rank_cnt,msgSize, rseed);
        //return new BTreeBcast(pComm,ranks,rank_cnt,msgSize);
      }




    }




template< typename T>
    inline TreeReduce<T> * TreeReduce<T>::Create(const MPI_Comm & pComm, Int * ranks, Int rank_cnt, Int msgSize, double rseed){
      //get communicator size
      Int nprocs = 0;
      MPI_Comm_size(pComm, &nprocs);


      if(nprocs<=FTREE_LIMIT){
#if ( _DEBUGlevel_ >= 1 ) || defined(REDUCE_VERBOSE)
statusOFS<<"FLAT TREE USED"<<endl;
#endif
        return new FTreeReduce<T>(pComm,ranks,rank_cnt,msgSize);
      }
      else{
#if ( _DEBUGlevel_ >= 1 ) || defined(REDUCE_VERBOSE)
statusOFS<<"BINARY TREE USED"<<endl;
#endif
        return new ModBTreeReduce<T>(pComm,ranks,rank_cnt,msgSize, rseed);
        //return new BTreeReduce<T>(pComm,ranks,rank_cnt,msgSize);
      }
    }








}

#endif
