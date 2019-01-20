#ifndef _PEXSI_TREE_HPP_
#define _PEXSI_TREE_HPP_

#include "pexsi/environment.hpp"
#include "pexsi/timer.h"

#include <vector>
#include <map>
#include <algorithm>
#include <string>
//#include <random>

// options to switch from a flat bcast/reduce tree to a binary tree

#ifndef FTREE_LIMIT
#define FTREE_LIMIT 16
#endif



namespace PEXSI{


extern std::map< MPI_Comm , std::vector<int> > commGlobRanks;

#ifdef NEW_BCAST

template< typename T>
  class TreeBcast2{
  protected:

    T * myData_;
    MPI_Request recvRequest_;
    NumVec<char> myRecvBuffer_;

    NumVec<MPI_Request> myRequests_;
    NumVec<MPI_Status> myStatuses_;

    bool done_;
    bool fwded_;
    //    bool isAllocated_;

    Int myRoot_;
    MPI_Comm comm_;
    vector<Int> myDests_;
    Int myRank_;
    Int msgSize_;
    bool isReady_;
    Int mainRoot_;
    Int tag_;
    Int numRecv_;

#ifdef COMM_PROFILE_BCAST
  protected:
    Int myGRoot_;
    Int myGRank_;
    //vector<int> Granks_;
  public:
    void SetGlobalComm(const MPI_Comm & pGComm){
      if(commGlobRanks.count(comm_)==0){
        MPI_Group group2 = MPI_GROUP_NULL;
        MPI_Comm_group(pGComm, &group2);
        MPI_Group group1 = MPI_GROUP_NULL;
        MPI_Comm_group(comm_, &group1);

        int size;
        MPI_Comm_size(comm_,&size);
        vector<int> globRanks(size);
        vector<int> Lranks(size);
        for(int i = 0; i<size;++i){Lranks[i]=i;}
        MPI_Group_translate_ranks(group1, size, &Lranks[0],group2, &globRanks[0]);
        commGlobRanks[comm_] = globRanks;
      }
      myGRoot_ = commGlobRanks[comm_][myRoot_];
      myGRank_ = commGlobRanks[comm_][myRank_];
      //Granks_.resize(myDests_.size());
      //for(int i = 0; i<myDests_.size();++i){
      //  Granks_[i] = globRanks[myDests_[i]];
      //}

      //statusOFS<<myDests_<<std::endl;
      //statusOFS<<Granks_<<std::endl;
    }
#endif



  protected:
    virtual void buildTree(Int * ranks, Int rank_cnt)=0;





  public:

    TreeBcast2(){
      comm_ = MPI_COMM_WORLD;
      myRank_=-1;
      myRoot_ = -1;
      msgSize_ = -1;
      numRecv_ = -1;
      tag_=-1;
      mainRoot_=-1;
      isReady_ = false;
      myData_ = NULL;
      recvRequest_ = MPI_REQUEST_NULL;
      fwded_=false;
      //      isAllocated_=false;
      done_ = false;
    }


    TreeBcast2(const MPI_Comm & pComm, Int * ranks, Int rank_cnt,Int msgSize):TreeBcast2(){
      comm_ = pComm;
      int tmprank;
      MPI_Comm_rank(comm_,&tmprank);
      myRank_ = tmprank;
      myRoot_ = -1;
      msgSize_ = msgSize;
      numRecv_ = 0;
      tag_=-1;
      mainRoot_=ranks[0];
      isReady_ = false;
    }


    virtual TreeBcast2 * clone() const = 0;

    TreeBcast2(const TreeBcast2 & Tree){
      this->Copy(Tree);
    }

    virtual void Copy(const TreeBcast2 & Tree){
      this->comm_ = Tree.comm_;
      this->myRank_ = Tree.myRank_;
      this->myRoot_ = Tree.myRoot_;
      this->msgSize_ = Tree.msgSize_;

      this->numRecv_ = Tree.numRecv_;
      this->tag_= Tree.tag_;
      this->mainRoot_= Tree.mainRoot_;
      this->isReady_ = Tree.isReady_;
      this->myDests_ = Tree.myDests_;


      this->recvRequest_ = Tree.recvRequest_;
      this->myRecvBuffer_ = Tree.myRecvBuffer_;
      this->myRequests_ = Tree.myRequests_;
      this->myStatuses_ = Tree.myStatuses_;
      this->myData_ = Tree.myData_;
      if(Tree.myData_==(T*)&Tree.myRecvBuffer_[0]){
        this->myData_=(T*)&this->myRecvBuffer_[0];
      }




      this->fwded_= Tree.fwded_;
      //      this->isAllocated_= Tree.isAllocated_;
      this->done_= Tree.done_;
    }

    void Reset(){
      assert(done_);
      CleanupBuffers();
      done_=false;
      myData_ = NULL;
      recvRequest_ = MPI_REQUEST_NULL;
      fwded_=false;
      //      isAllocated_=false;
      isReady_=false;
      numRecv_ = 0;
    }

    //bool IsAllocated(){return isAllocated_;}

    virtual ~TreeBcast2(){
      CleanupBuffers();
    }


    static TreeBcast2<T> * Create(const MPI_Comm & pComm, Int * ranks, Int rank_cnt, Int msgSize,double rseed);

    virtual inline Int GetNumRecvMsg(){return numRecv_;}
    virtual inline Int GetNumMsgToSend(){return GetDestCount();}
    inline void SetDataReady(bool rdy){
      isReady_=rdy;
      //numRecv_ = rdy?1:0;
    }
    inline void SetTag(Int tag){ tag_ = tag;}
    inline int GetTag(){ return tag_;}

    bool IsDone(){return done_;}
    bool IsDataReady(){return isReady_;}
    bool IsDataReceived(){return numRecv_==1;}

    Int * GetDests(){ return &myDests_[0];}
    Int GetDest(Int i){ return myDests_[i];}
    Int GetDestCount(){ return myDests_.size();}
    Int GetRoot(){ return myRoot_;}

    bool IsRoot(){ return myRoot_==myRank_;}
    Int GetMsgSize(){ return msgSize_;}

    void ForwardMessage( ){
      if(myRequests_.m()!=GetDestCount()){
        myRequests_.Resize(GetDestCount());
        SetValue(myRequests_,MPI_REQUEST_NULL);
      }
      for( Int idxRecv = 0; idxRecv < myDests_.size(); ++idxRecv ){
        Int iProc = myDests_[idxRecv];
        // Use Isend to send to multiple targets
        MPI_Isend( myData_, msgSize_, MPI_BYTE,
            iProc, tag_,comm_, &myRequests_[idxRecv] );

#if ( _DEBUGlevel_ >= 1 ) || defined(BCAST_VERBOSE)
        statusOFS<<myRank_<<" FWD to "<<iProc<<" on tag "<<tag_<<std::endl;
#endif
#ifdef COMM_PROFILE_BCAST
        //        statusOFS<<idxRecv<<std::endl;
        //        statusOFS<<myDests_<<std::endl;
        //        statusOFS<<Granks_<<std::endl;
        //PROFILE_COMM(myGRank_,Granks_[idxRecv],tag_,msgSize_);
        PROFILE_COMM(myGRank_,commGlobRanks[comm_][iProc],tag_,msgSize_);
#endif
      } // for (iProc)
      fwded_ = true;
    }

    void CleanupBuffers(){
      myRequests_.Clear();
      myStatuses_.Clear();
      myRecvBuffer_.Clear();
    }


    void SetLocalBuffer(T * locBuffer){
      if(myData_!=NULL && myData_!=locBuffer){
        if(numRecv_>0){
          CopyLocalBuffer(locBuffer);
        }
        if(!fwded_){
          myRecvBuffer_.Clear();
        }
      }

      myData_ = locBuffer;
    }

    //async wait and forward
    virtual bool Progress(){
      if(done_){
        return true;
      }

      bool done = false;

      if (myRank_==myRoot_){
        if(isReady_){
          if(!fwded_){
#if ( _DEBUGlevel_ >= 1 ) || defined(BCAST_VERBOSE)
            statusOFS<<myRank_<<" FORWARDING on tag "<<tag_<<std::endl;
#endif
            ForwardMessage();
          }
          else{

            if(myStatuses_.m()!=GetDestCount()){
              myStatuses_.Resize(GetDestCount());
              recvRequest_ = MPI_REQUEST_NULL;
            }
            //test the send requests
            int flag = 0;
            int reqCnt = GetDestCount();
            if(reqCnt>0){
              assert(reqCnt == myRequests_.m());
              MPI_Testall(reqCnt,&myRequests_[0],&flag,&myStatuses_[0]);
              done = flag==1;
            }
            else{
              done=true;
            }
          }
        }
      }
      else{
        bool received = (numRecv_==1);

        if(!received){
          if(recvRequest_ == MPI_REQUEST_NULL ){
#if ( _DEBUGlevel_ >= 1 ) || defined(BCAST_VERBOSE)
            statusOFS<<myRank_<<" POSTING RECV on tag "<<tag_<<std::endl;
#endif
            //post the recv
            PostRecv();
          }
          else{

            if(myStatuses_.m()!=GetDestCount()){
              myStatuses_.Resize(GetDestCount());
              recvRequest_ = MPI_REQUEST_NULL;
            }
#if ( _DEBUGlevel_ >= 1 ) || defined(BCAST_VERBOSE)
            statusOFS<<myRank_<<" TESTING RECV on tag "<<tag_<<std::endl;
#endif
            //test
            int flag = 0;
            MPI_Status stat;
            int test = MPI_Test(&recvRequest_,&flag,&stat);
            assert(test==MPI_SUCCESS);

            if(flag==1){
              numRecv_=1;
              received = true;

              if(!fwded_){
#if ( _DEBUGlevel_ >= 1 ) || defined(BCAST_VERBOSE)
                statusOFS<<myRank_<<" FORWARDING on tag "<<tag_<<std::endl;
#endif
                ForwardMessage();
              }
            }
          }
        }
        else {
          assert(fwded_);
          //test the send requests
          int flag = 0;
          int reqCnt = GetDestCount();
          if(reqCnt>0){
            assert(reqCnt == myRequests_.m());
            MPI_Testall(reqCnt,&myRequests_[0],&flag,&myStatuses_[0]);
            done = flag==1;
          }
          else{
            done=true;
          }
        }
      }

      if(done){
        //free the unnecessary arrays
        myRequests_.Clear();
        myStatuses_.Clear();
#if ( _DEBUGlevel_ >= 1 ) || defined(BCAST_VERBOSE)
        statusOFS<<myRank_<<" EVERYTHING COMPLETED on tag "<<tag_<<std::endl;
#endif
      }

      done_ = done;

      return done;
    }

    //blocking wait
    void Wait(){
      if(!done_){
        while(!Progress());
      }
    }

    T * GetLocalBuffer(){
      return myData_;
    }

    virtual void PostRecv()
    {
      if(this->numRecv_<1 && this->recvRequest_==MPI_REQUEST_NULL && myRank_!=myRoot_){

        if(myData_==NULL){
          myRecvBuffer_.Resize(msgSize_);
          myData_ = (T*)&myRecvBuffer_[0];
        }
        MPI_Irecv( (char*)this->myData_, this->msgSize_, MPI_BYTE,
            this->myRoot_, this->tag_,this->comm_, &this->recvRequest_ );
      }
    }



    void CopyLocalBuffer(T* destBuffer){
      std::copy((char*)myData_,(char*)myData_+GetMsgSize(),(char*)destBuffer);
    }




  };


template< typename T>
  class FTreeBcast2: public TreeBcast2<T>{
  protected:
    virtual void buildTree(Int * ranks, Int rank_cnt){

      Int idxStart = 0;
      Int idxEnd = rank_cnt;



      this->myRoot_ = ranks[0];

      if(this->myRank_==this->myRoot_){
        this->myDests_.insert(this->myDests_.begin(),&ranks[1],&ranks[0]+rank_cnt);
      }

#if (defined(BCAST_VERBOSE))
      statusOFS<<"My root is "<<this->myRoot_<<std::endl;
      statusOFS<<"My dests are ";
      for(int i =0;i<this->myDests_.size();++i){statusOFS<<this->myDests_[i]<<" ";}
      statusOFS<<std::endl;
#endif
    }



  public:
    FTreeBcast2(const MPI_Comm & pComm, Int * ranks, Int rank_cnt, Int msgSize):TreeBcast2<T>(pComm,ranks,rank_cnt,msgSize){
      //build the binary tree;
      buildTree(ranks,rank_cnt);
    }


    virtual FTreeBcast2 * clone() const{
      FTreeBcast2 * out = new FTreeBcast2(*this);
      return out;
    }
  };

template< typename T>
  class BTreeBcast2: public TreeBcast2<T>{
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

#if (defined(BCAST_VERBOSE))
      statusOFS<<"My root is "<<myRoot_<<std::endl;
      statusOFS<<"My dests are ";
      for(int i =0;i<this->myDests_.size();++i){statusOFS<<this->myDests_[i]<<" ";}
      statusOFS<<std::endl;
#endif
    }



  public:
    BTreeBcast2(const MPI_Comm & pComm, Int * ranks, Int rank_cnt, Int msgSize):TreeBcast2<T>(pComm,ranks,rank_cnt,msgSize){
      //build the binary tree;
      buildTree(ranks,rank_cnt);
    }

    virtual BTreeBcast2<T> * clone() const{
      BTreeBcast2<T> * out = new BTreeBcast2<T>(*this);
      return out;
    }

  };



template< typename T>
  class ModBTreeBcast2: public TreeBcast2<T>{
  protected:
    double rseed_;

    virtual void buildTree(Int * ranks, Int rank_cnt){

      Int idxStart = 0;
      Int idxEnd = rank_cnt;

      //sort the ranks with the modulo like operation
      if(rank_cnt>1){
        //Int new_idx = (int)((rand()+1.0) * (double)rank_cnt / ((double)RAND_MAX+1.0));

        //        srand(ranks[0]+rank_cnt);
        Int new_idx = (Int)rseed_ % (rank_cnt - 1) + 1;
        //Int new_idx = (int)((rank_cnt - 0) * ( (double)this->rseed_ / (double)RAND_MAX ) + 0);// (this->rseed_)%(rank_cnt-1)+1;
        //statusOFS<<"NEW IDX: "<<new_idx<<endl;



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

#if (defined(REDUCE_VERBOSE))
      statusOFS<<"My root is "<<myRoot_<<std::endl;
      statusOFS<<"My dests are ";
      for(int i =0;i<this->myDests_.size();++i){statusOFS<<this->myDests_[i]<<" ";}
      statusOFS<<std::endl;
#endif
    }



  public:
    ModBTreeBcast2(const MPI_Comm & pComm, Int * ranks, Int rank_cnt, Int msgSize, double rseed):TreeBcast2<T>(pComm,ranks,rank_cnt,msgSize){
      //build the binary tree;
      rseed_ = rseed;
      buildTree(ranks,rank_cnt);
    }

    //virtual void Copy(const ModBTreeBcast & Tree){
    //  comm_ = Tree.comm_;
    //  myRank_ = Tree.myRank_;
    //  myRoot_ = Tree.myRoot_;
    //  msgSize_ = Tree.msgSize_;

    //  numRecv_ = Tree.numRecv_;
    //  tag_= Tree.tag_;
    //  mainRoot_= Tree.mainRoot_;
    //  isReady_ = Tree.isReady_;
    //  myDests_ = Tree.myDests_;

    //  rseed_ = Tree.rseed_;
    //  myRank_ = Tree.myRank_;
    //  myRoot_ = Tree.myRoot_;
    //  msgSize_ = Tree.msgSize_;

    //  numRecv_ = Tree.numRecv_;
    //  tag_= Tree.tag_;
    //  mainRoot_= Tree.mainRoot_;
    //  isReady_ = Tree.isReady_;
    //  myDests_ = Tree.myDests_;
    //}

    virtual ModBTreeBcast2 * clone() const{
      ModBTreeBcast2 * out = new ModBTreeBcast2(*this);
      return out;
    }

  };



#endif

































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


#if defined(COMM_PROFILE_BCAST) || defined(COMM_PROFILE)
protected:
  Int myGRank_;
  Int myGRoot_;
  //vector<int> Granks_;
public:
  void SetGlobalComm(const MPI_Comm & pGComm){
    if(commGlobRanks.count(comm_)==0){
      MPI_Group group2 = MPI_GROUP_NULL;
      MPI_Comm_group(pGComm, &group2);
      MPI_Group group1 = MPI_GROUP_NULL;
      MPI_Comm_group(comm_, &group1);

      int size;
      MPI_Comm_size(comm_,&size);
      vector<int> globRanks(size);
      vector<int> Lranks(size);
      for(int i = 0; i<size;++i){Lranks[i]=i;}
      MPI_Group_translate_ranks(group1, size, &Lranks[0],group2, &globRanks[0]);
      commGlobRanks[comm_] = globRanks;
    }
    myGRoot_ = commGlobRanks[comm_][myRoot_];
    myGRank_ = commGlobRanks[comm_][myRank_];
    //   Granks_.resize(myDests_.size());
    //   for(int i = 0; i<myDests_.size();++i){
    //     Granks_[i] = globRanks[myDests_[i]];
    //   }
    //statusOFS<<myDests_<<std::endl;
    //statusOFS<<Granks_<<std::endl;
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
    int tmprank;
    MPI_Comm_rank(comm_,&tmprank);
    myRank_ = tmprank;
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

    tag_= Tree.tag_;
    mainRoot_= Tree.mainRoot_;
    myDests_ = Tree.myDests_;

    //numRecv_ = Tree.numRecv_;
    //isReady_ = Tree.isReady_;
    isReady_ = false;
    numRecv_ = 0;
  }

  virtual TreeBcast * clone() const = 0;

  void Reset(){
    //statusOFS<<"RESET CALLED"<<std::endl;
    this->numRecv_ = 0;
    this->isReady_ = false;
  }



  static TreeBcast * Create(const MPI_Comm & pComm, Int * ranks, Int rank_cnt, Int msgSize,double rseed);

  virtual inline Int GetNumRecvMsg(){return numRecv_;}
  virtual inline Int GetNumMsgToRecv(){return 1;}
  inline void SetDataReady(bool rdy){ isReady_=rdy; }
  inline void SetTag(Int tag){ tag_ = tag;}
  inline int GetTag(){ return tag_;}


  Int * GetDests(){ return &myDests_[0];}
  Int GetDest(Int i){ return myDests_[i];}
  Int GetDestCount(){ return myDests_.size();}
  Int GetRoot(){ return myRoot_;}
  Int GetMsgSize(){ return msgSize_;}

  void ForwardMessage( char * data, size_t size, int tag, MPI_Request * requests ){
    tag_ = tag;
    for( Int idxRecv = 0; idxRecv < myDests_.size(); ++idxRecv ){
      Int iProc = myDests_[idxRecv];
      // Use Isend to send to multiple targets
      MPI_Isend( data, size, MPI_BYTE,
          iProc, tag,comm_, &requests[2*iProc+1] );

#if defined(COMM_PROFILE_BCAST) || defined(COMM_PROFILE)
      //        statusOFS<<idxRecv<<std::endl;
      //        statusOFS<<Granks_<<std::endl;
      //PROFILE_COMM(myGRank_,Granks_[idxRecv],tag,msgSize_);
      PROFILE_COMM(myGRank_,commGlobRanks[comm_][iProc],tag_,msgSize_);
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

#if (defined(BCAST_VERBOSE))
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

#if (defined(BCAST_VERBOSE))
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
      //Int new_idx = (Int)rseed_ % (rank_cnt - 1) + 1;
//      Int new_idx = (int)((rank_cnt - 0) * ( (double)this->rseed_ / (double)RAND_MAX ) + 0);// (this->rseed_)%(rank_cnt-1)+1;
      Int new_idx = (int)(this->rseed_)%(rank_cnt-1)+1;
      //statusOFS<<"NEW IDX: "<<new_idx<<endl;

      assert(new_idx<rank_cnt && new_idx>=1);

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
      assert(idxStart<rank_cnt && idxStart>=0);
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
      assert(idxStartL<rank_cnt && idxStartL>=0);
      assert(idxStartH<rank_cnt && idxStartH>=0);

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

#if (defined(REDUCE_VERBOSE))
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
    ((TreeBcast*)this)->Copy(*((const TreeBcast*)&Tree));
    ///      comm_ = Tree.comm_;
    ///      myRank_ = Tree.myRank_;
    ///      myRoot_ = Tree.myRoot_;
    ///      msgSize_ = Tree.msgSize_;
    ///
    ///      numRecv_ = Tree.numRecv_;
    ///      tag_= Tree.tag_;
    ///      mainRoot_= Tree.mainRoot_;
    ///      isReady_ = Tree.isReady_;
    ///      myDests_ = Tree.myDests_;
    ///
    ///      myRank_ = Tree.myRank_;
    ///      myRoot_ = Tree.myRoot_;
    ///      msgSize_ = Tree.msgSize_;
    ///
    ///      numRecv_ = Tree.numRecv_;
    ///      tag_= Tree.tag_;
    ///      mainRoot_= Tree.mainRoot_;
    ///      isReady_ = Tree.isReady_;
    ///      myDests_ = Tree.myDests_;


    rseed_ = Tree.rseed_;

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

#if (defined(REDUCE_VERBOSE))
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
      numRoots = std::min( rank_cnt, numRoots + (Int)pow(2.0,level));
      Int numNextRoots = std::min(rank_cnt,numRoots + (Int)pow(2.0,(level+1)));
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

#if (defined(BCAST_VERBOSE))
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
  NumVec<char> myLocalBuffer_;
  NumVec<char> myRecvBuffers_;
  NumVec<T *> remoteData_;
  NumVec<MPI_Request> myRequests_;
  NumVec<MPI_Status> myStatuses_;
  NumVec<int> recvIdx_;

  bool fwded_;
  bool done_;
  bool isAllocated_;
  Int numRecvPosted_;

public:
  TreeReduce(const MPI_Comm & pComm, Int * ranks, Int rank_cnt, Int msgSize):TreeBcast(pComm,ranks,rank_cnt,msgSize){
    myData_ = NULL;
    sendRequest_ = MPI_REQUEST_NULL;
    fwded_=false;
    done_=false;
    isAllocated_=false;
    numRecvPosted_= 0;
  }


  virtual TreeReduce * clone() const = 0;

  TreeReduce(const TreeReduce & Tree){
    this->Copy(Tree);
  }

  virtual void Copy(const TreeReduce & Tree){
    ((TreeBcast*)this)->Copy(*(const TreeBcast*)&Tree);

    //      this->comm_ = Tree.comm_;
    //      this->myRank_ = Tree.myRank_;
    //      this->myRoot_ = Tree.myRoot_;
    //      this->msgSize_ = Tree.msgSize_;
    //      this->numRecv_ = Tree.numRecv_;
    //      this->tag_= Tree.tag_;
    //      this->mainRoot_= Tree.mainRoot_;
    //      this->isReady_ = Tree.isReady_;
    //      this->myDests_ = Tree.myDests_;


    this->myData_ = NULL;
    this->sendRequest_ = MPI_REQUEST_NULL;
    this->fwded_= false;
    this->done_= false;
    this->isAllocated_= Tree.isAllocated_;
    this->numRecvPosted_= 0;

    //this->myLocalBuffer_.resize(Tree.myLocalBuffer_.size());
    //this->remoteData_ = Tree.remoteData_;
    //this->recvIdx_ = Tree.recvIdx_;

    CleanupBuffers();
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


  }

  void Reset(){
    //        assert(done_ || myDests_.size()==0);
    CleanupBuffers();
    done_=false;

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

      //statusOFS<<"DOING SUM"<<std::endl;
      //gdb_lock();
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

    if(done_){
      return true;
    }
    if(!isAllocated_){
      return false;
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


    done_ = retVal;
    return retVal;
  }

  //blocking wait
  void Wait(){
    if(!done_){
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
  virtual void Reduce( Int idxRecv, Int idReq){
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
      PROFILE_COMM(myGRank_,myGRoot_,tag_,0);
#endif
    }
    else{
      MPI_Isend( (char*)myData_, msgSize_, MPI_BYTE,
          iProc, tag_,comm_, &sendRequest_ );
#ifdef COMM_PROFILE
      PROFILE_COMM(myGRank_,myGRoot_,tag_,msgSize_);
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

#if (defined(REDUCE_VERBOSE))
    statusOFS<<"My root is "<<this->myRoot_<<std::endl;
    statusOFS<<"My dests are ";
    for(int i =0;i<this->myDests_.size();++i){statusOFS<<this->myDests_[i]<<" ";}
    statusOFS<<std::endl;
#endif
  }

  virtual void Reduce( ){
    //add thing to my data
    blas::Axpy(this->msgSize_/sizeof(T), ONE<T>(), this->remoteData_[0], 1, this->myData_, 1 );


#if (defined(REDUCE_DEBUG))
    statusOFS << std::endl << /*"["<<snode.Index<<"]*/" Recv contrib"<< std::endl;
    for(int i = 0; i < this->msgSize_/sizeof(T); ++i){
      statusOFS<< this->remoteData_[0][i]<< " ";
      if(i%3==0){statusOFS<<std::endl;}
    }
    statusOFS<<std::endl;

    statusOFS << std::endl << /*"["<<snode.Index<<"]*/" Reduce buffer now is"<< std::endl;
    for(int i = 0; i < this->msgSize_/sizeof(T); ++i){
      statusOFS<< this->myData_[i]<< " ";
      if(i%3==0){statusOFS<<std::endl;}
    }
    statusOFS<<std::endl;
#endif

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

#if (defined(REDUCE_VERBOSE))
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

      //Int new_idx = (int)((rank_cnt - 0) * ( (double)this->rseed_ / (double)RAND_MAX ) + 0);// (this->rseed_)%(rank_cnt-1)+1;
      //Int new_idx = (Int)rseed_ % (rank_cnt - 1) + 1;
//      Int new_idx = (int)((rank_cnt - 0) * ( (double)this->rseed_ / (double)RAND_MAX ) + 0);// (this->rseed_)%(rank_cnt-1)+1;
      Int new_idx = (int)(this->rseed_)%(rank_cnt-1)+1;

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

#if (defined(REDUCE_VERBOSE))
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
    ((TreeReduce<T>*)this)->Copy(*((const TreeReduce<T>*)&Tree));
    //this->comm_ = Tree.comm_;
    //this->myRank_ = Tree.myRank_;
    //this->myRoot_ = Tree.myRoot_;
    //this->msgSize_ = Tree.msgSize_;

    //this->numRecv_ = Tree.numRecv_;
    //this->tag_= Tree.tag_;
    //this->mainRoot_= Tree.mainRoot_;
    //this->isReady_ = Tree.isReady_;
    //this->myDests_ = Tree.myDests_;


    //this->myData_ = Tree.myData_;
    //this->sendRequest_ = Tree.sendRequest_;
    //this->fwded_= Tree.fwded_;
    //this->isAllocated_= Tree.isAllocated_;
    //this->numRecvPosted_= Tree.numRecvPosted_;

    //this->myLocalBuffer_ = Tree.myLocalBuffer_;
    //this->myRecvBuffers_ = Tree.myRecvBuffers_;
    //this->remoteData_ = Tree.remoteData_;
    //this->myRequests_ = Tree.myRequests_;
    //this->myStatuses_ = Tree.myStatuses_;
    //this->recvIdx_ = Tree.recvIdx_;
    this->rseed_ = Tree.rseed_;
  }




  virtual ModBTreeReduce * clone() const{
    ModBTreeReduce * out = new ModBTreeReduce(*this);
    return out;
  }

};


template< typename T>
class PalmTreeReduce: public TreeReduce<T>{
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
          if(r==this->myRank_){
            this->myRoot_ = p;
          }

          if(p==this->myRank_){
            this->myDests_.push_back(r);
          }
        }
      }
    }

#if (defined(BCAST_VERBOSE))
    statusOFS<<"My root is "<<this->myRoot_<<std::endl;
    statusOFS<<"My dests are ";
    for(int i =0;i<this->myDests_.size();++i){statusOFS<<this->myDests_[i]<<" ";}
    statusOFS<<std::endl;
#endif
  }



public:
  PalmTreeReduce(const MPI_Comm & pComm, Int * ranks, Int rank_cnt, Int msgSize):TreeReduce<T>(pComm,ranks,rank_cnt,msgSize){
    //build the binary tree;
    buildTree(ranks,rank_cnt);
  }



  virtual void Copy(const PalmTreeReduce & Tree){
    ((TreeReduce<T>*)this)->Copy(*((const TreeReduce<T>*)&Tree));
    //this->comm_ = Tree.comm_;
    //this->myRank_ = Tree.myRank_;
    //this->myRoot_ = Tree.myRoot_;
    //this->msgSize_ = Tree.msgSize_;

    //this->numRecv_ = Tree.numRecv_;
    //this->tag_= Tree.tag_;
    //this->mainRoot_= Tree.mainRoot_;
    //this->isReady_ = Tree.isReady_;
    //this->myDests_ = Tree.myDests_;


    //this->myData_ = Tree.myData_;
    //this->sendRequest_ = Tree.sendRequest_;
    //this->fwded_= Tree.fwded_;
    //this->isAllocated_= Tree.isAllocated_;
    //this->numRecvPosted_= Tree.numRecvPosted_;

    //this->myLocalBuffer_ = Tree.myLocalBuffer_;
    //this->myRecvBuffers_ = Tree.myRecvBuffers_;
    //this->remoteData_ = Tree.remoteData_;
    //this->myRequests_ = Tree.myRequests_;
    //this->myStatuses_ = Tree.myStatuses_;
    //this->recvIdx_ = Tree.recvIdx_;
    //this->rseed_ = Tree.rseed_;
  }




  virtual PalmTreeReduce * clone() const{
    PalmTreeReduce * out = new PalmTreeReduce(*this);
    return out;
  }

};



inline TreeBcast * TreeBcast::Create(const MPI_Comm & pComm, Int * ranks, Int rank_cnt, Int msgSize, double rseed){
  //get communicator size
  int nprocs = 0;
  MPI_Comm_size(pComm, &nprocs);


#if defined(FTREE)
  return new FTreeBcast(pComm,ranks,rank_cnt,msgSize);
#elif defined(MODBTREE)
  return new ModBTreeBcast(pComm,ranks,rank_cnt,msgSize, rseed);
#elif defined(BTREE)
  return new BTreeBcast(pComm,ranks,rank_cnt,msgSize);
#elif defined(PALMTREE)
  return new PalmTreeBcast(pComm,ranks,rank_cnt,msgSize);
#endif


  //      return new PalmTreeBcast(pComm,ranks,rank_cnt,msgSize);
  //      return new ModBTreeBcast(pComm,ranks,rank_cnt,msgSize, rseed);
  //      return new RandBTreeBcast(pComm,ranks,rank_cnt,msgSize);

  if(nprocs<=FTREE_LIMIT){
    return new FTreeBcast(pComm,ranks,rank_cnt,msgSize);
  }
  else{
    return new ModBTreeBcast(pComm,ranks,rank_cnt,msgSize, rseed);
  }




}




template< typename T>
inline TreeReduce<T> * TreeReduce<T>::Create(const MPI_Comm & pComm, Int * ranks, Int rank_cnt, Int msgSize, double rseed){
  //get communicator size
  int nprocs = 0;
  MPI_Comm_size(pComm, &nprocs);

#if defined(FTREE)
  return new FTreeReduce<T>(pComm,ranks,rank_cnt,msgSize);
#elif defined(MODBTREE)
  return new ModBTreeReduce<T>(pComm,ranks,rank_cnt,msgSize, rseed);
#elif defined(BTREE)
  return new BTreeReduce<T>(pComm,ranks,rank_cnt,msgSize);
#elif defined(PALMTREE)
  return new PalmTreeReduce<T>(pComm,ranks,rank_cnt,msgSize);
#endif


  if(nprocs<=FTREE_LIMIT){
#if ( _DEBUGlevel_ >= 1 ) || defined(REDUCE_VERBOSE)
    statusOFS<<"FLAT TREE USED"<<std::endl;
#endif
    return new FTreeReduce<T>(pComm,ranks,rank_cnt,msgSize);
  }
  else{
#if ( _DEBUGlevel_ >= 1 ) || defined(REDUCE_VERBOSE)
    statusOFS<<"BINARY TREE USED"<<std::endl;
#endif
    return new ModBTreeReduce<T>(pComm,ranks,rank_cnt,msgSize, rseed);
    //return new BTreeReduce<T>(pComm,ranks,rank_cnt,msgSize);
  }
}


#ifdef NEW_BCAST
template< typename T>
inline TreeBcast2<T> * TreeBcast2<T>::Create(const MPI_Comm & pComm, Int * ranks, Int rank_cnt, Int msgSize, double rseed){
  //get communicator size
  int nprocs = 0;
  MPI_Comm_size(pComm, &nprocs);




#if defined(FTREE)
  return new FTreeBcast2<T>(pComm,ranks,rank_cnt,msgSize);
#elif defined(MODBTREE)
  return new ModBTreeBcast2<T>(pComm,ranks,rank_cnt,msgSize,rseed);
#elif defined(BTREE)
  return new BTreeBcast2<T>(pComm,ranks,rank_cnt,msgSize);
#endif


  if(nprocs<=FTREE_LIMIT){
#if ( _DEBUGlevel_ >= 1 ) || defined(REDUCE_VERBOSE)
    statusOFS<<"FLAT TREE USED"<<endl;
#endif

    return new FTreeBcast2<T>(pComm,ranks,rank_cnt,msgSize);

  }
  else{
#if ( _DEBUGlevel_ >= 1 ) || defined(REDUCE_VERBOSE)
    statusOFS<<"BINARY TREE USED"<<endl;
#endif
    return new ModBTreeBcast2<T>(pComm,ranks,rank_cnt,msgSize, rseed);
    //        //return new BTreeReduce<T>(pComm,ranks,rank_cnt,msgSize);
  }
}
#endif






}

#endif
