#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <cstdlib>

using namespace std;

#define AVG 1
#define TAG 2

#define TAG_COUNT 10 

int main(int argc, char * argv[]){
  string filename,dummy;
  int sender, receiver, tag, size;

  map< int, double > agg;
  //map< int, double > count;

  int op = 0;
  int nproc = 0;
  int tag_only = -1;
  int tag_count = -1;
  if(argc>2){
    op = atoi(argv[1]);
    nproc = atoi(argv[2]);
  }
  if(argc>3){
    tag_only = atoi(argv[3]);
  }

  vector<int> count(nproc,0);

  while(cin >> filename){
    cerr<<filename<<endl;

    ifstream fin(filename.c_str());
    fin >> dummy >> dummy >> dummy >> dummy;
    while(fin >> sender >> receiver >> tag >> size){

      if(tag_only!=-1){
        tag = tag % TAG_COUNT;
        if( tag != tag_only){ continue;}
      }

      map<int, double>::iterator it=agg.find(sender);
      if(it==agg.end()){
        agg.insert( std::pair<int,double>(sender,double(size)) );
      }
      else{
        it->second += double(size);
      }
      if(op==AVG){
        count[sender-1]++;
      }


    }
    fin.close();


  }

  for(map<int, double>::iterator it = agg.begin(); it != agg.end(); it++){
    cout << it->first << " " << ((op==AVG)?it->second/(double)count[it->first-1]:it->second)<<endl;
  }

  return 0;
}
