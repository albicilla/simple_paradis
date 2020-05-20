#include <cassert>
#include <iostream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <random>
#include <sstream>
#include <thread>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <climits>

using namespace std;

#define INF 1e9
#define FOR(i,a,b) for(int i=a;i<b;i++)
#define rep(i,b) FOR(i,0,b)
#define dump(x) cerr<<#x<<"="<<x<<endl
#define ALL(a) (a).begin(),(a).end()
#define EACH(e,v) for(auto& e:v)
#define SORT(v) sort(ALL(v))
#define PERM(v) SORT(v);for(bool c##p=1;c##p;c##p=next_permutation(ALL(v)))
#define printArr(x) for(auto itr:x)cerr<<"="<<itr<<endl

#include <omp.h>

class CppClock{
     public:
     std::chrono::system_clock::time_point start;
     CppClock(){
         start = std::chrono::system_clock::now();
     }
     double get(){
         auto end = std::chrono::system_clock::now();
         double elapsed = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count());
         return elapsed;
     }
 };
using ll = long long;

template< class _Type>
inline void _swap(_Type &a, _Type&b)  {
  _Type temp = b;
  b = a;
  a = temp;
}





const int kisuu=256;
const int MaxThreadNum=224;
const long long MaxDataSize=10000000000;
const long long MaxDataNum=4294967295;
const int MaxKisuu=256;
int threadNum;
int NumRange=0;
ll unsortCnt=0;

std::vector<ll> Dataset;
ll Datasize;



//Sort 対象データセット
void makeDataset(){
    
    Dataset.clear();

    std::uniform_int_distribution<long long> dist(0.0,INT_MAX);

    random_device seed_gen;
    default_random_engine engine(seed_gen());

    Dataset.resize(Datasize);
    
    for(int i=0;i<Datasize;i++){
       int result = dist(engine);
       Dataset[i]=(result);
    }
    
   
}

template<class D>
bool issorted(std::vector<D>&arr){
    for(int i=0;i<arr.size()-1;i++){
        if(arr[i]>arr[i+1]){
            std::cerr<<"i="<<i<<std::endl;
            return false;
        }
    }
    return true;
}




static const int kRadixBits = 8;
static const size_t kInsertSortThreshold = 0;
static const int kRadixMask = (1 << kRadixBits) - 1;
static const int kRadixBin = 1 << kRadixBits;



template<class D>
inline int determineDigitBucket(int stage,D num){
  return ((num>>(8*stage))&kRadixMask);
}


void report_num_threads(int level)
{
    #pragma omp single
    {
        printf("Level %d: number of threads in the team - %d\n",
               level, omp_get_num_threads());
    }
}

template<class T>
bool compare(const T &x,const T &y){
  return x < y;
}
template <class RandomIt>
inline void insert_sort_core_(RandomIt s, RandomIt e)
{
    for (RandomIt i = s + 1; i < e; ++i) {
        if (compare(*i, *(i - 1))) {
            RandomIt j;
            auto tmp = *i;
            *i = *(i - 1);
            for (j = i - 1; j > s && compare(tmp, *(j - 1)); --j) {
                *j = *(j - 1);
            }
            *j = tmp;
        }
    }
}



double createHistTime=0.0;
double paradisParmutateTime=0.0;
double paradisRepairTime=0.0;
double createGhGt=0.0;
double quickSortTime = 0.0;
double MiTime = 0.0;
int localhistCall=0;
int recum=0;
int CallProcesses=0;
int swapNum=0;
int needRepairNum=0;
int prefer_insert=0;
template<class D,int kth_byte>
inline void RadixSort(vector<D>& arr,int elenum,int start,int processes=1){
    int cnt[MaxKisuu];
     
    rep(i,kisuu){cnt[i]=0;}    

    //step1
    int part=elenum/processes;
    int res=elenum%(processes);

    //threadNum個のヒストグラムを作る 
    int **localHists;
    localHists = (int**)malloc(sizeof(int *)*processes);
    for(int i=0;i<processes;i++){
      localHists[i]=(int*)malloc(sizeof(int)*kisuu);
    }
    localhistCall+=elenum;

    int gh[MaxKisuu],gt[MaxKisuu],starts[MaxKisuu],ends[MaxKisuu];
    int ph[MaxThreadNum][MaxKisuu];
    int pt[MaxThreadNum][MaxKisuu];
    
    //indicesをpartitionしてカウント
    int nvalue = elenum;

    // sum i(Ci) > 0 初めは全要素
    int SumCi=elenum;
    int roop=0;
    //for paradis repair
    int pfp[processes+1];
    int var_p=processes;

#pragma omp parallel num_threads(processes)   
      {

          int th = omp_get_thread_num();
	  #pragma omp for
          rep(i,kisuu){
	    rep(th,processes)localHists[th][i]=0;
	  }
	  #pragma omp barrier
	  #pragma omp for
	  for(int i=start;i<start+elenum;i++){
	    assert(arr[i]>=(D)0);
	    D digit=determineDigitBucket(kth_byte,arr[i]);
	    localHists[th][digit]++;
	  }
	  #pragma omp barrier
          #pragma omp for      
          for(int i=0;i<kisuu;i++){
            for(int j=0;j<processes;j++){
             cnt[i]+=localHists[j][i];
            }
          }
    #pragma omp barrier
    #pragma omp for
    for(int i=0;i<processes;i++){
	free(localHists[i]);
    }
    #pragma omp single
    free(localHists);


    gh[0]=start;
    gt[0]=start+cnt[0];
    starts[0]=gh[0];
    //step2
    #pragma omp single
    for(int i=1;i<kRadixBin;i++){
        //calc ghi
        gh[i]=gh[i-1]+cnt[i-1];
        //calc gti;
        gt[i]=gh[i]+cnt[i];

        starts[i]=gh[i];
    }

    #pragma omp barrier
    //step3
    while(SumCi!=0){
      

#pragma omp for
      for(int ii=0;ii<var_p;ii++)
        {
	    int pID=omp_get_thread_num();

	
	    for(int i=0;i<kisuu;i++){
	      int part=(int)(gt[i]-gh[i])/var_p;
	      int res=(int)(gt[i]-gh[i])%(var_p);
	      
	      if(pID<var_p-1){
		ph[pID][i]=part*pID+gh[i];
		pt[pID][i]=part*(pID+1)+gh[i];
	      }else{
		ph[pID][i]=part*pID+gh[i];
		pt[pID][i]=part*(pID+1)+gh[i]+res;
	      }
	    }

	    for(int i=0;i<kisuu;i++){
	      int head=ph[pID][i];
	      while(head<pt[pID][i]){
		D v=arr[head];
		int k=determineDigitBucket(kth_byte,v);
		while(k!=i&&ph[pID][k]<pt[pID][k]){
		  _swap(v,arr[ph[pID][k]++]);
		  k=determineDigitBucket(kth_byte,v);
		  //swapNum++;
		}
		if(k==i){
		  arr[head++]=arr[ph[pID][i]];
		  arr[ph[pID][i]++]=v;
		}else{
		  arr[head++]=v;
		}
	      }
	    }	    
	}//end of omp pfp
#pragma omp single
      {
	int pfpN=kisuu/var_p;
	int pfpM=kisuu%var_p;
	pfp[0]=0;
	int pfpMR=0;
	for(int i=1;i<var_p+1;i++){
	  if(pfpMR<pfpM)pfpMR++;
	  pfp[i]=i*pfpN+pfpMR;
	}
      }
#pragma omp single
      SumCi=0;      
#pragma omp barrier
#pragma omp for
	for(int k=0;k<var_p;k++){
		for(int i=pfp[k];i<pfp[k+1];i++){
		  int tail=gt[i];
		  {
		    for(int pID=0;pID<var_p;pID++){
		      int head=ph[pID][i];
		      while(head<pt[pID][i]&&head<tail){
			D v=arr[head++];
			if(determineDigitBucket(kth_byte,v)!=i){
			  //fix < to <=
			  while(head<=tail){
                  D w=arr[--tail];
                  if(determineDigitBucket(kth_byte,w)==i){
			      arr[head-1]=w;
			      arr[tail]=v;
			      break;
			    }
			  }

			}
		      }
		    }
		  }
		  gh[i]=tail;
		}
	}
#pragma omp barrier
#pragma omp single
	{
	int prevSumCi=SumCi;
	SumCi=0;
	for(int i=0;i<kisuu;i++){
	  SumCi+=(gt[i]-gh[i]);
	}
	}

 #pragma omp barrier
      
    }//end of while
    }//end of omp2

     
    if(kth_byte>0){
      {
#pragma omp parallel  num_threads(processes) 
#pragma omp single
	{

	for(int i=0;i<kisuu;i++){
	  int nextStageThreads=1;
	  nextStageThreads=processes*(cnt[i]*(log(cnt[i])/log(kRadixBin))/(elenum*(log(elenum)/log(kRadixBin))));
      if(cnt[i]>64){
	  #pragma omp task	   
	    RadixSort<D,(kth_byte > 0 ? (kth_byte - 1) : 0)>(arr,cnt[i],starts[i],max((int)nextStageThreads,(int)1));
	  }
	  else if(cnt[i]>1){
	    //if elements less than 64 call insert sort
	    insert_sort_core_(arr.begin()+starts[i],arr.begin()+starts[i]+cnt[i]);
	  }
	}
	#pragma omp taskwait
	}
      }
    }
}

signed main(int argc, char** argv){

    if (argc < 2) {  // ./a.out 4 って書くと4スレッドになる
        printf("ERROR! set thread num\n");
        return 1;  // 指定漏れがあった場合は自殺
    }
    

    threadNum = atoi(argv[1]);

    if(threadNum>MaxThreadNum){
        printf("ERROR! too many thread num\n");
        return 1;  // 指定漏れがあった場合は自殺
    }
    if (argc < 3) {  
        printf("ERROR! set datasize\n");
        return 1;  // 指定漏れがあった場合は自殺
    }
        
    Datasize = atoi(argv[2]);
    if(Datasize>MaxDataSize){
        printf("ERROR! too many datasize\n");
        return 1;  // 指定漏れがあった場合は自殺
    }
    cout<<"radix="<<kisuu<<endl;
    cout<<"thread::hardware_concurrency()="<<thread::hardware_concurrency()<<endl;

   

    cout<<"creating dataset..."<<flush;

    makeDataset();
    vector<ll> D(Dataset);

    cout<<" finish!"<<endl;cout<<endl;
    
    cout<<"std::sort() is running..."<<flush;
    std::chrono::system_clock::time_point start,end;
    start = std::chrono::system_clock::now();
    
    std::sort(D.begin(),D.end());
    
    end = std::chrono::system_clock::now();
    double elapsed = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/1000.0);

    cout<<" finish!"<<endl;

    printf("std::sort time %lf[ms]\n",elapsed);
    cout<<endl;
    
       
    cout<<"PARADIS is running..."<<flush;
    start = std::chrono::system_clock::now();
    //sortしたい目的の配列,levelの数,次のlevelに渡すindexの配列,levelの深さ
    omp_set_nested(1);
    RadixSort<ll,3>(Dataset,Datasize,0,threadNum);
    end = std::chrono::system_clock::now();


    
    if(!issorted(Dataset)){
        std::cerr<<"Not sorted"<<std::endl;
    }else{
      //std::cout<<"sorted!! good job!"<<std::endl;
      cout<<" finish!"<<endl;
    }
    

    elapsed = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/1000.0);

    printf("paradis time %lf[ms]\n",elapsed);
}

