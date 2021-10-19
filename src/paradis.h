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
#include <omp.h>

#define FOR(i,a,b) for(int i=a;i<b;i++)
#define rep(i,b) FOR(i,0,b)

const int MaxThreadNum=224;
const long long MaxDataSize=10000000000;
const long long MaxDataNum=4294967295;
const int MaxKisuu=256;

std::vector<int> Dataset;
long long Datasize;

static const int kRadixBits = 8;
static const size_t kInsertSortThreshold = 0;
static const int kRadixMask = (1 << kRadixBits) - 1;
static const int kRadixBin = 1 << kRadixBits;

template<class D>
inline int determineDigitBucket(int stage,D num){
  return ((num>>(8*stage))&kRadixMask);
}


template< class _Type>
inline void _swap(_Type &a, _Type&b)  {
  _Type temp = b;
  b = a;
  a = temp;
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


template<int kth_byte,class RandomIt>
inline void PARADIS_core(RandomIt s,RandomIt t,RandomIt begin_itr,int processes=1){
    long long cnt[MaxKisuu]={0};

    long long elenum=distance(s,t);
    long long start=distance(begin_itr,s);
    //assert(start>=0);assert(elenum>=0);
    //step1
    //assert(processes>0);
    long long part=elenum/processes;
    long long res=elenum%processes;

    long long localHists[MaxThreadNum][MaxKisuu];
    long long gh[MaxKisuu],gt[MaxKisuu],starts[MaxKisuu],ends[MaxKisuu];
    long long ph[MaxThreadNum][MaxKisuu];
    long long pt[MaxThreadNum][MaxKisuu];

    long long SumCi=elenum;
    long long pfp[processes+1];
    int var_p=processes;

#pragma omp parallel num_threads(processes)
    {
        int th=omp_get_thread_num();
        #pragma omp for
        rep(i,kRadixBin){
            rep(t,processes)localHists[t][i]=0;
        }
        #pragma omp barrier
        #pragma omp for
        for(int i=start;i<start+elenum;i++){
            int digit=determineDigitBucket(kth_byte,*(begin_itr+i));
            localHists[th][digit]++;
        }
        #pragma omp barrier
        #pragma omp for
        for(int i=0;i<kRadixBin;i++){
            for(int j=0;j<processes;j++){
                cnt[i]+=localHists[j][i];
            }
        }
        #pragma omp barrier
        #pragma omp single
        {
            gh[0]=start;
            gt[0]=gh[0]+cnt[0];
            starts[0]=gh[0];
        }
        //step2
        #pragma omp single
        for(int i=1;i<kRadixBin;i++){
            //calc ghi
            gh[i]=gh[i-1]+cnt[i-1];
            //calc gti
            gt[i]=gh[i]+cnt[i];
            starts[i]=gh[i];
        }
        #pragma omp barrier
        //step3
        while(SumCi!=0){
            #pragma omp for
            for(int ii=0;ii<processes;ii++){
                int pID=omp_get_thread_num();
                for(int i=0;i<kRadixBin;i++){
                    long long part=(long long)(gt[i]-gh[i])/(long long)var_p;
                    long long res=(long long)(gt[i]-gh[i])%(long long)(var_p);
                    if(pID<var_p-1){
                        ph[pID][i]=part*pID+gh[i];
                        pt[pID][i]=part*(pID+1LL)+gh[i];
                    }else{
                        ph[pID][i]=part*pID+gh[i];
                        pt[pID][i]=part*(pID+1LL)+gh[i]+res;
                    }
                }
            

                for(int i=0;i<kRadixBin;i++){
                    long long head=ph[pID][i];
                    while(head<pt[pID][i]){
                        auto v=*(begin_itr+head);
                        int k=determineDigitBucket(kth_byte,v);
                        while(k!=i&&ph[pID][k]<pt[pID][k]){
                            _swap(v,*(begin_itr+(int)ph[pID][k]));ph[pID][k]++;
                            k=determineDigitBucket(kth_byte,v);
                        }
                        if(k==i){
                            *(begin_itr+head)=*(begin_itr+ph[pID][i]);head++;
                            *(begin_itr+ph[pID][i])=v;ph[pID][i]++;
                        }else{
                            *(begin_itr+head)=v;head++;
                        }
                    }
                }
            }//end of omp permute
            #pragma omp single
            {
                SumCi=0;
                long long pfpN=kRadixBin/var_p;
                long long pfpM=kRadixBin%var_p;
                pfp[0]=0LL;
                long long pfpMR=0LL;
                for(long long i=1LL;i<var_p+1LL;i++){
                    if(pfpMR<pfpM)pfpMR++;
                    pfp[i]=i*pfpN+pfpMR;
                }
            }
            #pragma omp barrier
            #pragma omp for
            for(int k=0;k<processes;k++){
                for(long long i=pfp[k];i<pfp[k+1];i++){
                    long long tail=gt[i];
                    {
                        for(int pID=0;pID<processes;pID++){
                            long long head=ph[pID][i];
                            while(head<pt[pID][i]&&head<tail){
                                int v=*(begin_itr+head);head++;
                                if(determineDigitBucket(kth_byte,v)!=i){
                                    while(head<=tail){
                                        tail--;
                                        int w=*(begin_itr+tail);
                                        if(determineDigitBucket(kth_byte,w)==i){
                                            *(begin_itr+(head-1))=w;
                                            *(begin_itr+tail)=v;
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
                SumCi-0;
                for(int i=0;i<kRadixBin;i++){
                    SumCi+=(gt[i]-gh[i]);
                }
            }
            #pragma omp barrier
        

        
        }//end of while
    }//end of omp2

    if(kth_byte>0){
#pragma omp parallel num_threads(processes)
        #pragma omp single
        {
            for(int i=0;i<kRadixBin;i++){
                int nextStageThreads=1;
                nextStageThreads=processes*(cnt[i]*(log(cnt[i])/log(kRadixBin))/(elenum*(log(elenum)/log(kRadixBin))));
                if(cnt[i]>64LL){
                    #pragma omp task
		  PARADIS_core<(kth_byte > 0 ? (kth_byte - 1) : 0)>(begin_itr+starts[i],begin_itr+(starts[i]+cnt[i]),begin_itr,std::max(nextStageThreads,1));
                }else if(cnt[i]>1){
                    insert_sort_core_(begin_itr+starts[i],begin_itr+(starts[i]+cnt[i]));
                    //std::sort(begin_itr+starts[i],begin_itr+(starts[i]+cnt[i]));
                }
            }
            #pragma omp taskwait
        }
    }    
}

template<class RandomIt>
inline void PARADIS(RandomIt s,RandomIt t,int threadNum){
    const size_t vsize=sizeof(typename std::iterator_traits<RandomIt>::value_type);
    PARADIS_core<vsize-1>(s,t,s,threadNum);
}

