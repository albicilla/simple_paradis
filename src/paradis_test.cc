#include "paradis.h"


//Sort 対象データセット
void makeDataset(long long Datasize,std::vector<int>& Dataset){ 
  Dataset.clear();

  std::uniform_int_distribution<int> dist(0.0,INT_MAX);

  std::random_device seed_gen;
  std::default_random_engine engine(seed_gen());

  Dataset.resize(Datasize);
    
  for(int i=0;i<Datasize;i++){
    int result = dist(engine);
    Dataset[i]=(result);
  }  
}


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



int main(int argc,char ** argv){
  if(argc<2){
    printf("ERROR! set thread num\n");
    return 1;
  }
  
  int threadNum = atoi(argv[1]);
  
  if(threadNum>MaxThreadNum){
    printf("ERROR! too many thread num\n");
    return 1;
  }
  
  if(argc < 3){
    printf("ERROR! set datasize\n");
    return 1;
  }

  long long Datasize = atoi(argv[2]);
  if(Datasize>MaxDataSize){
    printf("ERROR! too many datasize\n");
    return 1;
  }
  
  std::vector<int> Dataset;

  makeDataset(Datasize,Dataset);

  std::cout<<"PARADIS is running..."<<std::flush;
  auto start = std::chrono::system_clock::now();
  omp_set_nested(1);
  PARADIS(Dataset.begin(),Dataset.end(),threadNum);
  auto end = std::chrono::system_clock::now();


    
  if(!std::is_sorted(Dataset.begin(),Dataset.end())){
    std::cerr<<"Not sorted"<<std::endl;
  }else{
    std::cout<<" finish!"<<std::endl;
  }
    

  double elapsed = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/1000.0);

  printf("paradis time %lf[ms]\n",elapsed);
  
}
