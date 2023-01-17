#include "itensor/all.h"
#include "itensor/util/print_macro.h"
#include "parallel_tdvp.h"

using namespace itensor;
using std::vector;

int main(int argc, char* argv[])
{
      Environment env(argc,argv);

         parallelDebugWait(env);
 
  MPS psi;
  MPS psi0;
  MPO H;
  int N = 50; //number of sites
  Real tstep = 0.02; //time step (smaller is generally more accurate)
  Real ttotal = 1.0; //total time to evolve
  Real cutoff = 1E-8; //truncation error cutoff when restoring MPS form
  int steps=ttotal/tstep;
    Sweeps sweeps;
	 Args argsMPS=itensor::Args("Cutoff=",cutoff,"Truncate", true,"DoNormalize",false,"SVDMethod=","gesdd");
 auto sites = SpinHalf(N);
  std::ofstream file;
  //Define a site set object "sites" which lets us
  //easily obtain Site indices defining our Hilbert space
  //and S=1/2 single-site operators
         if(env.firstNode())
        {
    
     sweeps = Sweeps(1);
    sweeps.cutoff() = cutoff;
    sweeps.niter() = 10;
    


  //Make initial MPS psi to be in the Neel state
  auto state = InitState(sites);
  for(auto j : range1(N))
    {
      state.set(j,j%2==1?"Up":"Dn");
    }
  psi = MPS(state);
  auto ampo = AutoMPO(sites);
  for(int j = 1; j < N; ++j)
    {
      ampo += 0.5,"S+",j,"S-",j+1;
      ampo += 0.5,"S-",j,"S+",j+1;
      ampo +=     "Sz",j,"Sz",j+1;
    }
  H = toMPO(ampo);


  //Save initial state;
   psi0 = psi;
 // write to file
file.open("data/overlap_p2tdvp.csv", std::ofstream::trunc);
  file<<"#time,real,imag"<<std::endl;	
}
 

  
 env.broadcast(sites,H,psi,sweeps);
 	 Partition P;

      std::vector<ITensor> Vs;
     splitWavefunction(env,psi,P,Vs, argsMPS);
 auto PH = computeHEnvironment(env,P,psi,Vs,H);
  //Time evolve, overwriting psi when done
         Observer obs;
	 bool first=true;
   for(int i =0; i<steps; i++)
    {
//       //Print(innerC(psi,psi0));
      
      
 MPS psi_new=psi;
   gather_vector(env,P, psi_new, Vs);
    if(env.firstNode())
       	    	     {
  Print(innerC(psi_new,psi0));
      std::complex<double> overlap=innerC(psi_new,psi0);
       file<<i*tstep<<","<<overlap.real() << ","<< overlap.imag()<<std::endl;    
   std::cout<< i*tstep<<std::endl;
 		     }
 Act(env,P,psi,Vs,PH,sweeps,obs,-tstep*Cplx_i, argsMPS, first);
//   printfln("Maximum MPS bond dimension after time evolution is %d",maxLinkDim(psi));
      	   	if(first)
      	  {
   
      	    first=false;
      	  }
    }

 

  return 0;
}
