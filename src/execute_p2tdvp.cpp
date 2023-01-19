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
  int N = 10; //number of sites
  Real tstep = 0.01; //time step (smaller is generally more accurate)
  Real ttotal = 1.0; //total time to evolve
  Real cutoff = 1E-9; //truncation error cutoff when restoring MPS form
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
  // make a non product state the initial state
  auto psi_in = MPS(state);
  auto sweeps_gs = Sweeps(5); //number of sweeps is 5
  sweeps_gs.maxdim() = 10,20,100,100,200;
  sweeps.cutoff() = 1E-10;
  auto ampo_in = AutoMPO(sites);
  for(int j = 1; j < N; ++j)
    {
      ampo_in += 0.5,"S+",j,"S-",j+1;
      ampo_in += 0.5,"S-",j,"S+",j+1;
      
    }
  auto H_in = toMPO(ampo_in);
  auto [energy,psi_gs] = dmrg(H_in,psi_in,sweeps_gs);

  println("Ground State Energy = ",energy);


  auto ampo = AutoMPO(sites);
  for(int j = 1; j < N; ++j)
    {
      ampo += 0.5,"S+",j,"S-",j+1;
      ampo += 0.5,"S-",j,"S+",j+1;
      ampo +=     "Sz",j,"Sz",j+1;
    }
  H = toMPO(ampo);

 psi=psi_gs;
  //Save initial state;
   psi0 = psi;
 // write to file
file.open("data/correlfunc_p2tdvp.csv", std::ofstream::trunc);
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
		       int site_1=N/2;
		       int site_2=site_1+1;
		      auto op_1 = op(sites,"Sz",site_1);
		       auto op_2 = op(sites,"Sz",site_2);
		       //'gauge' the MPS to site i
		       //any 'position' between i and j, inclusive, would work here
		       psi_new.position(site_1); 

		       //Create the bra/dual version of the MPS psi
		       auto psidag_new = dag(psi_new);

		       //Prime the link indices to make them distinct from
		       //the original ket links
		       psidag_new.prime("Link");

		       //index linking i-1 to i:
		       auto li_1 = leftLinkIndex(psi_new,site_1);

		       auto C = prime(psi_new(site_1),li_1)*op_1;
		       C *= prime(psidag_new(site_1),"Site");
		       for(int k = site_1+1; k < site_2; ++k)
			 {
			   C *= psi_new(k);
			   C *= psidag_new(k);
			 }
		       //index linking j to j+1:
		       auto li_2 = rightLinkIndex(psi_new,site_2);

		       C *= prime(psi_new(site_2),li_2)*op_2;
		       C *= prime(psidag_new(site_2),"Site");

		       auto result = eltC(C);
		       Print(innerC(psi_new,psi0));
		       Print(result);
		      
		       file<<i*tstep<<","<<result.real() << ","<< result.imag()<<std::endl;    
		       std::cout<< i*tstep<<std::endl;
   printfln("Maximum MPS bond dimension after time evolution is %d",maxLinkDim(psi_new));
 		     }
 Act(env,P,psi,Vs,PH,sweeps,obs,-tstep*Cplx_i, argsMPS, first);

      	   	if(first)
      	  {
   
      	    first=false;
      	  }
    }

 

  return 0;
}
