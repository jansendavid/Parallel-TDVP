#include "itensor/all.h"
#include "itensor/util/print_macro.h"

using namespace itensor;
using std::vector;

int main()
{

  Real tstep = 0.01; //time step (smaller is generally more accurate)
  Real ttotal = 6.0; //total time to evolve
  Real cutoff = 1E-10; //truncation error cutoff when restoring MPS form
  int steps=ttotal/tstep;


  //Make initial MPS psi to be in the Neel state

  auto f_mps = h5_open("initial_state.h5",'r'); //open HDF5 file in read 'r' mode
  auto psi = h5_read<MPS>(f_mps,"MPS_initial");

  int N = length(psi); //number of sites
  //Define a site set object "sites" which lets us
  //easily obtain Site indices defining our Hilbert space
  //and S=1/2 single-site operators
  auto sites = SpinHalf(siteInds(psi));
  auto state = InitState(sites);
  //Create a std::vector (dynamically sizeable array)
  //to hold the Trotter gates
  auto gates = vector<BondGate>();

  //Create the gates exp(-i*tstep/2*hterm)
  //and add them to gates
  for(int b = 1; b <= N-1; ++b)
    {
      auto hterm = op(sites,"Sz",b)*op(sites,"Sz",b+1);
      hterm += 0.5*op(sites,"S+",b)*op(sites,"S-",b+1);
      hterm += 0.5*op(sites,"S-",b)*op(sites,"S+",b+1);

      auto g = BondGate(sites,b,b+1,BondGate::tReal,tstep/2.,hterm);
      gates.push_back(g);
    }
  //Create the gates exp(-i*tstep/2*hterm) in reverse
  //order (to get a second order Trotter breakup which
  //does a time step of "tstep") and add them to gates
  for(int b = N-1; b >= 1; --b)
    {
      auto hterm = op(sites,"Sz",b)*op(sites,"Sz",b+1);
      hterm += 0.5*op(sites,"S+",b)*op(sites,"S-",b+1);
      hterm += 0.5*op(sites,"S-",b)*op(sites,"S+",b+1);

      auto g = BondGate(sites,b,b+1,BondGate::tReal,tstep/2.,hterm);
      gates.push_back(g);
    }

  //Save initial state;
  auto psi0 = psi;
  // write to file
  std::ofstream file;
  file.open("data/correlfunc_trotter.csv", std::ofstream::trunc);
  file<<"#time,real,imag"<<std::endl;
  //Time evolve, overwriting psi when done
  for(int i =0; i<steps; i++)
    {


      //Make the operators you want to measure
      int site_1=N/2;
		       int site_2=site_1+1;
		       auto op_1 = op(sites,"Sz",site_1);
		      auto op_2 = op(sites,"Sz",site_2);

      //'gauge' the MPS to site i
      //any 'position' between i and j, inclusive, would work here
      psi.position(site_1); 

      //Create the bra/dual version of the MPS psi
      auto psidag = dag(psi);

      //Prime the link indices to make them distinct from
      //the original ket links
      psidag.prime("Link");

      //index linking i-1 to i:
      auto li_1 = leftLinkIndex(psi,site_1);

      auto C = prime(psi(site_1),li_1)*op_1;
      C *= prime(psidag(site_1),"Site");
      for(int k = site_1+1; k < site_2; ++k)
	{
	  C *= psi(k);
	  C *= psidag(k);
	}
      //index linking j to j+1:
      auto li_2 = rightLinkIndex(psi,site_2);

      C *= prime(psi(site_2),li_2)*op_2;
      C *= prime(psidag(site_2),"Site");

      auto result = eltC(C); //or eltC(C) if expecting complex

      Print(innerC(psi,psi0));

      file<<i*tstep<<","<<result.real() << ","<< result.imag()<<std::endl;    
      gateTEvol(gates,tstep,tstep,psi,{"Cutoff=",cutoff,"Verbose=",true,  "MaxDim=",1000});

  printfln("Maximum MPS bond dimension after time evolution is %d",maxLinkDim(psi));

  std::cout<< "time "<<i*tstep<<std::endl;
    }

 

  return 0;
}
