#include "itensor/all.h"
#include "itensor/util/print_macro.h"

using namespace itensor;
using std::vector;

int main()
{
  int N = 50; //number of sites
  Real tstep = 0.02; //time step (smaller is generally more accurate)
  Real ttotal = 1.0; //total time to evolve
  Real cutoff = 1E-8; //truncation error cutoff when restoring MPS form
  int steps=ttotal/tstep;
  //Define a site set object "sites" which lets us
  //easily obtain Site indices defining our Hilbert space
  //and S=1/2 single-site operators
  auto sites = SpinHalf(N);

  //Make initial MPS psi to be in the Neel state
  auto state = InitState(sites);
  for(auto j : range1(N))
    {
      state.set(j,j%2==1?"Up":"Dn");
    }
  auto psi = MPS(state);

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
  file.open("data/overlap_trotter.csv", std::ofstream::trunc);
  file<<"#time,real,imag"<<std::endl;
  //Time evolve, overwriting psi when done
  for(int i =0; i<steps; i++)
    {
      Print(innerC(psi,psi0));
      std::complex<double> overlap=innerC(psi,psi0);
      file<<i*tstep<<","<<overlap.real() << ","<< overlap.imag()<<std::endl;    
      gateTEvol(gates,tstep,tstep,psi,{"Cutoff=",cutoff,"Verbose=",true});

  printfln("Maximum MPS bond dimension after time evolution is %d",maxLinkDim(psi));

  std::cout<< "time "<<i*tstep<<std::endl;
    }

 

  return 0;
}
