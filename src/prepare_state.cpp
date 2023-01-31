#include "itensor/all.h"
#include "itensor/util/print_macro.h"

using namespace itensor;
using std::vector;

int main()
{
  int N = 800; //number of sites



  auto sites = SpinHalf(N);

  //Make initial MPS psi to be in the Neel state
  auto state = InitState(sites);

  auto sweeps_gs = Sweeps(5); //number of sweeps is 5
  sweeps_gs.maxdim() = 10,20,100,100,200;
  sweeps_gs.cutoff() = 1E-10;
  auto ampo_in = AutoMPO(sites);
  for(auto j : range1(N))
    {
      state.set(j,j%2==1?"Up":"Dn");
    }
auto psi_in = MPS(state);
  //auto psi_in = randomMPS(state);
  for(int j = 1; j < N; ++j)
    {
      ampo_in += 0.5,"S+",j,"S-",j+1;
      ampo_in += 0.5,"S-",j,"S+",j+1;
      
    }
  auto H_in = toMPO(ampo_in);
  auto [energy,psi_gs] = dmrg(H_in,psi_in,sweeps_gs);

  println("Ground State Energy = ",energy);
  auto f_mps = h5_open("initial_state.h5",'w'); //open HDF5 file in write 'w' mode
  h5_write(f_mps,"MPS_initial",psi_gs); 
  close(f_mps);


  return 0;
}
