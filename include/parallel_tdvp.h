#pragma once
#include "itensor/mps/dmrg.h"
#include "itensor/util/parallel.h"
#include "partition.h"
#include "parallel_dmrg.h"


namespace itensor {

/* function to gather the maximum value of a double distributed onto different processes */
double
gatherMax_val(Environment const& env, double& value)
    { 
    if(env.nnodes() == 1) return value;
    const int root = 0;
    
    if(env.firstNode())
        { 
        for(int i = 1; i < env.nnodes(); i++)
            {
            MailBox mailbox(env,i);
            double tmp;
            mailbox.receive(tmp);
	    if(tmp<0)
	      {std::cout<< "something wrong in gatherMax"<<std::endl;}
	    value=std::max(tmp, value);

            }
        }
    else
        {
        MailBox mailbox(env,root);
        mailbox.send(value);
	value=-1;
        }
    return value;
    }
/* function to gather the MPS distributed onto different processes */
/* This function can definetly be optimized, alternatively, one can contract parts of the MPS on different processes firsr */
  void gather_vector( Environment const& env,  Partition & P, MPS& psi,    std::vector<ITensor> & Vs)
  {

   
        
	   for(int b = 2; b <= P.Nb(); ++b)
                {
		  // for(auto j : range1(P.begin(b),P.end(b)))
		  //{
		   ITensor X;
		   bool edge=true;
		         for(auto j : range1(P.begin(b),P.end(b)))
		  {
		    //	    std::cout<< "node "<< env.rank()<<std::endl;	    
	    	   if(env.rank()+1==b)
	    	     {

		       MailBox mbox(env, 0);	  
		       	 
	    	  X=psi.A(j);
 if(j==P.begin(b+1)-1)
		    {
		      X*=Vs.at( b);
	
		    }
		  else{
	
		  }
				  	    mbox.send(X);

	    	 
	    	  }
		   	
	   
	    	   	   if(env.firstNode())
	    	     {
		   
		         MailBox mbox(env, b-1);
	
	      mbox.receive(X);
	 

		  if(j==P.begin(2))
		    {
		      psi.Aref(j)=Vs.at(b-1)*X;

	
		    }
		  else{
		     psi.Aref(j)=X;

		  }

	      
	     
		   	
	    }
		  }
	}
  }



  /* Function to time evolve a MPS */

   template < typename HamT, typename step_type>
void 
Act(Environment const& env,
            Partition const& P,
            MPS & psi,
            std::vector<ITensor> & Vs,
            HamT & PH,
            Sweeps const& sweeps,
            Observer & obs,
	    step_type t,
    Args args, bool Testing)
    {

    using EdgeType = stdx::decay_t<decltype(PH.L())>;
    //Maximum number of Davidson iterations at boundaries
    auto boundary_niter = args.getInt("BoundaryIter",sweeps.niter(1)+6);
    bool DoNormalize=args.getBool("DoNormalize",false);
    auto b = env.rank()+1;
    auto jl = P.begin(b);
    auto jr = P.end(b);

    Real energy = 0.;
  ITensor phi0,phi1;
    auto psw = ParallelSweeper(jl,jr,b);

    MailBox mboxL,mboxR;
    if(not env.firstNode()) mboxL = MailBox(env,env.rank()-1);
    if(not env.lastNode())  mboxR = MailBox(env,env.rank()+1);
    if(Testing)
      {
    psi.leftLim(jl);
    psi.rightLim(jr);
    Testing=false;
    psi.position(psw.j);
      }


    //Include MPINode number to use in the observer
    args.add("BlockStart",jl);
    args.add("BlockEnd",jr);
    args.add("MPINode",env.rank()+1);
    int jmid = (jr-jl)/2.;
    const int numCenter = 2;
    //args.getInt("NumCenter",2);
     if(numCenter != 2)
       {std::cout<< "does not work for tdvp1 yet!!"<<std::endl;
	 return;}

    for(int sw = 1; sw <= sweeps.nsweep(); ++sw)
        {
 args.add("Sweep",sw);
        args.add("Cutoff",sweeps.cutoff(sw));
        args.add("MinDim",sweeps.mindim(sw));
        args.add("MaxDim",sweeps.maxdim(sw));
        args.add("Noise",sweeps.noise(sw));
        args.add("MaxIter",sweeps.niter(sw));

         for(psw.newSweep(); psw.doingFull(); ++psw)
             {
            auto j = psw.j;
            auto dir = psw.dir();

	        PH.numCenter(numCenter);
            PH.position(j,psi);

             phi1 = psi.A(j)*psi.A(j+1);


	      applyExp(PH,phi1,t/2,args);
      
/* 		 if(DoNormalize) */
	      
/* 		   { */
/* phi1 /= norm(phi1); */
/* 		   } */
	     auto energy=0;
            auto spec = psi.svdBond(j,phi1,dir,PH,args);
            
            if(env.rank()+1 == env.nnodes()/2
            && dir == Fromright
            && j == jmid)
                {
                printfln("Truncation error for sweep %d (node %d) at site %d is terr=%.12f",
                         sw,b,j,spec.truncerr());
                }

            args.add("AtBond",j);
            args.add("AtBoundary",false);
            args.add("Energy",energy);
            args.add("Direction",dir);
            obs.measure(args);

	     if((psw.ha == 1 && psw.j+numCenter-1 != length(psi) && dir==Fromleft) || (psw.ha == 2 && psw.j != 1 && dir==Fromright))

                {
                auto b1 = (psw.ha == 1 ? psw.j+1 : psw.j);
	
		 phi0 = psi.A(b1);

			 	      PH.numCenter(numCenter-1);
				       PH.position(b1,psi);
 
				       applyExp(PH,phi0,-t/2,args);
				
                 /* if(DoNormalize) */
		 /*   {		       	          phi0 /= norm(phi0); */
		 /*   } */
              
		        psi.ref(b1) = phi0;
              
		}
	     else if((psw.ha == 1 && psw.j != 1  && dir==Fromright) || (psw.ha == 2 && psw.j+numCenter-1 != length(psi) && dir==Fromleft))
                 {
		   auto b1 = (psw.ha == 1 ?  psw.j : psw.j+1);
	
		    phi0 = psi(b1);
        
          
 
		     PH.numCenter(numCenter-1);
		       PH.position(b1,psi);
 
		           applyExp(PH,phi0,-t/2,args);
		
                 /* if(DoNormalize) */
		 /*   {	    phi0 /= norm(phi0); */
		 /*   } */
                
      
		        psi.ref(b1) = phi0;
  
         
	       	}
		
            if(psw.atRight() && dir==Fromleft && bool(mboxR))
                {
	

     				auto n = j+1;

  PH.numCenter(numCenter);
          PH.position(n,psi);

                  Boundary<EdgeType> B; 
                  mboxR.receive(B); 
                 psi.Aref(n+1) = B.A; 
                  PH.R(B.HH); 
                  B = Boundary<EdgeType>(); //to save memory 

               auto& V = Vs.at(b); 

                 auto phi = psi.A(n)*V*psi.A(n+1); 
	
 			applyExp(PH,phi,t/2,args); 

 /* if(DoNormalize)  */
 /* 		  {  */
 /* 		     				phi /= norm(phi);  */
 /* 		  }  */
		
 				auto spec = isvd(phi,psi.Aref(n),V,psi.Aref(n+1), args); 
			
				////////////////////////////////////////////
				// time evolve twice for small bond dimension on bonds

                 phi = psi.A(n)*V*psi.A(n+1); 
	
 			applyExp(PH,phi,t/2,args); 

 /* if(DoNormalize)  */
 /* 		  {  */
 /* 		     				phi /= norm(phi);  */
 /* 		  }  */
		
 			        spec = isvd(phi,psi.Aref(n),V,psi.Aref(n+1), args); 
  /////////////////////////////////////////////////
		
                B.HH = PH.L();
                B.UU = psi.A(n)*V;
                B.A = psi.A(n+1);
                B.energy = energy;
                mboxR.send(B);

                psi.Aref(n+1) *= V;
                psi.rightLim(n+1);





	
		phi0 = psi.A(n);
	
			 	      PH.numCenter(numCenter-1);
				       PH.position(n,psi);
 
				       applyExp(PH,phi0,-t/2,args);
				       //
                /* if(DoNormalize) */
		/*   { */
		/*       		        	          phi0 /= norm(phi0); */
                
		/*   } */
               
		        psi.ref(n) = phi0;
                  } 
              else if(psw.atLeft() && dir==Fromright && bool(mboxL)) 
                  { 


 			    PH.numCenter(numCenter); 
		
			    
                 auto n = j-1; 
	
   PH.numCenter(numCenter); 
                 PH.position(n,psi); 

                Boundary<EdgeType> B;
                B.A = psi.A(n+1);
                B.HH = PH.R();
                mboxL.send(B);

                mboxL.receive(B);
                PH.L(B.HH);
                PH.shift(n,Fromleft,B.UU);
                psi.Aref(n+1) = B.A;
                energy = B.energy;
		phi0 = psi.A(n+1);
      PH.numCenter(numCenter-1);
      PH.position(n+1,psi);
      
      applyExp(PH,phi0,-t/2,args);
       /* if(DoNormalize) */
       /* 		  { */
       /* 		        phi0 /= norm(phi0); */
       /* 		  } */
                // if(numCenter == 2)
                //     {
		        psi.ref(n+1) = phi0;
		// printfln("Node %d done with boundary step, energy %.5f -> %.5f",b,prev_energy,energy);
                  } 
              } 

    //     if(obs.checkDone(args)) break;
           }

    //  printfln("Block %d final energy = %.12f",b,energy);

    return;

    } 

  ///////////////////////////////////////////////////////////////////////////
  // Next follows a DMRG style interface, I have not tested in yet
/* Real  */
/* parallel_tdvp(Environment const& env, */
/*               MPS & psi, */
/*               MPO const& H, */
/*               Sweeps const& sweeps, */
/*               Observer & obs, */
/*               Args args = Args::global()) */
/*     { */
/*     Partition P; */
/*     std::vector<ITensor> Vs; */
/*     splitWavefunction(env,psi,P,Vs,args); */
/*     //  gather_vector(env,P, psi); */
/*     auto PH = computeHEnvironment(env,P,psi,Vs,H,args); */
/*     double dt=0.1; */
/*       return ptdvpWorker(env,P,psi,Vs,PH,sweeps,obs,dt, args); */
/*     } */

/* // */
/* // parallel_tdvp with single MPO or IQMPO */
/* // */

/* Real  */
/* parallel_tdvp(Environment const& env, */
/*               MPS & psi, */
/*               MPO const& H, */
/*               Sweeps const& sweeps, */
/*               Args const& args = Args::global()) */
/*     { */
/*     Observer obs; */
/*     return parallel_tdvp(env,psi,H,sweeps,obs,args); */
/*     } */



/* // */
/* // parallel_tdvp with an (implicit) sum of MPOs or IQMPOs */
/* // and an observer object */
/* // */

/* Real  */
/* parallel_tdvp(Environment const& env, */
/*               MPS & psi, */
/*               std::vector<MPO> const& Hset, */
/*               Sweeps const& sweeps, */
/*               Observer & obs, */
/*               Args args = Args::global()) */
/*     { */
/*     Partition P; */

/*      std::vector<ITensor> Vs; */
/*     splitWavefunction(env,psi,P,Vs,args); */
/*     auto PH = computeHEnvironment(env,P,psi,Vs,Hset,args); */
/*     double dt=0.1; */
/*     return ptdvpWorker(env,P,psi,Vs,PH,sweeps,obs,dt, args); */
/*     } */

/* // */
/* // parallel_tdvp with an (implicit) sum of MPOs or IQMPOs */
/* // */

/* Real  */
/* parallel_tdvp(Environment const& env, */
/*               MPS & psi, */
/*               std::vector<MPO> const& Hset, */
/*               Sweeps const& sweeps, */
/*               Args const& args = Args::global()) */
/*     { */
      
/*     Observer obs; */
/*     return parallel_tdvp(env,psi,Hset,sweeps,obs,args); */
/*     } */


/*  template < typename HamT, typename step_type> */
/* Real  */
/* ptdvpWorker(Environment const& env, */
/*             Partition const& P, */
/*             MPS & psi, */
/*             std::vector<ITensor> & Vs, */
/*             HamT & PH, */
/*             Sweeps const& sweeps, */
/*             Observer & obs, */
/* 	    step_type t, */
/*             Args args) */
/*     { */

/*     using EdgeType = stdx::decay_t<decltype(PH.L())>; */
/*     //Maximum number of Davidson iterations at boundaries */
/*     auto boundary_niter = args.getInt("BoundaryIter",sweeps.niter(1)+6); */
/*     bool DoNormalize=args.getBool("DoNormalize",false); */
/*     auto b = env.rank()+1; */
/*     auto jl = P.begin(b); */
/*     auto jr = P.end(b); */

/*     Real energy = 0.; */
/*   ITensor phi0,phi1; */
/*     auto psw = ParallelSweeper(jl,jr,b); */

/*     MailBox mboxL,mboxR; */
/*     if(not env.firstNode()) mboxL = MailBox(env,env.rank()-1); */
/*     if(not env.lastNode())  mboxR = MailBox(env,env.rank()+1); */

/*     psi.leftLim(jl); */
/*     psi.rightLim(jr); */
/*     psi.position(psw.j); */

/*     //Include MPINode number to use in the observer */
/*     args.add("BlockStart",jl); */
/*     args.add("BlockEnd",jr); */
/*     args.add("MPINode",env.rank()+1); */
/*     int jmid = (jr-jl)/2.; */
/*     const int numCenter = 2; */
/*     //args.getInt("NumCenter",2); */
/*      if(numCenter != 2) */
/*        {std::cout<< "does not work for tdvp1 yet!!"<<std::endl; */
/* 	 return 0;} */
/*     // if(numCenter != 1) */
/*     //     args.add("Truncate",args.getBool("Truncate",true)); */
/*     // else */
/*     //     args.add("Truncate",args.getBool("Truncate",false)); */
/*     //  args.add("DoNormalize",true); */
/*     for(int sw = 1; sw <= sweeps.nsweep(); ++sw) */
/*         { */
/*         args.add("Sweep",sw); */
/*         args.add("Cutoff",sweeps.cutoff(sw)); */
/*         args.add("MinDim",sweeps.mindim(sw)); */
/*         args.add("MaxDim",sweeps.maxdim(sw)); */
/*         args.add("Noise",sweeps.noise(sw)); */
/*         args.add("MaxIter",sweeps.niter(sw)); */

/*         if(!PH.doWrite() */
/*            && args.defined("WriteM") */
/*            && sweeps.maxdim(sw) >= args.getInt("WriteM")) */
/*             { */
/*             printfln("\nNode %d turning on write to disk, write_dir = %s", */
/*                      b,args.getString("WriteDir","./")); */
/*             //psi.doWrite(true); */
/*             PH.doWrite(true); */
/*             } */

/*            printfln("Doing sweep %d for node %d (maxm=%d, cutoff=%.0E, mindim=%d)",sw,b,sweeps.maxdim(sw),sweeps.cutoff(sw),sweeps.mindim(sw)); */

/*          for(psw.newSweep(); psw.doingFull(); ++psw) */
/*              { */
/*             auto j = psw.j; */
/*             auto dir = psw.dir(); */
/*             //printfln("%d j = %d (%d,%d) %s",b,j,jl,jr,dir==Fromleft?"Fromleft":"Fromright"); */
/* 	    // if(env.firstNode()) */
/* 	    //   { */
/* 	    //	    std::cout<< "node rank " << env.rank()+1<<" act on site  "<< j<< " and "<< j+1<< " and ha = "<<psw.ha <<std::endl; */
/* 		//	      } */
/* 	      //   if(env.lastNode()) */
/* 	      // { */
/* 	      // 	std::cout<< "node las "<< j <<std::endl; */
/* 	      // } */
/* 	        PH.numCenter(numCenter); */
/*             PH.position(j,psi); */

/*              phi1 = psi.A(j)*psi.A(j+1); */

/* 	     //  energy = davidson(PH,phi1,args); */
/* 	      applyExp(PH,phi1,t/2,args); */
/*             //if(env.rank()+1 == 1) printfln("%s j = %d energy = %.10f",dir==Fromleft?"->":"<-",j,energy); */
/* 		 if(DoNormalize) */
	      
/* 		   { */
/* phi1 /= norm(phi1); */
/* 		   } */
/* 	     auto energy=0; */
/*             auto spec = psi.svdBond(j,phi1,dir,PH,args); */
            
/*             if(env.rank()+1 == env.nnodes()/2  */
/*             && dir == Fromright  */
/*             && j == jmid) */
/*                 { */
/*                 printfln("Truncation error for sweep %d (node %d) at site %d is terr=%.12f", */
/*                          sw,b,j,spec.truncerr()); */
/*                 } */

/*             args.add("AtBond",j); */
/*             args.add("AtBoundary",false); */
/*             args.add("Energy",energy); */
/*             args.add("Direction",dir); */
/*             obs.measure(args); */

/* 	     if((psw.ha == 1 && psw.j+numCenter-1 != length(psi) && dir==Fromleft) || (psw.ha == 2 && psw.j != 1 && dir==Fromright)) */

/*                 { */
/*                 auto b1 = (psw.ha == 1 ? psw.j+1 : psw.j); */
/* 		//	std::cout<< " XXX   rank "<<env.rank()+1<< " and j " << psw.j<<"  and "<< b1 << std::endl; */
/* 		//     if(numCenter == 2) */
/* 		//  { */
/* 		 phi0 = psi.A(b1); */
/* 		    //  } */
/* 		 //			 std::cout<< "j "<< j << " b1 "<<b1 <<std::endl; */
/* 			 // PH.position(b1,psi); */
/* 			 	      PH.numCenter(numCenter-1); */
/* 				       PH.position(b1,psi); */
 
/* 				       applyExp(PH,phi0,-t/2,args); */
/* 				       //					  energy = davidson(PH,phi0,args); */
/*                  if(DoNormalize) */
/* 		   {		       	          phi0 /= norm(phi0); */
/* 		   } */
/*                 // if(numCenter == 2) */
/*                 //     { */
/* 		        psi.ref(b1) = phi0; */
/*                 //     } */
/* 		} */
/* 	     else if((psw.ha == 1 && psw.j != 1  && dir==Fromright) || (psw.ha == 2 && psw.j+numCenter-1 != length(psi) && dir==Fromleft)) */
/*                  { */
/* 		   auto b1 = (psw.ha == 1 ?  psw.j : psw.j+1); */
/* 		   //  	std::cout<< " YYY  rank "<<env.rank()+1<< " and j " << psw.j<<"  and "<< b1 << std::endl; */
/* 		//          if(numCenter == 2) */
/*                 //     { */
/* 		    phi0 = psi(b1); */
/*                 //     } */
          
 
/* 		     PH.numCenter(numCenter-1); */
/* 		       PH.position(b1,psi); */
 
/* 		           applyExp(PH,phi0,-t/2,args); */
/* 			   //energy = davidson(PH,phi0,args); */
/*                  if(DoNormalize) */
/* 		   {	    phi0 /= norm(phi0); */
/* 		   } */
                
/*                 // if(numCenter == 2) */
/*                 //     { */
/* 		        psi.ref(b1) = phi0; */
/*                 //     } */
  
/* 	       	} */
		
/*             if(psw.atRight() && dir==Fromleft && bool(mboxR)) */
/*                 { */
/* 		  //	  std::cout<< "WWW rank "<<env.rank()+1<< " j "<<j << " n " << j+1<<std::endl; */
/*     //             auto prev_energy = energy; */
/*     //             printfln("Node %d communicating with right, boundary_niter=%d",b,boundary_niter); */
/*      				auto n = j+1; */
/* 					std::cout<< " n hete "<< n << " at  "<<env.rank() <<std::endl; */
/*   PH.numCenter(numCenter); */
/*         PH.position(n,psi); */

/*                  Boundary<EdgeType> B; */
/*                  mboxR.receive(B); */
/*                  psi.Aref(n+1) = B.A; */
/*                  PH.R(B.HH); */
/*                  B = Boundary<EdgeType>(); //to save memory */

/*                  auto& V = Vs.at(b); */

/*                 auto phi = psi.A(n)*V*psi.A(n+1); */
/* 		//   */
/* 		//	std::cout<< "Pri 11111"<< V<<std::endl; */
/* 		//energy = davidson(PH,phi,{args,"MaxIter=",boundary_niter}); */
/* 			applyExp(PH,phi,t/2,args); */
/* 			applyExp(PH,phi,t/2,args); */
/* 		//energy = davidson(PH,phi,args); */
/* 			//			std::cout<< "norm s "<< norm(phi)<<std::endl;  */
/* if(DoNormalize) */
/* 		  { */
/* 				phi /= norm(phi); */
/* 		  } */
/* 			//	std::cout<< "norm s "<< norm(phi)<<std::endl; */
/* 				auto spec = isvd(phi,psi.Aref(n),V,psi.Aref(n+1), args); */
/* 				//	std::cout<< "Pri "<< V<<std::endl; */

		
/*                 B.HH = PH.L(); */
/*                 B.UU = psi.A(n)*V; */
/*                 B.A = psi.A(n+1); */
/*                 B.energy = energy; */
/*                 mboxR.send(B); */

/*                 psi.Aref(n+1) *= V; */
/*                 psi.rightLim(n+1); */





/* 		// new part */
/* 		phi0 = psi.A(n); */
/* 		    //  } */
/* 		 //			 std::cout<< "j "<< j << " b1 "<<b1 <<std::endl; */
/* 			 // PH.position(b1,psi); */
/* 			 	      PH.numCenter(numCenter-1); */
/* 				       PH.position(n,psi); */
 
/* 				       applyExp(PH,phi0,-t/2,args); */
/* 				       //					  energy = davidson(PH,phi0,args); */
/*                 if(DoNormalize) */
/* 		  { */
/* 		       		        	          phi0 /= norm(phi0); */
                
/* 		  } */
/*                 // if(numCenter == 2) */
/*                 //     { */
/* 		        psi.ref(n) = phi0; */




/*                 // args.add("AtBond",n); */
/*                 // args.add("AtBoundary",true); */
/*                 // args.add("Energy",energy); */
/*                 // obs.measure(args); */

/* 		//  printfln("Node %d done with boundary step, energy %.5f -> %.5f",b,prev_energy,energy); */
/*                  } */
/*              else if(psw.atLeft() && dir==Fromright && bool(mboxL)) */
/*                  { */
/* 		   //	  std::cout<< "LLL rank "<<env.rank()+1<< " j "<<j << " n " << j-1<<std::endl; */
/* 			    PH.numCenter(numCenter); */
/* 			    //	      PH.position(n,psi); */
/*     //             auto prev_energy = energy; */
/*     //             printfln("Node %d communicating with left, boundary_niter=%d",b,boundary_niter); */
			    
/*                 auto n = j-1; */
	
/*   PH.numCenter(numCenter); */
/*                 PH.position(n,psi); */

/*                 Boundary<EdgeType> B; */
/*                 B.A = psi.A(n+1); */
/*                 B.HH = PH.R(); */
/*                 mboxL.send(B); */

/*                 mboxL.receive(B); */
/*                 PH.L(B.HH); */
/*                 PH.shift(n,Fromleft,B.UU); */
/*                 psi.Aref(n+1) = B.A; */
/*                 energy = B.energy; */
/* 		phi0 = psi.A(n+1); */
/*       PH.numCenter(numCenter-1); */
/*       PH.position(n+1,psi); */
      
/*       applyExp(PH,phi0,-t/2,args); */
/*        if(DoNormalize) */
/* 		  { */
/* 	       	          phi0 /= norm(phi0); */
/* 		  } */
/*                 // if(numCenter == 2) */
/*                 //     { */
/* 		        psi.ref(n+1) = phi0; */
/* 		// printfln("Node %d done with boundary step, energy %.5f -> %.5f",b,prev_energy,energy); */
/*                  } */
/*              } */

/*     //     if(obs.checkDone(args)) break; */
/*            } */

/*     printfln("Block %d final energy = %.12f",b,energy); */

/*     return energy; */

/*     } // ptdvpWorker */






} //namespace itensor


