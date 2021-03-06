
/* This is a code to compute correlators
*/

#include <Grid/Grid.h>
#include <Grid/algorithms/iterative/BlockConjugateGradient.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;


//////////////////////////////////////////////
// Load Fermion into propagator 
//////////////////////////////////////////////

void FermToProp_s(LatticeStaggeredPropagator &Qprop, LatticeStaggeredFermion &psi , const int c)
{
  const int nc = 3 ;  // later use the dimension 

  for(int i = 0; i < nc ; ++i)
    {
      pokeColour(Qprop, peekColour(psi, i), i, c);
    }

}

//////////////////////////////////////////////////
//  Apply anti-periodic boundary conditions in time
////////////////////////////////////////////////////

void anti_periodic( LatticeGaugeField & Umu , int nt)
{
  int mu = 3 ;  // time directiom
  GridBase *grid = Umu._grid;
  // conformable(grid,g._grid);
  LatticeColourMatrix U(grid);

  U= PeekIndex<LorentzIndex>(Umu,mu);

  // code hacked from Grid/Grid/qcd/action/fermion/FermionOperatorImpl.h`
  Lattice<iScalar<vInteger> > t(grid); LatticeCoordinate(t,3);
  LatticeComplex phases(grid); phases=1.0;

  --nt ;
  phases = where( t ==(Integer)nt, phases,-phases);
  U = phases * U ;

  PokeIndex<LorentzIndex>(Umu,U,mu);

  cout << "Anti-periodic boundary conditions in time applied" << endl ; 
}


#include "hybrid_ops.h"

int main (int argc, char ** argv)
{
  typedef typename ImprovedStaggeredFermionR::FermionField FermionField; 
  typedef typename ImprovedStaggeredFermionR::ComplexField ComplexField; 
  typename ImprovedStaggeredFermionR::ImplParams params; 


  Grid_init(&argc,&argv);

  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();
  GridCartesian               Grid(latt_size,simd_layout,mpi_layout);
  GridRedBlackCartesian     RBGrid(&Grid);

  for(int mu=0;mu<Nd;mu++) 
    {
      cout << "Lattice[ " << mu << " ]= " << latt_size[mu] << endl ;
    }


  std::vector<int> seeds({1,2,3,4});
  GridParallelRNG          pRNG(&Grid);  pRNG.SeedFixedIntegers(seeds);

  FermionField result(&Grid); result=zero;
///////////////////////////////////////////////////////////////
//		Staggered Sign Functions
///////////////////////////////////////////////////////////////

  LatticeInteger coor[4] = {&Grid, &Grid, &Grid, &Grid}; 
  LatticeInteger n[4] = {&Grid, &Grid, &Grid, &Grid};
  LatticeComplex signs[4] =  {&Grid, &Grid, &Grid, &Grid}; 
  LatticeComplex One(&Grid), minusOne(&Grid); One = 1; minusOne = -1;

  for(int m=0; m<4; m++) {
    LatticeCoordinate(coor[m], m);  
  }

  n[0] = coor[0]; n[1] = n[0] + coor[1]; n[2] = n[1] + coor[2]; n[3] = n[2] + coor[3];
  
  for(int i=0; i<4; i++) {
    signs[i] = where((mod(n[i],2)== (Integer) 1), minusOne, One) ;
  } 

// Decided to split the phases up a bit, here are the standard staggered sign functions, corresponding
// to the gamma matrices. the phase from the inversion of M is signs[3].  
//////////////////////////////////////////////////////////////

  //SU3::HotConfiguration(pRNG,Umu); 
  //SU3::ColdConfiguration(Umu); // Umu = 1
   LatticeGaugeField Umu(&Grid); 


    FieldMetaData header;
     
     std::cout <<GridLogMessage<<"**************************************"<<std::endl;
     std::cout <<GridLogMessage<<"** Reading back ILDG conf    *********"<<std::endl;
     std::cout <<GridLogMessage<<"**************************************"<<std::endl;
     emptyUserRecord record;
     std::string file("/home/gray/grid_code/grid-staggered/examples/ckpoint_ildg.4000");
     
     IldgReader _IldgReader;
     _IldgReader.open(file);
     _IldgReader.readConfiguration(Umu,header);
     _IldgReader.close();
   

  int t_dir = 3;
  int nt = latt_size[t_dir];

  anti_periodic(Umu, nt);

  const int g_trans = 0;
  if( g_trans == 1)
    {
      LatticeColourMatrix   g(&Grid); // Gauge xform
      SU3::RandomGaugeTransform(pRNG,Umu,g); // Unit gauge
      cout << "Random Gauge Transform applied "  << endl ; 
    }
  else
    {
   cout << "NO Gauge Transform applied "  << endl ; 
    }

  // This uses the Milc conventions. See gridStaggInvert.cc in the MILC code.
  // Naive staggered action
  RealD u0=1.0;
  RealD c1= 2.0 ;
  RealD c2= 0.0 ;

  // Naik coefficients (See reference in Grid library)
  // RealD c1=9.0/8.0; RealD c2=-1.0/24.0;

  RealD mass=0.1;

  ImprovedStaggeredFermionR Ds(Umu,Umu,Grid,RBGrid,2.0*mass,c1,c2,u0);

  MdagMLinearOperator<ImprovedStaggeredFermionR,FermionField> HermOp(Ds);
  ConjugateGradient<FermionField> CG(1.0e-6,10000);
  //LinearOperatorBase<FermionField> LinOp(Ds);

  //  ./Grid/qcd/QCD.h
  LatticeStaggeredFermion local_src(&Grid) ;
  LatticeStaggeredFermion out(&Grid) ;
  LatticeStaggeredPropagator Qprop[2] = {&Grid, &Grid}  ;

  Qprop[1] = Qprop[0] = zero ;
  cout << "Inversion for mass "   << mass << endl ; 
  
  enum dirs {x=0, y=1, z=2, t=3};
//  dirs dir = x; // just a bit more readable

  // Compute the staggered quark propagators
  for(int k=0; k<2; k++) {

    for(int ic = 0 ; ic < 3 ; ++ic)
      {
        cout << "---------------------------------------------" << endl ; 
        cout << "Inversion for colour " << ic << endl ; 
        cout << "---------------------------------------------" << endl ; 

        // create point source
     
        std::vector<int> site({0,0,0,0});
        ColourVector cv = zero;
        cv()()(ic)=1.0;  
        local_src = zero;
        pokeSite(cv,local_src,site);
      
        if(k) local_src = onemp_local(Umu, local_src, signs, z); 
        Ds.Mdag(local_src, out) ; // apply Mdagger
        local_src = out;

        // invert 
        out = zero ;  // intial guess
        CG(HermOp, local_src, out);

	if(k) out = onemp_local(Umu, out, signs, z);
        // add solution to propagator structure
        FermToProp_s(Qprop[k], out , ic  ) ; 
      }
  }
  //  -----  Use the quark propagator to compute the correlator

 
  std::vector<TComplex> corr(nt) ;

  // contract the quark propagators
  LatticeComplex  c(&Grid);

  
  c = trace(adj(Qprop[0]) * Qprop[1]) ; 
  //c = c * signs[2];	// phase from inversion of M 
  
  //  this correlator over the lattice is summed over the spatial
  //   lattice at each timeslice t.
  cout << "\nTp = " << Tp  << "\n"; 
  sliceSum(c, corr, Tp);
  

  // output the correlators
  cout << "\n\n1-+ HYBRID meson \n\n";

//  cout << "\nCorrelator in spin component " << z << endl; 
  for(int tt = 0 ; tt < nt ; ++tt) {  
    cout << tt << " "  <<  real(corr[tt])  << "  " << imag(corr[tt]) << endl ;
  }
  
Grid_finalize();

}
