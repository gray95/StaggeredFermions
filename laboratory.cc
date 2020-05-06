/***************************************************************************************************


    Perhaps a better starting point:
    typedef typename ImprovedStaggeredFermionR::FermionField FermionField; 


    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_wilson_cg_unprec.cc
    Starting point:    Test_staggered_cg_unprec.cc 

	THIS IS A SCRIPT TO TEST / TRY OUT VARIOUS BITS OF GRID IN A (hopefully) CONTROLLED MANNER
  

****************************************************************************************************/
#include <Grid/Grid.h>
#include <Grid/algorithms/iterative/BlockConjugateGradient.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

template<class d>
struct scal {
  d internal;
};

Gamma::Algebra Gmu [] = {
  Gamma::Algebra::GammaX,
  Gamma::Algebra::GammaY,
  Gamma::Algebra::GammaZ,
  Gamma::Algebra::GammaT
  };


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


  std::vector<int> seeds({1,2,3,4});
  GridParallelRNG          pRNG(&Grid);  pRNG.SeedFixedIntegers(seeds);

  FermionField result(&Grid); result=zero;

  //  create a hot su3 configuration
  LatticeGaugeField Umu(&Grid); 
  LatticeRealD  a(&Grid), b(&Grid), res(&Grid); 
  LatticeInteger coor[3] = {&Grid, &Grid, &Grid};
  LatticeInteger test(&Grid); 
  LatticeCoordinate(coor[0], 0); LatticeCoordinate(coor[1], 1); LatticeCoordinate(coor[2], 2);
  test = coor[0] + coor[1] + coor[2] ;  
  test = Cshift(test, 0, 1);
  cout << test; 
// End of the Grid
  Grid_finalize();
}
