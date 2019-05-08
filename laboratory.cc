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
  LatticeRealD test(&Grid), a(&Grid), b(&Grid); test = zero;
  LatticeInteger coor(&Grid);  //  SU3::HotConfiguration(pRNG,Umu);
  LatticeCoordinate(coor, 1);
  SU3::ColdConfiguration(Umu); // Umu = 1  
  a = -1; b = 1;
  test = where((mod(coor,2)==1), a, b);
  cout << test; 
/*
  const int g_trans =  0;
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


  //  cout << "DEBUG " << Umu << endl ; 

  double volume=1;
  for(int mu=0;mu<Nd;mu++){
    volume=volume*latt_size[mu];
  }  
  
  RealD mass=0.1;
  RealD c1=9.0/8.0;
  RealD c2=-1.0/24.0;
  RealD u0=1.0;
  ImprovedStaggeredFermionR Ds(Umu,Umu,Grid,RBGrid,mass,c1,c2,u0);

  MdagMLinearOperator<ImprovedStaggeredFermionR,FermionField> HermOp(Ds);
  ConjugateGradient<FermionField> CG(1.0e-6,10000);
  //LinearOperatorBase<FermionField> LinOp(Ds);

  //  ./Grid/qcd/QCD.h
  LatticeStaggeredFermion local_src(&Grid) ;
  LatticeStaggeredFermion out(&Grid) ;
  LatticeStaggeredPropagator Qprop(&Grid)  ;

  Qprop = zero ;

  cout << "Inversion for mass "   << mass << endl ; 

  // Compute the staggered quark propator
  for(int ic = 0 ; ic < 3 ; ++ic)
    {
      cout << "---------------------------------------------" << endl ; 
      cout << "Inversion for colour " << ic << endl ; 
      cout << "---------------------------------------------" << endl ; 

      // create point source

      // tests/core/Test_staggered5Dvec.cc
      // src is FermiionFiel
      std::vector<int> site({0,0,0,0});
      ColourVector cv = zero;
      cv()()(ic)=1.0;  //  add ic craig
      local_src = zero;
      pokeSite(cv,local_src,site);

      // invert 
       out = zero ;  // intial guess

      CG(HermOp,local_src,out);
      //  cout << "out = " << out << endl ; 
      // Apply M^\dagger
      Ds.Mdag(out, local_src) ;
      //Ds.AdjOp(out, local_src);
      out = local_src ;

      //cout << "Colour Vector: " << cv << endl;
      //cout << "Site: " << site << endl;
      //cout << "out: " << out << endl;

      // add solution to propagator structure
      FermToProp_s(Qprop, out , ic  ) ; 
    }
* *************************************************************
 *
 * 			DIAGNOSTICS/DEBUGGING
 *
 **************************************************************

  //  -----  Use the quark propagator to compute the pion correlator
  //
  // cout << "Qprop = " << Qprop << endl ; 

  int t_dir = 3;
  int nt =latt_size[t_dir];

  // pion correlator
  std::vector<TComplex> corr(nt)  ;

  // contract the quark propagators
  LatticeComplex  c(&Grid)  ;

  c = trace(Qprop * adj(Qprop)) ; 

  //  this correlator over the lattice is summed over the spatial
  //   lattice at each timeslice t.
  cout << "Tp = " << Tp  << endl; 
  sliceSum(c, corr, Tp);

  // output the correlators
#if 1
  for(int tt = 0 ; tt < nt ; ++tt)
    {
      double ttt = real(corr[tt]) ;
      cout << tt << " "  <<  ttt  << endl ;
    }
#endif
  //  cout << "corr = " << corr << endl ; 
*/
  // End of the Grid
  Grid_finalize();
}
