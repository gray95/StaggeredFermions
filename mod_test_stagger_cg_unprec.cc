    /*************************************************************************************


    Perhaps a better starting point:
    typedef typename ImprovedStaggeredFermionR::FermionField FermionField; 


    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_wilson_cg_unprec.cc
    Starting point:    Test_staggered_cg_unprec.cc 


    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
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

void FermToProp_s(LatticeStaggeredPropagator Qprop, LatticeStaggeredFermion psi , const int c)
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

  FermionField src(&Grid); random(pRNG,src);
  RealD nrm = norm2(src);
  FermionField result(&Grid); result=zero;

  //  create a hot su3 configuration
  LatticeGaugeField Umu(&Grid); 
  SU3::HotConfiguration(pRNG,Umu);

  double volume=1;
  for(int mu=0;mu<Nd;mu++){
    volume=volume*latt_size[mu];
  }  
  
  RealD mass=0.1;		// Are these the right numerical values for HISQ action?
  RealD c1=9.0/8.0;
  RealD c2=-1.0/24.0;
  RealD u0=1.0;
  ImprovedStaggeredFermionR Ds(Umu, Umu, Grid, RBGrid, mass, c1, c2, u0, params);

  MdagMLinearOperator<ImprovedStaggeredFermionR,FermionField> HermOp(Ds);
  ConjugateGradient<FermionField> CG(1.0e-8,10000);


  //  ./Grid/qcd/QCD.h
  LatticeStaggeredFermion local_src(&Grid) ;
  LatticeStaggeredFermion out(&Grid) ;
  LatticeStaggeredPropagator Qprop(&Grid) ;

  

  // Compute the staggered quark propagator
  for(int ic = 0 ; ic < 3 ; ++ic)
    {

      // create point source

      // tests/core/Test_staggered5Dvec.cc
      // src is FermionField
      std::vector<int> site({0,0,0,0});	// 5D or 4D ?
      ColourVector cv = zero;
      cv()()(ic)=1.0;  //  add ic craig
      src = zero;
      pokeSite(cv,local_src,site);

      // invert 
      out = zero ;  // intial guess

      CG(HermOp,local_src,out);

      // add solution to propagator structure
      FermToProp_s(Qprop, out , ic  ) ; 

    }

  // 
  //  -----  Use the quark propagator to compute the pion correlator
  //

  int t_dir = 3;
  int nt =latt_size[t_dir];

  // pion correlator
  std::vector<TComplex> corr(nt);

  // contract the quark propagators
  LatticeComplex  c(&Grid) ;
  c = trace(Qprop * adj(Qprop)) ;

  //  this correlator over the lattice is summed over the spatial
  //   lattice at each timeslice t.
  sliceSum(c, corr, Tp);

  // output the correlators
  for(int tt = 0 ; tt < nt ; ++tt)
    {
      cout << tt << " "  << real(corr[tt])  << endl ;
    }

  // End of the Grid
  Grid_finalize();
}

    /***************************************************************************
     * Convert a propagator to a fermion & back.
     **************************************************************************/
/*    LatticeFermion ferm(&UGrid);
    LatticePropagator prop(&UGrid), ref(&UGrid);
    gaussian(rng4, prop);

    // Define variables for sanity checking a single site.
    typename SpinColourVector::scalar_object fermSite;
    typename SpinColourMatrix::scalar_object propSite;
    std::vector<int> site(Nd, 0);

    for (int s = 0; s < Ns; ++s)
    for (int c = 0; c < Nc; ++c)
    {
        ref = prop;
        PropToFerm<WilsonImplR>(ferm, prop, s, c);
        FermToProp<WilsonImplR>(prop, ferm, s, c);

        std::cout << "Spin = " << s << ", Colour = " << c << std::endl;
        ref -= prop;
        std::cout << "Prop->Ferm->Prop test, diff = " << Grid::sqrt(norm2(ref)) << std::endl;

        peekSite(fermSite, ferm, site);
        peekSite(propSite, prop, site);
        for (int s2 = 0; s2 < Ns; ++s2)
        for (int c2 = 0; c2 < Nc; ++c2)
        {
            if (propSite()(s2, s)(c2, c) != fermSite()(s2)(c2))
            {
                std::cout << propSite()(s2, s)(c2, c) << " != "
                          << fermSite()(s2)(c2) << " for spin = " << s2
                          << ", col = " << c2 << std::endl;
            }
        }
    }

*/
