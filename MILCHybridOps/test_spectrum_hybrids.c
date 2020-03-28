/* Code to implement the lightest hybrid multiplet in MILC
 *
 * 	{ 1mp, 0mp, 1mm, 2mp } - staggered, local ops
 *
 * 	void mult_(J^PC)(taste)_field();
 *
 */


/* Craig wrote this - (renamed from 'mult_1mp_lrho_field') */
void mult_1mpi_field( int pdir, su3_vector *src_in, su3_vector *dest );

// declared & defined in this file
void mult_1mm5_field(int pdir,  su3_vector *src_in, su3_vector *dest) ;
void mult_0mpi_field( su3_vector *src_in, su3_vector *dest ) ;
void mult_2mpi_field(int pdir, su3_vector *src_in, su3_vector *dest) ;



//##############################################################################
/* "Multiply by" a one-minus-minus hybrid operator 
    quark operator is "pion5"(ie do nothing), gluon operator is magnetic field.
    1--: \psi-bar gamma_5 B_i \psi
   "pdir" is the polarization direction of the meson */
void mult_1mm5_field( int pdir, su3_vector *src_in, su3_vector *dest ){

  char debug[][4] = { "XUP", "YUP", "ZUP", "TUP" };

  register int dir,i;
  register site *s;

  static int do_init = 1 ;
  if( do_init )
    {
	int dummy = init_hybrids() ;
	do_init = 0 ;
    }

  su3_vector *cg_p = create_v_field();
  su3_vector *src = create_v_field();

  /** copy source so src_in and dest can be the same ***/
  copy_v_field(src, src_in);

  /* set destination to zero */
  FORALLSITES(i,s){ clearvec( dest + i ); }

  mult_by_field_strength( (pdir+1)%3, (pdir+2)%3, src, cg_p );
  
  FORALLSITES(i,s){
    add_su3_vector( dest + i,  cg_p + i , dest + i ) ;
  }

  destroy_v_field(cg_p);
  destroy_v_field(src);

} /* end mult_1mm5_field */



//##############################################################################
/* "Multiply by" a one-minus-plus operator 
    quark operator is "rho0", gluon operator is magnetic field B.
    Cross product is spin one.
    Z component = "rho_X*F_XZ + rho_Y*F_YZ.
   "pdir" is the polarization direction of the meson */
void mult_1mpi_field( int pdir, su3_vector *src_in, su3_vector *dest )
{
  char debug[][4] = {"XUP", "YUP", "ZUP", "TUP"  };

  register int dir,i;
  register site *s;
  int dir_a , dir_b ; 

  static int do_init = 1 ;

  if( do_init ){
    int dummy = init_hybrids() ;
    do_init = 0 ;
  }

  switch (pdir) {
    case XUP:
      dir_a = YUP; dir_b = ZUP;
      break;
    case YUP:
      dir_a = ZUP; dir_b = XUP;
      break;
    case ZUP:
      dir_a = XUP; dir_b = YUP;
      break;
    default:
      if(this_node==0)printf("pdir = %d not coded up \n" , pdir);
      terminate(0);
  }

  /* use cg_p and ttt as temporary storage */
  su3_vector *ttt = create_v_field();
  su3_vector *cg_p = create_v_field();
  su3_vector *src = create_v_field();

  int r0[4] = { 0 , 0 ,0 , 0 }  ;

  /** copy source so src_in and dest can be the same ***/
  copy_v_field(src, src_in);

  /* set destination to zero */
  FORALLSITES(i,s){ clearvec( dest + i ); }

  mult_by_field_strength( pdir, dir_a, src, cg_p );
  mult_rhoi_field(dir_a, r0, cg_p, ttt ) ;

  FORALLSITES(i,s){
    add_su3_vector( dest + i,  ttt + i , dest + i ) ;
  }

  mult_by_field_strength( pdir, dir_b, src, cg_p );
  mult_rhoi_field(dir_b, r0, cg_p, ttt ) ;

  FORALLSITES(i,s){
    sub_su3_vector( dest + i,  ttt + i , dest + i ) ;
  }

  destroy_v_field(ttt);
  destroy_v_field(cg_p);
  destroy_v_field(src);

} /* end mult_1mpi_field */



//##########################################################################
/* "Multiply by" a zero-minus-plus hybrid operator 
    quark operator is "rho", gluon operator is magnetic field.
    0-+ = \psibar gamma_i \psi B_i */
void mult_0mpi_field( su3_vector *src, su3_vector *dest ){

  char debug[][4] = { "XUP", "YUP", "ZUP", "TUP" } ;

  register int dir,i;
  register site *s;

  static int do_init = 1 ;
  if( do_init ){
    int dummy = init_hybrids() ;
	do_init = 0 ;
  }

  /* use cg_p and ttt as temporary storage */
  su3_vector *ttt = create_v_field();
  su3_vector *cg_p = create_v_field();
  su3_vector *src = create_v_field();

  int r0[4] = { 0 , 0 ,0 , 0 }  ;

  /* copy source so src_in and dest can be the same */
  copy_v_field(src, src_in);

  /* set destination to zero */
  FORALLSITES(i,s){
    clearvec( dest + i );
  }

  /* Loop over spatial directions */
  for(dir=XUP; dir<=ZUP ; dir++){
         mult_by_field_strength( (dir+1)%3, (dir+2)%3, src, cg_p ); 
	 mult_rhoi_field(dir, r0, cg_p, ttt ) ;

	 FORALLSITES(i,s){
	  add_su3_vector( dest + i,  ttt + i , dest + i ) ;
	}
  }

  destroy_v_field(ttt);
  destroy_v_field(cg_p);
  destroy_v_field(src);

} /* end mult_0mpi_field */



//########################################################################
/* "Multiply by" a two-minus-plus hybrid operator. 
    Quark operator is "rhoi", gluon operator is magnetic field.
    2-+ = abs(epsilon_ijk) \psibar gamma_j \psi B_k */
void mult_2mpi_field( int pdir, su3_vector *src_in, su3_vector *dest ){

  char debug[][4] = {"XUP", "YUP", "ZUP", "TUP"  };

  register int dir,i;
  register site *s;
  int dir_a , dir_b ; 

  static int do_init = 1 ;

  if( do_init )
    {
      int dummy = init_hybrids() ;
      do_init = 0 ;
    }

  switch (pdir) {
    case XUP:
      dir_a = YUP; dir_b = ZUP;
      break;
    case YUP:
      dir_a = ZUP; dir_b = XUP;
      break;
    case ZUP:
      dir_a = XUP; dir_b = YUP;
      break;
    default:
      if(this_node==0)printf("pdir = %d not coded up \n" , pdir);
      terminate(0);
  }


  /* use cg_p and ttt as temporary storage */
  su3_vector *ttt = create_v_field();
  su3_vector *cg_p = create_v_field();
  su3_vector *src = create_v_field();

  int r0[4] = { 0 , 0 ,0 , 0 }  ;

  /** copy source so src_in and dest can be the same ***/
  copy_v_field(src, src_in);

  /* set destination to zero */
  FORALLSITES(i,s){
    clearvec( dest + i );
  }

  mult_by_field_strength( dir_a, pdir, src, cg_p );
  mult_rhoi_field(dir_a, r0, cg_p, ttt ) ;
  FORALLSITES(i,s){
    add_su3_vector( dest + i,  ttt + i , dest + i ) ;
  }

  mult_by_field_strength( dir_b, pdir, src, cg_p );
  mult_rhoi_field(dir_b, r0, cg_p, ttt ) ;
  FORALLSITES(i,s){
    add_su3_vector( dest + i,  ttt + i , dest + i ) ;
  }

  destroy_v_field(ttt);
  destroy_v_field(cg_p);
  destroy_v_field(src);

} /* end mult_2mpi_field */


#if 0
  node0_printf("DEBUG tup %d %s\n", TUP, debug[TUP] );
  node0_printf("DEBUG xup %d %s\n", XUP, debug[XUP] );
  node0_printf("DEBUG yup %d %s\n", YUP, debug[YUP] );
  node0_printf("DEBUG zup %d %s\n", ZUP, debug[ZUP] );

  node0_printf("DEBUG pdir = %s\n", debug[pdir] );
#endif


