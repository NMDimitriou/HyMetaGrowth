#include "declarations.cuh"


void scaleGrowth(InitialData& id)
{
	//we have e^(st) the exponential phase
	// Ae^(st) = 2A => Ae^(st) = Ae^(ln2)
	// the douling time is: td = ln2/s
}


/*
#define N      14000    // number of data points to fit 
#define TMAX   (14.0) // time variable in [0,TMAX] 

struct data {
  size_t n;
  double * t ;
  double * y ;
  int y0;
  double s;
};

int expb2_f (const gsl_vector * x, void *data, gsl_vector * f)
{
  size_t n = ((struct data *)data)->n;
  double *t = ((struct data *)data)->t ;
  double *y = ((struct data *)data)->y ;

  int y0   = ((struct data *)data)->y0 ;
  double s = ((struct data *)data)->s  ;

  double sp = gsl_vector_get (x, 0);

  size_t i;

  for (i = 0; i < n; i++)
    {
      // Model Yi = 2^ (sp * s * ti) //
      double Yi = y0 * pow( 2.0f, sp * s * t[i] );
      gsl_vector_set (f, i, Yi - y[i]);
    }

  return GSL_SUCCESS;
}


int expb2_df (const gsl_vector * x, void *data, gsl_matrix * J)
{
  size_t n = ((struct data *)data)->n;
  double *t = ((struct data *)data)->t;

  int y0   = ((struct data *)data)->y0 ;
  double s = ((struct data *)data)->s  ;

  double sp = gsl_vector_get (x, 0);

  size_t i;

  for (i = 0; i < n; i++)
    {
      // Jacobian matrix J(i,j) = dfi / dxj, /
      // where fi = (Yi - yi)/sigma[i],     /
      //       Yi = Yi = 2^ (sp * s * ti)  /
      // and the xj is the parameter (sp) /
      double df = y0*s*t[i]*log(2.0f)*pow(2.0f, sp * s * t[i]);
      gsl_matrix_set (J, i, 0, df);
    }

  return GSL_SUCCESS;
}

void callback(const size_t iter, void *params,
         const gsl_multifit_nlinear_workspace *w)
{
  gsl_vector *f = gsl_multifit_nlinear_residual(w);
  gsl_vector *x = gsl_multifit_nlinear_position(w);
  double rcond;

  // compute reciprocal condition number of J(x) /
  gsl_multifit_nlinear_rcond(&rcond, w);

  fprintf(stderr, "iter %2zu: scale = %.4f, cond(J) = %8.4f, |f(x)| = %.4f\n",
          iter,
          gsl_vector_get(x, 0),
          1.0 / rcond,
          gsl_blas_dnrm2(f));
}



void scaleGrowth(DataArray& host, InitialData& id)
{
  const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
  gsl_multifit_nlinear_workspace *w;
  gsl_multifit_nlinear_fdf fdf;
  gsl_multifit_nlinear_parameters fdf_params =
    gsl_multifit_nlinear_default_parameters();
  const size_t n = N;
  const size_t p = 1;

  gsl_vector *f;
  gsl_matrix *J;
  gsl_matrix *covar = gsl_matrix_alloc (p, p);
  double t[N], y[N], weights[N];
  int y0=0;
  //count cells;
  char filename[4096];
  sprintf(filename, "IC/coord_idx_%s_D0.txt", host.dat);
  int nsp; //cell spatial index
  ifstream read(filename);
  while(read>>nsp)
  {
  	y0 ++;
  }
  double s = id.s;
  struct data d = { n, t, y, y0, s};
  double x_init[1] = { 1.0 }; // starting values /
  gsl_vector_view x = gsl_vector_view_array (x_init, p);
  gsl_vector_view wts = gsl_vector_view_array(weights, n);
  gsl_rng * r;
  double chisq, chisq0;
  int status, info;
  size_t i;

  const double xtol = 1e-8;
  const double gtol = 1e-8;
  const double ftol = 0.0;

  gsl_rng_env_setup();
  r = gsl_rng_alloc(gsl_rng_default);

  // define the function to be minimized /
  fdf.f  = expb2_f;
  fdf.df = expb2_df;  // set to NULL for finite-difference Jacobian /
  fdf.fvv = NULL;     // not using geodesic acceleration /
  fdf.n = n;
  fdf.p = p;
  fdf.params = &d;


  // this is the data to be fitted /
  for (i = 0; i < n; i++)
    {
      double ti = i * TMAX / (n - 1.0);
      double yi = y0* exp (id.s * ti);
      //double si = 0.1 * yi;
      //double dy = gsl_ran_gaussian(r, si);

      t[i] = ti;
      y[i] = yi;
      weights[i] = 1.0;// / (si * si);
      if(i==0) {   printf ("cells: %g at t = %g\n", y[i], ti);}
    };

  // allocate workspace with default parameters /
  w = gsl_multifit_nlinear_alloc (T, &fdf_params, n, p);

  // initialize solver with starting point and weights /
  gsl_multifit_nlinear_winit (&x.vector, &wts.vector, &fdf, w);

  // compute initial cost function /
  f = gsl_multifit_nlinear_residual(w);
  gsl_blas_ddot(f, f, &chisq0);

  printf("--------------------------------------------------- \n");
  printf("Initializing estimation of growth scale coefficient \n");
  // solve the system with a maximum of 100 iterations /
  status = gsl_multifit_nlinear_driver(100, xtol, gtol, ftol,
                                       callback, NULL, &info, w);

  // compute covariance of best fit parameters /
  J = gsl_multifit_nlinear_jac(w);
  gsl_multifit_nlinear_covar (J, 0.0, covar);

  // compute final cost /
  gsl_blas_ddot(f, f, &chisq);

#define FIT(i) gsl_vector_get(w->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

  fprintf(stderr, "summary from method '%s/%s'\n",
          gsl_multifit_nlinear_name(w),
          gsl_multifit_nlinear_trs_name(w));
  fprintf(stderr, "number of iterations: %zu\n",
          gsl_multifit_nlinear_niter(w));
  fprintf(stderr, "function evaluations: %zu\n", fdf.nevalf);
  fprintf(stderr, "Jacobian evaluations: %zu\n", fdf.nevaldf);
  fprintf(stderr, "reason for stopping: %s\n",
          (info == 1) ? "small step size" : "small gradient");
  fprintf(stderr, "initial |f(x)| = %f\n", sqrt(chisq0));
  fprintf(stderr, "final   |f(x)| = %f\n", sqrt(chisq));

  {
    double dof = n - p;
    double c = GSL_MAX_DBL(1, sqrt(chisq / dof));

    fprintf(stderr, "chisq/dof = %g\n", chisq / dof);

    fprintf (stderr, "scale coeff = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
  }

  fprintf (stderr, "status = %s\n", gsl_strerror (status));
  id.scale = FIT(0);

  gsl_multifit_nlinear_free (w);
  gsl_matrix_free (covar);
  gsl_rng_free (r);

  
  printf("\nGrowth scale coefficient was found to be: %f\n", id.scale);
  printf("--------------------------------------------------- \n");
  

}
*/
