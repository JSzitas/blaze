


get_hessian <- function( fun, pars, parscale = rep(1, length(pars)), eps = 1e-03 ) {

  diff_par = rep(NA,length(pars))

  for( i in  seq_along(pars) ) {
    eps = eps/parscale[i]
    # get differences for 1 eps positive
    df_1 <- f_min_grad( dpar )



  }

  for (i = 0; i < npar; i++) {
    eps = OS->ndeps[i]/(OS->parscale[i]);
    dpar[i] = dpar[i] + eps;
    fmingr(npar, dpar, df1, (void *)OS);
    dpar[i] = dpar[i] - 2 * eps;
    fmingr(npar, dpar, df2, (void *)OS);
    for (j = 0; j < npar; j++)
      REAL(ans)[i * npar + j] = (OS->fnscale) * (df1[j] - df2[j])/
        (2 * eps * (OS->parscale[i]) * (OS->parscale[j]));
    dpar[i] = dpar[i] + eps;
  }
  // now symmetrize
  for (i = 0; i < npar; i++)
    for (j = 0; j < i; j++) {
      double tmp =
        0.5 * (REAL(ans)[i * npar + j] + REAL(ans)[j * npar + i]);
      REAL(ans)[i * npar + j] = REAL(ans)[j * npar + i] = tmp;
    }
}
