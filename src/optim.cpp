
SEXP getListElement(SEXP list, char *str)
{
  if (!isNewList(list))
    error(_("invalid argument type"));
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  int i;

  for (i = 0; i < length(list); i++)
    if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
    return elmt;
}

static double * vect(int n)
{
  return (double *)R_alloc(n, sizeof(double));
}

typedef struct opt_struct
{
  SEXP R_fcall;    /* function */
  SEXP R_gcall;    /* gradient */
  SEXP R_env;      /* where to evaluate the calls */
  double* ndeps;   /* tolerances for numerical derivatives */
  double fnscale;  /* scaling for objective */
  double* parscale;/* scaling for parameters */
  int usebounds;
  double* lower, *upper;
  SEXP names;	     /* names for par */
} opt_struct, *OptStruct;



static double fminfn(int n, double *p, void *ex)
{
  SEXP s, x;
  int i;
  double val;
  OptStruct OS = (OptStruct) ex;
  PROTECT_INDEX ipx;

  PROTECT(x = allocVector(REALSXP, n));
  if(!isNull(OS->names)) setAttrib(x, R_NamesSymbol, OS->names);
  for (i = 0; i < n; i++) {
    if (!R_FINITE(p[i])) error(_("non-finite value supplied by optim"));
    REAL(x)[i] = p[i] * (OS->parscale[i]);
  }
  SETCADR(OS->R_fcall, x);
  PROTECT_WITH_INDEX(s = eval(OS->R_fcall, OS->R_env), &ipx);
  REPROTECT(s = coerceVector(s, REALSXP), ipx);
  if (LENGTH(s) != 1)
    error(_("objective function in optim evaluates to length %d not 1"),
          LENGTH(s));
  val = REAL(s)[0]/(OS->fnscale);
  UNPROTECT(2);
  return val;
}


void
vmmin(int n0, double *b, double *Fmin, optimfn fminfn, optimgr fmingr,
      int maxit, int trace, int *mask,
      double abstol, double reltol, int nREPORT, void *ex,
      int *fncount, int *grcount, int *fail)
{
  Rboolean accpoint, enough;
  double *g, *t, *X, *c, **B;
  int   count, funcount, gradcount;
  double f, gradproj;
  int   i, j, ilast, iter = 0;
  double s, steplength;
  double D1, D2;
  int   n, *l;

  if (maxit <= 0) {
    *fail = 0;
    *Fmin = fminfn(n0, b, ex);
    *fncount = *grcount = 0;
    return;
  }

  if (nREPORT <= 0)
    error(_("REPORT must be > 0 (method = \"BFGS\")"));
  l = (int *) R_alloc(n0, sizeof(int));
  n = 0;
  for (i = 0; i < n0; i++) if (mask[i]) l[n++] = i;
  g = vect(n0);
  t = vect(n);
  X = vect(n);
  c = vect(n);
  B = Lmatrix(n);
  f = fminfn(n0, b, ex);
  if (!R_FINITE(f))
    error(_("initial value in 'vmmin' is not finite"));
  if (trace) Rprintf("initial  value %f \n", f);
  *Fmin = f;
  funcount = gradcount = 1;
  fmingr(n0, b, g, ex);
  iter++;
  ilast = gradcount;

  do {
    if (ilast == gradcount) {
      for (i = 0; i < n; i++) {
        for (j = 0; j < i; j++) B[i][j] = 0.0;
        B[i][i] = 1.0;
      }
    }
    for (i = 0; i < n; i++) {
      X[i] = b[l[i]];
      c[i] = g[l[i]];
    }
    gradproj = 0.0;
    for (i = 0; i < n; i++) {
      s = 0.0;
      for (j = 0; j <= i; j++) s -= B[i][j] * g[l[j]];
      for (j = i + 1; j < n; j++) s -= B[j][i] * g[l[j]];
      t[i] = s;
      gradproj += s * g[l[i]];
    }

    if (gradproj < 0.0) {	/* search direction is downhill */
      steplength = 1.0;
      accpoint = FALSE;
      do {
        count = 0;
        for (i = 0; i < n; i++) {
          b[l[i]] = X[i] + steplength * t[i];
          if (reltest + X[i] == reltest + b[l[i]]) /* no change */
            count++;
        }
        if (count < n) {
          f = fminfn(n0, b, ex);
          funcount++;
          accpoint = R_FINITE(f) &&
            (f <= *Fmin + gradproj * steplength * acctol);
          if (!accpoint) {
            steplength *= stepredn;
          }
        }
      } while (!(count == n || accpoint));
      enough = (f > abstol) &&
        fabs(f - *Fmin) > reltol * (fabs(*Fmin) + reltol);
      /* stop if value if small or if relative change is low */
      if (!enough) {
        count = n;
        *Fmin = f;
      }
      if (count < n) {/* making progress */
      *Fmin = f;
        fmingr(n0, b, g, ex);
        gradcount++;
        iter++;
        D1 = 0.0;
        for (i = 0; i < n; i++) {
          t[i] = steplength * t[i];
          c[i] = g[l[i]] - c[i];
          D1 += t[i] * c[i];
        }
        if (D1 > 0) {
          D2 = 0.0;
          for (i = 0; i < n; i++) {
            s = 0.0;
            for (j = 0; j <= i; j++)
              s += B[i][j] * c[j];
            for (j = i + 1; j < n; j++)
              s += B[j][i] * c[j];
            X[i] = s;
            D2 += s * c[i];
          }
          D2 = 1.0 + D2 / D1;
          for (i = 0; i < n; i++) {
            for (j = 0; j <= i; j++)
              B[i][j] += (D2 * t[i] * t[j]
                            - X[i] * t[j] - t[i] * X[j]) / D1;
          }
        } else {	/* D1 < 0 */
      ilast = gradcount;
        }
      } else {	/* no progress */
      if (ilast < gradcount) {
        count = 0;
        ilast = gradcount;
      }
      }
    } else {		/* uphill search */
      count = 0;
      if (ilast == gradcount) count = n;
      else ilast = gradcount;
      /* Resets unless has just been reset */
    }
    if (trace && (iter % nREPORT == 0))
      Rprintf("iter%4d value %f\n", iter, f);
    if (iter >= maxit) break;
    if (gradcount - ilast > 2 * n)
      ilast = gradcount;	/* periodic restart */
  } while (count != n || ilast != gradcount);
  if (trace) {
    Rprintf("final  value %f \n", *Fmin);
    if (iter < maxit) Rprintf("converged\n");
    else Rprintf("stopped after %i iterations\n", iter);
  }
  *fail = (iter < maxit) ? 0 : 1;
  *fncount = funcount;
  *grcount = gradcount;
}


/* par fn gr method options */
SEXP optim(SEXP call, SEXP op, SEXP args, SEXP rho)
{
  SEXP par, fn, gr, method, options, tmp, slower, supper;
  SEXP res, value, counts, conv;
  int i, npar=0, *mask, trace, maxit, fncount = 0, grcount = 0, nREPORT, tmax;
  int ifail = 0;
  double *dpar, *opar, val = 0.0, abstol, reltol, temp;
  const char *tn;
  OptStruct OS;
  PROTECT_INDEX par_index;

  args = CDR(args);
  OS = (OptStruct) R_alloc(1, sizeof(opt_struct));
  OS->usebounds = 0;
  OS->R_env = rho;
  par = CAR(args);
  OS->names = getAttrib(par, R_NamesSymbol);
  args = CDR(args); fn = CAR(args);
  if (!isFunction(fn)) error(_("'fn' is not a function"));
  args = CDR(args); gr = CAR(args);
  args = CDR(args); method = CAR(args);
  if (!isString(method)|| LENGTH(method) != 1)
    error(_("invalid '%s' argument"), "method");
  tn = CHAR(STRING_ELT(method, 0));
  args = CDR(args); options = CAR(args);
  PROTECT(OS->R_fcall = lang2(fn, R_NilValue));
  PROTECT_WITH_INDEX(par = coerceVector(par, REALSXP), &par_index);
  if (MAYBE_REFERENCED(par))
    REPROTECT(par = duplicate(par), par_index);
  npar = LENGTH(par);
  dpar = vect(npar);
  opar = vect(npar);
  trace = asInteger(getListElement(options, "trace"));
  OS->fnscale = asReal(getListElement(options, "fnscale"));
  tmp = getListElement(options, "parscale");
  if (LENGTH(tmp) != npar)
    error(_("'parscale' is of the wrong length"));
  PROTECT(tmp = coerceVector(tmp, REALSXP));
  OS->parscale = vect(npar);
  for (i = 0; i < npar; i++) OS->parscale[i] = REAL(tmp)[i];
  UNPROTECT(1);
  for (i = 0; i < npar; i++)
    dpar[i] = REAL(par)[i] / (OS->parscale[i]);
  PROTECT(res = allocVector(VECSXP, 5));
  SEXP names;
  PROTECT(names = allocVector(STRSXP, 5));
  SET_STRING_ELT(names, 0, mkChar("par"));
  SET_STRING_ELT(names, 1, mkChar("value"));
  SET_STRING_ELT(names, 2, mkChar("counts"));
  SET_STRING_ELT(names, 3, mkChar("convergence"));
  SET_STRING_ELT(names, 4, mkChar("message"));
  setAttrib(res, R_NamesSymbol, names);
  UNPROTECT(1);
  PROTECT(value = allocVector(REALSXP, 1));
  PROTECT(counts = allocVector(INTSXP, 2));
  SEXP countnames;
  PROTECT(countnames = allocVector(STRSXP, 2));
  SET_STRING_ELT(countnames, 0, mkChar("function"));
  SET_STRING_ELT(countnames, 1, mkChar("gradient"));
  setAttrib(counts, R_NamesSymbol, countnames);
  UNPROTECT(1);
  PROTECT(conv = allocVector(INTSXP, 1));
  abstol = asReal(getListElement(options, "abstol"));
  reltol = asReal(getListElement(options, "reltol"));
  maxit = asInteger(getListElement(options, "maxit"));
  if (maxit == NA_INTEGER) error(_("'maxit' is not an integer"));

  if (strcmp(tn, "Nelder-Mead") == 0) {
    double alpha, beta, gamm;

    alpha = asReal(getListElement(options, "alpha"));
    beta = asReal(getListElement(options, "beta"));
    gamm = asReal(getListElement(options, "gamma"));
    nmmin(npar, dpar, opar, &val, fminfn, &ifail, abstol, reltol,
          (void *)OS, alpha, beta, gamm, trace, &fncount, maxit);
    for (i = 0; i < npar; i++)
      REAL(par)[i] = opar[i] * (OS->parscale[i]);
    grcount = NA_INTEGER;

  }
  else if (strcmp(tn, "SANN") == 0) {
    tmax = asInteger(getListElement(options, "tmax"));
    temp = asReal(getListElement(options, "temp"));
    if (trace) trace = asInteger(getListElement(options, "REPORT"));
    if (tmax == NA_INTEGER || tmax < 1) // PR#15194
      error(_("'tmax' is not a positive integer"));
    if (!isNull(gr)) {
      if (!isFunction(gr)) error(_("'gr' is not a function"));
      PROTECT(OS->R_gcall = lang2(gr, R_NilValue));
    } else {
      PROTECT(OS->R_gcall = R_NilValue); /* for balance */
    }
    samin (npar, dpar, &val, fminfn, maxit, tmax, temp, trace, (void *)OS);
    for (i = 0; i < npar; i++)
      REAL(par)[i] = dpar[i] * (OS->parscale[i]);
    fncount = npar > 0 ? maxit : 1;
    grcount = NA_INTEGER;
    UNPROTECT(1);  /* OS->R_gcall */

  } else if (strcmp(tn, "BFGS") == 0) {
    SEXP ndeps;

    nREPORT = asInteger(getListElement(options, "REPORT"));
    if (!isNull(gr)) {
      if (!isFunction(gr)) error(_("'gr' is not a function"));
      PROTECT(OS->R_gcall = lang2(gr, R_NilValue));
    } else {
      PROTECT(OS->R_gcall = R_NilValue); /* for balance */
ndeps = getListElement(options, "ndeps");
if (LENGTH(ndeps) != npar)
  error(_("'ndeps' is of the wrong length"));
OS->ndeps = vect(npar);
PROTECT(ndeps = coerceVector(ndeps, REALSXP));
for (i = 0; i < npar; i++) OS->ndeps[i] = REAL(ndeps)[i];
UNPROTECT(1);
    }
    mask = (int *) R_alloc(npar, sizeof(int));
    for (i = 0; i < npar; i++) mask[i] = 1;
    vmmin(npar, dpar, &val, fminfn, fmingr, maxit, trace, mask, abstol,
          reltol, nREPORT, (void *)OS, &fncount, &grcount, &ifail);
    for (i = 0; i < npar; i++)
      REAL(par)[i] = dpar[i] * (OS->parscale[i]);
    UNPROTECT(1); /* OS->R_gcall */
  } else if (strcmp(tn, "CG") == 0) {
    int type;
    SEXP ndeps;

    type = asInteger(getListElement(options, "type"));
    if (!isNull(gr)) {
      if (!isFunction(gr)) error(_("'gr' is not a function"));
      PROTECT(OS->R_gcall = lang2(gr, R_NilValue));
    } else {
      PROTECT(OS->R_gcall = R_NilValue); /* for balance */
ndeps = getListElement(options, "ndeps");
if (LENGTH(ndeps) != npar)
  error(_("'ndeps' is of the wrong length"));
OS->ndeps = vect(npar);
PROTECT(ndeps = coerceVector(ndeps, REALSXP));
for (i = 0; i < npar; i++) OS->ndeps[i] = REAL(ndeps)[i];
UNPROTECT(1);
    }
    cgmin(npar, dpar, opar, &val, fminfn, fmingr, &ifail, abstol,
          reltol, (void *)OS, type, trace, &fncount, &grcount, maxit);
    for (i = 0; i < npar; i++)
      REAL(par)[i] = opar[i] * (OS->parscale[i]);
    UNPROTECT(1); /* OS->R_gcall */

  } else if (strcmp(tn, "L-BFGS-B") == 0) {
    SEXP ndeps, smsg;
    double *lower = vect(npar), *upper = vect(npar);
    int lmm, *nbd = (int *) R_alloc(npar, sizeof(int));
    double factr, pgtol;
    char msg[60];

    nREPORT = asInteger(getListElement(options, "REPORT"));
    factr = asReal(getListElement(options, "factr"));
    pgtol = asReal(getListElement(options, "pgtol"));
    lmm = asInteger(getListElement(options, "lmm"));
    if (!isNull(gr)) {
      if (!isFunction(gr)) error(_("'gr' is not a function"));
      PROTECT(OS->R_gcall = lang2(gr, R_NilValue));
    } else {
      PROTECT(OS->R_gcall = R_NilValue); /* for balance */
ndeps = getListElement(options, "ndeps");
if (LENGTH(ndeps) != npar)
  error(_("'ndeps' is of the wrong length"));
OS->ndeps = vect(npar);
PROTECT(ndeps = coerceVector(ndeps, REALSXP));
for (i = 0; i < npar; i++) OS->ndeps[i] = REAL(ndeps)[i];
UNPROTECT(1);
    }
    args = CDR(args); slower = CAR(args); /* coerce in calling code */
args = CDR(args); supper = CAR(args);
for (i = 0; i < npar; i++) {
  lower[i] = REAL(slower)[i] / (OS->parscale[i]);
  upper[i] = REAL(supper)[i] / (OS->parscale[i]);
  if (!R_FINITE(lower[i])) {
    if (!R_FINITE(upper[i])) nbd[i] = 0; else nbd[i] = 3;
  } else {
    if (!R_FINITE(upper[i])) nbd[i] = 1; else nbd[i] = 2;
  }
}
OS->usebounds = 1;
OS->lower = lower;
OS->upper = upper;
lbfgsb(npar, lmm, dpar, lower, upper, nbd, &val, fminfn, fmingr,
       &ifail, (void *)OS, factr, pgtol, &fncount, &grcount,
       maxit, msg, trace, nREPORT);
for (i = 0; i < npar; i++)
  REAL(par)[i] = dpar[i] * (OS->parscale[i]);
UNPROTECT(1); /* OS->R_gcall */
PROTECT(smsg = mkString(msg));
SET_VECTOR_ELT(res, 4, smsg);
UNPROTECT(1);
  } else
    error(_("unknown 'method'"));

  if(!isNull(OS->names)) setAttrib(par, R_NamesSymbol, OS->names);
  REAL(value)[0] = val * (OS->fnscale);
  SET_VECTOR_ELT(res, 0, par); SET_VECTOR_ELT(res, 1, value);
  INTEGER(counts)[0] = fncount; INTEGER(counts)[1] = grcount;
  SET_VECTOR_ELT(res, 2, counts);
  INTEGER(conv)[0] = ifail;
  SET_VECTOR_ELT(res, 3, conv);
  UNPROTECT(6);
  return res;
}


/* par fn gr options */
SEXP optimhess(SEXP call, SEXP op, SEXP args, SEXP rho)
{
  SEXP par, fn, gr, options, tmp, ndeps, ans;
  OptStruct OS;
  int npar, i , j;
  double *dpar, *df1, *df2, eps;

  args = CDR(args);
  OS = (OptStruct) R_alloc(1, sizeof(opt_struct));
  OS->usebounds = 0;
  OS->R_env = rho;
  par = CAR(args);
  npar = LENGTH(par);
  OS->names = getAttrib(par, R_NamesSymbol);
  args = CDR(args); fn = CAR(args);
  if (!isFunction(fn)) error(_("'fn' is not a function"));
  args = CDR(args); gr = CAR(args);
  args = CDR(args); options = CAR(args);
  OS->fnscale = asReal(getListElement(options, "fnscale"));
  tmp = getListElement(options, "parscale");
  if (LENGTH(tmp) != npar)
    error(_("'parscale' is of the wrong length"));
  PROTECT(tmp = coerceVector(tmp, REALSXP));
  OS->parscale = vect(npar);
  for (i = 0; i < npar; i++) OS->parscale[i] = REAL(tmp)[i];
  UNPROTECT(1);
  PROTECT(OS->R_fcall = lang2(fn, R_NilValue));
  PROTECT(par = coerceVector(par, REALSXP));
  if (!isNull(gr)) {
    if (!isFunction(gr)) error(_("'gr' is not a function"));
    PROTECT(OS->R_gcall = lang2(gr, R_NilValue));
  } else {
    PROTECT(OS->R_gcall = R_NilValue); /* for balance */
  }
  ndeps = getListElement(options, "ndeps");
  if (LENGTH(ndeps) != npar) error(_("'ndeps' is of the wrong length"));
  OS->ndeps = vect(npar);
  PROTECT(ndeps = coerceVector(ndeps, REALSXP));
  for (i = 0; i < npar; i++) OS->ndeps[i] = REAL(ndeps)[i];
  UNPROTECT(1);
  PROTECT(ans = allocMatrix(REALSXP, npar, npar));
  dpar = vect(npar);
  for (i = 0; i < npar; i++)
    dpar[i] = REAL(par)[i] / (OS->parscale[i]);
  df1 = vect(npar);
  df2 = vect(npar);
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
    SEXP nm = getAttrib(par, R_NamesSymbol);
  if(!isNull(nm)) {
    SEXP dm;
    PROTECT(dm = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(dm, 0, duplicate(nm));
    SET_VECTOR_ELT(dm, 1, duplicate(nm));
    setAttrib(ans, R_DimNamesSymbol, dm);
    UNPROTECT(1);
  }
  UNPROTECT(4);
  return ans;
}

