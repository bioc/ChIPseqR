/*
 * startScore.h
 *
 *  Created on: 01/12/2008
 *      Author: Peter Humburg
 */

#ifndef STARTSCORE_H_
#define STARTSCORE_H_

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Utils.h>

R_INLINE double _ipoisLRS(int, int, int, int);
R_INLINE double _dpoisLRS(double, double, double, double);
R_INLINE double _calcStat(int, int, int, int, int, int);
R_INLINE double _calcStat2(int, double, double, int, double, int, int, int);
R_INLINE double _ratioStat_pois(int, int, int, int, int, int, int, int, double, int);

SEXP startScore_pois(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

#endif /* STARTSCORE_H_ */
