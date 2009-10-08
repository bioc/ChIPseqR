/*
 * startScore.c
 *
 *  Created on: 01/12/2008
 *      Author: Peter Humburg
 */

#include "include/startScore.h"

/* Likelihood ratio statistic for two Poisson distributions
 * x0, x1 are the observed number of events from each distribution
 * t0, t1 are the corresponding observation times
 */
double _dpoisLRS(double x0, double x1, double t0, double t1){
	double stat;
	if(x0 == 0)
		stat = 2 * x1 * (log(t0 + t1) - log(t1));
	else if(x1 == 0)
		stat = 2 * x0 * (log(t0 + t1) - log(t0));
	else
		stat = 2.0*(x0*(log(x0) - log(t0)) + x1*(log(x1) - log(t1)) + (x0 + x1)*(log(t0 + t1) - log(x0 + x1)));
	if(stat < 0) stat = -1 * stat;
	return stat;
}

/* Calculate likelihood ratio statistic comparing the rates of two Poisson distributions.
 * This is an interface to _dpoisLRS allowing for integer arguments. */
double _ipoisLRS(int x0, int x1, int t0, int t1){
	return _dpoisLRS((double) x0, (double) x1, (double) t0, (double) t1);
}

/* Calculate (square root of) likelihood ratio statistic */
double _calcStat(int support, int background, int supLen, int bgLen, int sgn, int verbose){
	double stat;

	stat = _ipoisLRS(support, background, supLen, bgLen);
	if(verbose > 0) Rprintf("%f\n", stat);
	if((double)(support * bgLen)/(background * supLen) < 1) sgn = -1 * sgn;

	return sqrt(stat) * sgn;
}

/* Calculate (square root of) likelihood ratio statistic with capped rates*/
double _calcStat2(int support, double supRate, double supCut, int background, double bgRate, int supLen, int bgLen, int sgn){
	double stat, maxRate;

	if(supRate == 0 || support == 0) stat = NA_REAL;
	else{
		maxRate = qpois(supCut, supRate, TRUE, FALSE);
		if(support > maxRate) support = (int)maxRate;
		stat = _ipoisLRS(support, background, supLen, bgLen);

		if(supRate < bgRate) sgn = -1 * sgn;
		stat = sqrt(stat) * sgn;
	}

	return stat;
}

/* Compute likelihood ratio statistic based score for binding sites. Statistics for both support regions and the binding
 *  region are computed. The final statistic is the sum of the square roots. */
double _ratioStat_pois(int fwdSupCount, int bindSupCount, int revSupCount, int fwdBgCount, int revBgCount,
		int totalLen, int b, int supLen, double supCut, int verbose){

	int bgTotal = fwdBgCount + revBgCount;
	double stat, tmp_stat;

	int bgLen = totalLen - b - supLen;

	if(fwdBgCount == 0 || revBgCount == 0){
		stat = NA_REAL;
	}
	else{
		stat = 0;
		/* forward support region*/
		if(verbose > 0){
			Rprintf("forward\n");
			Rprintf("%d, %d, %d, %d\n", fwdSupCount, fwdBgCount, supLen, bgLen);
		}
		stat = _calcStat2(fwdSupCount, (double)revSupCount, supCut, fwdBgCount, (double) supLen*fwdBgCount/bgLen, supLen, bgLen, 1);
		if(verbose > 0) Rprintf("%f\n\n", stat);

		/* reverse support region*/
		if(verbose > 0){
			Rprintf("reverse\n");
			Rprintf("%d, %d, %d, %d\n", revSupCount, revBgCount, supLen, bgLen);
		}
		stat += _calcStat2(revSupCount, (double)fwdSupCount, supCut, revBgCount, (double) supLen*revBgCount/bgLen, supLen, bgLen, 1);
		if(verbose > 0) Rprintf("%f\n\n", _calcStat(revSupCount, revBgCount, supLen, bgLen, 1, 0));

		/* binding region*/
		if(verbose > 0){
			Rprintf("binding\n");
			Rprintf("%d, %d, %d, %d\n", bindSupCount, bgTotal, 2*b, 2*bgLen);
		}
		stat += _calcStat(bindSupCount, bgTotal, 2*b, 2*bgLen, -1, verbose);
		if(verbose > 0) {
			Rprintf("%f\n", _calcStat(bindSupCount, bgTotal, 2*b, 2*bgLen, -1, 0));
			Rprintf("%f\n", stat);
		}
	}

	return stat;
}

/* Calculate binding site score based on Poisson model */
SEXP startScore_pois(SEXP data, SEXP b, SEXP support, SEXP background, SEXP cutoff_bg, SEXP cutoff_sup){
	SEXP score;
	int len = nrows(data);
	int bgStart = 0, end, fwd, rev, fwdBgRate, revBgRate;

	int sup = INTEGER(support)[0];
	int bg = INTEGER(background)[0];
	int bind = INTEGER(b)[0];

	double supCut = REAL(cutoff_sup)[0], bgCut = REAL(cutoff_bg)[0], maxBg;

	PROTECT(score = allocVector(REALSXP, len - bind - 2*sup));

	/* initialise counts */
	int fwdSupCount = 0, revSupCount = 0, bindCount = 0, fwdBgCount = 0, revBgCount = 0, fwdBgCountOld, revBgCountOld;
	for(int i = 0; i < sup; i++){
		fwdSupCount += INTEGER(data)[i];
		revBgCount += INTEGER(data)[len + i];
		revSupCount += INTEGER(data)[len + sup + bind + i];
	}

	for(int i = sup; i < sup + bind; i++){
		bindCount += INTEGER(data)[i] + INTEGER(data)[i + len];
	}
	for(int i = sup+bind; i < bg; i++){
		fwdBgCount += INTEGER(data)[i];
		revBgCount += INTEGER(data)[i+ len];
	}
	fwdBgCountOld = fwdBgCount;
	revBgCountOld = revBgCount;
	fwdBgRate = fwdBgCount;
	revBgRate = revBgCount;

	/* calculate first score */
	REAL(score)[0] = _ratioStat_pois(fwdSupCount, bindCount, revSupCount, fwdBgCount, revBgCount, bg, bind, sup, supCut, 0);


	/* update counts and calculate remaining scores
	 * Scores are for binding sites starting at position i */
	for(int i = sup + 1; i < len - bind - sup; i++){
		R_CheckUserInterrupt();
		end = i + bind - 1;
		fwd = i - sup + 1;
		rev = end + sup - 1;

		/* update support regions*/
		fwdSupCount += INTEGER(data)[i-1] - INTEGER(data)[fwd-1];
		revSupCount += INTEGER(data)[rev + len] - INTEGER(data)[end + len];

		/* update binding site*/
		bindCount += INTEGER(data)[end] + INTEGER(data)[end + len] - INTEGER(data)[i-1] - INTEGER(data)[i-1 + len];

		/* update background */
		fwdBgCount += INTEGER(data)[fwd-1] - INTEGER(data)[end];
		revBgCount += INTEGER(data)[i-1 + len] - INTEGER(data)[rev + len];
		if(i + bind/2 + 1 > bgStart + bg/2 && bgStart + bg + 1 < len){
			bgStart++;
			fwdBgCount +=  INTEGER(data)[bgStart+bg] - INTEGER(data)[bgStart - 1];
	        revBgCount +=  INTEGER(data)[bgStart + bg + len] - INTEGER(data)[bgStart - 1 + len];
	        fwdBgRate = fwdBgCount;
	        revBgRate = revBgCount;
		}

		/* limit changes in background rate between adjacent windows */
		if(bgStart > bg){
			fwdBgCountOld += INTEGER(data)[bgStart - 1] - INTEGER(data)[bgStart - bg - 1];
			revBgCountOld +=  INTEGER(data)[bgStart -1 + len] - INTEGER(data)[bgStart - bg - 1 + len];
		}
		maxBg = qpois(bgCut, (double)fwdBgCountOld, TRUE, FALSE);
		if(fwdBgRate > maxBg){
			fwdBgRate = maxBg;
		}
		maxBg = qpois(bgCut, (double)revBgCountOld, TRUE, FALSE);
		if(revBgRate > maxBg){
			revBgRate = maxBg;
		}

		/* calculate score */
		REAL(score)[i - sup] = _ratioStat_pois(fwdSupCount, bindCount, revSupCount, fwdBgRate,
				revBgRate, bg, bind, sup, supCut,0);
	}

	UNPROTECT(1); /* score */
	return score;
}
