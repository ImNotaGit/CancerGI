#ifdef _OPENMP
#include <omp.h>
#endif
#include <Rcpp.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#define MIN(a,b) a<b?a:b 
#define MAX(a,b) a>b?a:b 
#define ABS(a) (a>0.0) ? a:-a

using namespace Rcpp;		// shorthand
using namespace std;		// shorthand
Rcpp::NumericMatrix binomialTest(Rcpp::IntegerVector scna1, Rcpp::IntegerMatrix scnaq){
	int numGenes = scnaq.nrow();
	int numSamples = scnaq.ncol();
	double numSamples2 =numSamples * numSamples;
	NumericMatrix binopscna(numGenes,9);
	IntegerVector scna2;
	IntegerVector counts(9);
	IntegerVector inx(numSamples);
	IntegerVector gene1cnt(3), gene2cnt(3);
	IntegerVector genecnt(9);
	NumericVector p(9);
	int xx;
	for (int gene2 = 0; gene2< numGenes; gene2++)
	{
		scna2 = scnaq(gene2,_);
		inx = 3 * scna2 + scna1;
		std::fill(counts.begin(), counts.end(), 0);
		std::fill(gene1cnt.begin(), gene1cnt.end(), 0);
		std::fill(gene2cnt.begin(), gene2cnt.end(), 0);

		for (int samp = 0; samp < numSamples; ++samp)
		{
			counts[inx[samp]]++;
		}
		
		for (int row = 0; row < 3; row++)
		{
			for (int col = 0; col < 3; col++)
			{
				xx = row *3 + col;
				gene1cnt(col) = gene1cnt(col) + counts(xx);
				gene2cnt(row) = gene2cnt(row) + counts(xx);
			}
		}
		for (int row = 0; row < 3; row++)
		{
			for (int col = 0; col < 3; col++)
			{
				genecnt(row*3 + col) = gene1cnt(col) * gene2cnt(row);
			}
		}
		for (int row = 0; row < 9; row++){
	        double exp = ((double)  genecnt(row))/numSamples2;
	        double obs = ((double) counts(row))/numSamples;
	        if( obs < exp ) {
    	      binopscna(gene2,row) = (-1)*(R::pbinom(counts(row), numSamples, exp, 1, 0));
	        } else {
    	      binopscna(gene2,row) = (R::pbinom(counts(row), numSamples, exp, 0, 0));
    	    }
       }
	}
	return binopscna;
}
