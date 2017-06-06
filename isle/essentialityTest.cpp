#ifdef _OPENMP
#include <omp.h>
#endif
#include <RcppArmadillo.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>      /* constants */
#define MIN(a,b) a<b?a:b 
#define MAX(a,b) a>b?a:b 
#define ABS(a) (a>0.0) ? a:-a

using namespace std;        // shorthand
using namespace arma;       // shorthand
using namespace Rcpp;       // shorthand
double  myRanksum(arma::vec x, arma::vec y) {
 // alternative hypothesis x > y 
 // x and y are assumed to be sorted. 
 int xlen = x.size();
 int ylen = y.size();
 if (!((xlen >0) & (ylen > 0)))
    return -1000;
 vec xy = join_vert(x,y);
 uvec sortxy = sort_index(xy); 
 double ranksumx = 0, ranksumy =0;
 double curr, next;
 int currCount =-1;
 double xinacc =0;
 vec temp;
 double xcnt=0; double ycnt =0;
 double acc =0;
 double pval;
 for (int i = 0; i < xlen + ylen; ++i)
 {
    curr = xy(sortxy(i));
    if(i < ((xlen + ylen) -1) )
        next = xy(sortxy(i +1));
    else 
        next = -1000;

    currCount++;
    if(sortxy(i) < xlen) 
        xcnt++;
    else 
        ycnt++;

    acc = acc + i + 1;
    
    if(curr != next){
        xinacc =   acc * xcnt/(xcnt + ycnt);   
        ranksumx = ranksumx + xinacc;
        acc = 0;
        xcnt =0; ycnt =0;
    }
     
 }
// if sample size is less than 10 use pwilcox otherwise use normal approximation
 if ((xlen <= 9) & (ylen <=9))
 {
 ranksumx= ranksumx - xlen*(xlen +1)/2.0;
    pval = R::pwilcox(ranksumx,  xlen,  ylen, 0, 0);
 } else{
    double mu =  xlen* (xlen + ylen +1) /2.0; 
    double sigma =  sqrt(ylen) * sqrt(mu/6.0);
    pval = R::pnorm(ranksumx,  mu,  sigma, 0, 0);

 }
 return pval;
}

using namespace Rcpp;       // shorthand
using namespace arma;       // shorthand
using namespace std;        // shorthand
arma::mat essentialityTest(arma::uvec scna1, arma::uvec mRNA1, arma::mat scnaq, arma::mat mRNAq, arma::mat ess, int threads=1){
   if ( threads > 0 )
    omp_set_num_threads( threads ); 
    int numGenes = scnaq.n_cols;
    mat out(numGenes,2);
    uvec grp1, grp2, grp3, grp4;
    grp1 = find(scna1 ==0);
    grp2 = find(mRNA1 ==0);
    grp3 = find(scna1 ==2);
    grp4 = find(mRNA1 ==2);
    int gene2;
  #pragma omp parallel for schedule(static)
    for ( gene2 = 0; gene2< numGenes; gene2++)
    {
        vec ess1 = ess.col(gene2);
        rowvec temp(2);
        temp( 0) = myRanksum(ess1(grp3), ess1(grp1));
        temp( 1) = myRanksum(ess1(grp4), ess1(grp2));
            out.row(gene2 ) = temp; 
    }
    return out;
}

using namespace Rcpp;       // shorthand
using namespace arma;       // shorthand
using namespace std;        // shorthand
arma::mat essentialityTestAll(arma::mat scnaq, arma::mat mRNAq, arma::mat ess, int threads=1){
   if ( threads > 0 )
    omp_set_num_threads( threads ); 
    int numGenes = scnaq.n_cols;
    mat out(numGenes*numGenes,2);
    int gene1, gene2;
    uvec grp1, grp2, grp3, grp4;
    for ( gene1 = 0; gene1 < numGenes; gene1++ )
    {
    vec mRNA1=mRNAq.col(gene1);
    vec scna1=scnaq.col(gene1);
    grp1 = find(scna1 ==0);
    grp2 = find(mRNA1 ==0);
    grp3 = find(scna1 ==2);
    grp4 = find(mRNA1 ==2);
  #pragma omp parallel for schedule(static)
    for ( gene2 = 0; gene2< numGenes; gene2++) if (gene1!=gene2)
    {
        vec ess1 = ess.col(gene2);
        rowvec temp(2);
        temp( 0) = myRanksum(ess1(grp3), ess1(grp1));
        temp( 1) = myRanksum(ess1(grp4), ess1(grp2));
            out.row(numGenes*(gene1)+gene2 ) = temp; 
    }}
    return out;
}
