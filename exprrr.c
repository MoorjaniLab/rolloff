#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <nicklib.h>
#include <getpars.h>

#include "badpairs.h"
#include "admutils.h"
#include "mcio.h"  
#include "mcmcpars.h"  
#include "egsubs.h"  
#include "../qpsubs.h"  

#define WVERSION   "120" 

/** 
 fit exponential

 New I/O (mcio.c) added
 Fits regression (AR model)  
*/


#define MAXFL  50   
#define MAXSTR  512
#define MAXPOPS  100
#define NEXP  1  


extern int packmode ;

int popcheck = YES ;
int fancynorm = YES ; 
int plotmode = NO ;
int outnum = -1 ;
int affmode = NO ;
int helpmode = NO ;
int pop2mode = NO  ;
int popallmode = NO  ;
int numexp = NEXP  ;
int nlag = 10 ;
char *iname = NULL ;
char *oname = NULL ;

int datacol = 1 ;
double loval = -1.0e20 ;
double hival = 1.0e20 ;
double addx  = 0 ;
double step = -1 ;

int popsizelimit = -1 ;
int xscratch[50000] ; 
double yscratch[100000] ;

char *trashdir = "/var/tmp" ;
int verbose = NO ;
int qtmode = NO ;
Indiv **indivmarkers;
SNP **snpmarkers ;
int numsnps, numindivs ; 

char  *genotypename = NULL ;
char  *snpname = NULL ;
char  *genooutfilename = NULL ;
char  *indoutfilename = NULL ;
char  *indivname = NULL ;
char *badsnpname = NULL ;
char *goodsnpname = NULL ;
char *badpairsname = NULL ;
char *markername = NULL ;
char *poplistname = NULL ;

char *outputname = NULL ;
FILE *ofile ;

int xchrom = -1 ;
int *fchrom, *lchrom ;

int details = NO ;
double fakespacing = 0.0 ;
double dthresh =  0.02 ;
int    heminum  = 50 ;
int    homonum  = 50 ;
double blgsize  = .001 ;  // morgans
double blglim   = .10  ;  // morgans
long seed = 0 ;

char  unknowngender = 'U' ;

void readcommands(int argc, char **argv) ;
int dohemiz(Indiv *ind1, Indiv *ind2, SNP **snpmarkers, int numsnps, SNP ***runlist, int *lrunlist, double dthresh)  ;
void pubrun(SNP ***runlist, SNP *slo, SNP *shi, int *pnruns, double dthresh) ;
void printruns(SNP ***runlist, int nruns, Indiv *ind1, Indiv *ind2) ;
double fitexp(double *xval, double *xfit, double *xexp, double *xco, int m, int n, int mode, int inititer) ;
double fitrrr(double *xval, double *xfit, double *xexp, double *xco, int m, int n, int mode, int inititer) ;

void dofit(double **xx, int n)  ;
void printhelp() ;

int main(int argc, char **argv)
{

  int i, j, k, g, t, ipop, jpop, kpop ; 
  int k1, k2 ;
  int len ;

  int ch1, ch2 ;

  int numvind, nignore, numrisks = 1 ;
  int numcols ;
  int maxbin ;
  double y, gdis, ysd ;
  double **xx ;  
  char *pop1, *pop2 ;
  char popstring[MAXSTR] ;   

  double *xfit, *xexp, *xco, *xval, *xbase ;
  int nexp, lenxfit ;

  ofile = stdout; 
  readcommands(argc, argv) ;
  if (helpmode) {
    printhelp() ;
    return 1 ; 
  } 

  if (iname==NULL) {
     printhelp() ;
     return 0 ;
  }
  if (oname != NULL) openit(oname, &ofile, "w") ;

  nexp = numexp ;
  maxbin = numlines(iname) ; 

  numcols = datacol+1 ; 
  xx = initarray_2Ddouble(numcols, maxbin,  0.0) ;

  maxbin = getxx(xx, maxbin, numcols, iname) ;

  xbase = xx[0] ;
  xval = xx[datacol] ;
  lenxfit = maxbin ; 
  for (;;) { 
   y = xbase[0] ; 
   if (y>=loval) break ;
   ++xbase ; 
   ++xval ; 
   --lenxfit ;
   if (lenxfit <= 0) fatalx("no data\n") ;
  }

  len = 0 ;
  for (k=0; k<lenxfit; ++k) {
   y = xbase[k] ; 
   if (y>hival) break ;
   ++len ;
  }
  if (len<50) fatalx("data too short %d\n", len) ;
  lenxfit = len-10 ;

  vsp(xval, xval, addx, lenxfit) ;
  ZALLOC(xfit, lenxfit, double) ;
  ZALLOC(xexp, nexp+1, double) ;
  ZALLOC(xco, nexp+1, double) ;
  k = numexp ;
   if (affmode) printf("fitting %d exponentials + affine\n", k) ;
   if (!affmode) printf("fitting %d exponentials\n", k) ;
    ysd = fitrrr(xval, xfit, xexp, xco, lenxfit, nlag, affmode, 5000) ;
    printf("exp, co nexp: %d\n", nexp ) ;
    printmatl(xexp, 1, nexp) ;
    printmatl(xco, 1, nexp+1) ;
    for (i=0; i<nexp; i++) { 
     xexp[i] = hlife(xexp[i]) ; 
    }
    printf("error sd: %12.6f\n", ysd) ;
    printf("halflife: ") ;
    printmat(xexp, 1, nexp) ;
    if (step>0) {  
     for (i=0; i<nexp; i++) { 
      y = xexp[i]*step;   
      y += 1.0e-20 ;
      xexp[i] = log(2.0) / y ;
     }
     printf("mean (generations): ") ;
     printmat(xexp, 1, nexp) ;
    }
    fprintf(ofile, "##fit: %s\n", iname) ;
    for (i=0; i<lenxfit; ++i)  { 
     fprintf(ofile, "%12.6f ", xbase[i]) ;
     fprintf(ofile, "%12.6f ", xval[i]) ;
     fprintf(ofile, "%12.6f ", xfit[i]) ;
     fprintf(ofile, "%12.6f ", xval[i] - xfit[i]) ;
     fprintf(ofile, "\n") ;
    }
  printf("##end of run\n") ;
  return 0 ;
}
void printhelp() 
{ 
 printf("exprrr: \n") ;
 printf (" -i iname   ## input\n") ;
 printf (" -o oname   ## output\n") ;
 printf (" -c col     ## data column  (0 is xval)\n") ;
 printf (" -l loval   ## lowest x value\n") ;
 printf (" -h hival   ## lowest x value\n") ;
 printf (" -x val     ## value to add to x\n") ;
 printf (" -r ran     ## seed for random generator\n")  ;
 printf (" -b nlag    ## for AR model\n")  ;
 printf (" -a         ## affine mode (add constant)\n") ;
 printf (" -V         ## verbose mode (add constant)\n") ;
 printf (" -m         ## print help menu and quit\n") ;

}


void readcommands(int argc, char **argv) 

{
  int i,haploid=0;
  char *parname = NULL ;
  int n ;

  while ((i = getopt (argc, argv, "c:i:n:r:s:o:l:h:x:b:Vam")) != -1) {

    switch (i)
      {

      case 'i':
	iname = strdup(optarg) ;
	break;

      case 'r':
	seed = atoi(optarg) ;
	break;

      case 'c':
	datacol = atoi(optarg) ;
	break;

      case 's':
	step = atof(optarg) ; // morgans
	break;

      case 'a':
	affmode = YES ;
	break; 

      case 'm':
	helpmode = YES ;
	break; 

      case 'o':
	oname = strdup(optarg) ;
	break; 

      case 'l':
	loval = atof(optarg) ;
	break; 

      case 'h':
	hival = atof(optarg) ;
	break; 

      case 'x':
	addx = atof(optarg) ;
	break; 

      case 'b':
	nlag = atoi(optarg) ;
	break; 


      case 'V':
	verbose = YES ;
	break; 

      default: 
	printf ("Usage: bad params.... \n") ;
	fatalx("bad params\n") ;
      }
  }
         
}

