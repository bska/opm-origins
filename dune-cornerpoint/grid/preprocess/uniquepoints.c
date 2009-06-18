#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <stdio.h>


#include "sparsetable.h"
#include "grdecl.h"
#include "uniquepoints.h"
#include "matalloc.h"

#define min(i,j) ((i)<(j) ? (i) : (j))
#define max(i,j) ((i)>(j) ? (i) : (j))
#define overlap(a1,a2,b1,b2) max(a1,b1) < min(a2,b2)


/* Compar function passed to qsort                       */
/*-------------------------------------------------------*/
static int compare(const void *a, const void *b)
{
  if (*(double*)a < *(double*) b) return -1;
  else return 1;
}

/* Creat sorted list of z-values in zcorn with actnum==1 */
/*-------------------------------------------------------*/
static int createSortedList(double *list, int n, int m, 
			    const double *z[], const int *a[])
{
  int i,j;
  double *ptr = list;
  for (i=0; i<n; ++i){
    for (j=0; j<m; ++j){
      if (a[j][i/2])  *ptr++ = z[j][i];
      /* else        fprintf(stderr, "skipping point in inactive cell\n"); */
    }
  }

  qsort(list, ptr-list, sizeof(double), compare);
  return ptr-list;
}


/* Remove points that are closer than tolerance in list  */
/* of increasing doubles                                 */
/*-------------------------------------------------------*/
static int uniquify(int n, double *list, double tolerance)
{
  if (n<1) return 0;
  int    i;
  int    pos = 0; 
  double val = list[pos++];/* Keep first value */
  
  for (i=1; i<n; ++i){
    if (list[i] - val > tolerance){
      val         = list[i];
      list[pos++] = val;
    }
  }

  /* Keep last value (one way or the other...) */
  if (list[n-1] - val > tolerance){
    list[pos-1] = list[n-1];
  }
  
  return pos;
}


/* Along single pillar: */
static int assignPointNumbers(int    begin, 
			       int    end, 
			       const double *zlist,
			       int    n, 
			       const double *zcorn,
			       const int    *actnum, 
			       int    *plist,
			       double tolerance)
{
  /* n     - number of cells */
  /* zlist - list of len unique z-values */
  /* start - number of unique z-values processed before. */

  int i, k;
  /* All points should now be within tolerance of a listed point. */


  const double *z = zcorn;
  const int    *a = actnum;
  int    *p = plist;

  k = begin;
  *p++ = INT_MIN; /* Padding to ease processing of faults */
  for (i=0; i<n; ++i){

    /* Skip inactive cells */
    if (!a[i/2]) {
      p[0] = p[-1];  /* Inactive cells are collapsed leaving void space.*/
      ++p;
      continue;
    }

    /* Find next k such that zlist[k] < z[i] < zlist[k+1] */
    while (k < end && zlist[k] + tolerance < z[i]){
      k++;
    }

    /* assert (k < len && z[i] - zlist[k] <= tolerance) */
    if (k == end || z[i] - zlist[k] > tolerance){
      fprintf(stderr, "Cannot associate  zcorn values with given list\n");
      fprintf(stderr, "of z-coordinates to given tolerance\n");      
      return 0;
    }

    *p++ = k;
  }
  *p++ = INT_MAX;/* Padding to ease processing of faults */


  return 1;
}


/*-------------------------------------------------------*/
static void igetvectors(const int dims[3], int i, int j, 
			const int *field, const int *v[])
{
  
  int im = max(1,       i  ) - 1;
  int ip = min(dims[0], i+1) - 1;
  int jm = max(1,       j  ) - 1;
  int jp = min(dims[1], j+1) - 1;
  
  v[0] = field + dims[2]*(im + dims[0]* jm);
  v[1] = field + dims[2]*(im + dims[0]* jp);
  v[2] = field + dims[2]*(ip + dims[0]* jm);
  v[3] = field + dims[2]*(ip + dims[0]* jp);
}

/*-------------------------------------------------------*/
static void dgetvectors(const int dims[3], int i, int j, const double *field, const double *v[])
{
  
  int im = max(1,       i  ) - 1;
  int ip = min(dims[0], i+1) - 1;
  int jm = max(1,       j  ) - 1;
  int jp = min(dims[1], j+1) - 1;
  
  v[0] = field + dims[2]*(im + dims[0]* jm);
  v[1] = field + dims[2]*(im + dims[0]* jp);
  v[2] = field + dims[2]*(ip + dims[0]* jm);
  v[3] = field + dims[2]*(ip + dims[0]* jp);
}

/* Assign point numbers p such that "zlist(p)==zcorn".   */
/* Assume that coordinate number is arranged in a        */
/* sequence such that the natural index is (k,i,j)       */
/*-------------------------------------------------------*/
int finduniquepoints(const struct grdecl *g,
		                            /* return values: */
		     int           *plist, /* list of point numbers on each pillar*/
		     sparse_table_t *ztab)
		      
{


  double *zlist = ztab->data; /* casting void* to double* */
  int     *zptr = ztab->ptr;

  int     i,j;

  int     d1[3]  = {2*g->dims[0], 2*g->dims[1], 2*g->dims[2]};
  int     len    = 0;
  double  *zout  = zlist;  
  int     pos    = 0;

  zptr[pos++] = zout - zlist;

  /* Loop over pillars, find unique points on each pillar */
  for (j=0; j < g->dims[1]+1; ++j){
    for (i=0; i < g->dims[0]+1; ++i){

      const int    *a[4];
      const double *z[4];
      
      /* Get positioned pointers for actnum and zcorn data */
      igetvectors(g->dims,   i,   j, g->actnum, a);
      dgetvectors(d1,      2*i, 2*j, g->zcorn,  z);

      len = createSortedList(     zout, d1[2], 4, z, a);
      len = uniquify        (len, zout, DBL_EPSILON);      

      /* Increment pointer to sparse table of unique zcorn values */
      zout        = zout + len;
      zptr[pos++] = zout - zlist;
    }
  }



  /* Loop over all vertical sets of zcorn values, assign point numbers */
  int *p = plist;
  for (j=0; j < 2*g->dims[1]; ++j){
    for (i=0; i < 2*g->dims[0]; ++i){
      
      /* pillar index */
      int pix = (i+1)/2 + (g->dims[0]+1)*((j+1)/2);
      
      /* cell column position */
      int cix = g->dims[2]*((i/2) + (j/2)*g->dims[0]);

      /* zcorn column position */
      int zix = 2*g->dims[2]*(i+2*g->dims[0]*j);
      
      const int    *a   = g->actnum + cix; 
      const double *z   = g->zcorn  + zix;

      if (!assignPointNumbers(zptr[pix], zptr[pix+1], zlist,
			      2*g->dims[2], z, a, p, DBL_EPSILON)){
	fprintf(stderr, "Something went wrong in assignPointNumbers");
	return 0;
      }

      p += 2 + 2*g->dims[2];
    }
  }
  return 1;
}


