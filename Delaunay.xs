#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

#include "ppport.h"

#include <string.h>

/* cross-platform macros to set floating point precision for doubles */
/* Mainly needed to restore defaults, since Triangle doesn't.        */
/* http://www.christian-seiler.de/projekte/fpmath/                   */

#include "xpfpa.h"

/* triangle.h needs these two defines                 */
/* REAL can be single or double, but I've used        */
/* double everywhere in the code here, so use double. */
/* VOID should be void. Triangle sets it to int to    */
/* "fool dumb compilers" but void looks like it's     */
/* more portable.                                     */
/* We also set REAL in the Build.PL compiler flags.   */

#define REAL double
#define VOID void

#include "triangle.h"

/* these let us refer to the struct in triangle.h    */

typedef struct triangulateio *Math__Geometry__Delaunay__TriangulateioPtr;
typedef struct triangulateio  Math__Geometry__Delaunay__Triangulateio;

/* for T_ARRAY typmap */

typedef double doubleArray;
typedef int       intArray;

/* for our custom T_ARRAY input code */

#define doubleArrayPtrContains double
#define intArrayPtrContains    int

/* memory allocation functions for T_ARRAY typemap INPUT */

doubleArray * doubleArrayPtr( int nelem ) {
    doubleArray * array;
    Newx(array, nelem, doubleArray);
    return array;
    }
intArray * intArrayPtr( int nelem ) {
    intArray * array;
    Newx(array, nelem, intArray);
    return array;
    }

/* TODO: */
/* Consider offering the -u option, the "user defined constraint" function feature.   */
/* Define a C function that calls a Perl function stub to be overidden. See perlcall. */


MODULE = Math::Geometry::Delaunay	PACKAGE = Math::Geometry::Delaunay	
PROTOTYPES: DISABLE

void
_triangulate(arg0, arg1, arg2, arg3)
    char * arg0
    Math::Geometry::Delaunay::TriangulateioPtr	arg1
    Math::Geometry::Delaunay::TriangulateioPtr	arg2
    Math::Geometry::Delaunay::TriangulateioPtr	arg3
    PREINIT:
      char * doVoronoi = 0;
      XPFPA_DECLARE() /* declares vars to stash the floating point config */
    CODE:
      doVoronoi = strchr(arg0, 'v');
      /* set normal double precision mode (triangle usually does this too) */
      XPFPA_SWITCH_DOUBLE()
      triangulate(arg0,arg1,arg2,arg3);
      /* set the floating point config back to whatever it was */
      XPFPA_RESTORE()

      /* Triangle just copies pointers from IN to OUT for these two lists of redundant info. */
      /* To enable safe garbage collection/DESTROY, we make a copy of the data for OUT.      */

      if (arg1->numberofholes   && arg1->holelist  ) {
          Newx(arg2->holelist, arg1->numberofholes   *  2, double);
          Copy(arg1->holelist, arg2->holelist , arg1->numberofholes   *  2, double);
          if (doVoronoi) {
              Newx(arg3->holelist, arg1->numberofholes   *  2, double);
              Copy(arg1->holelist, arg3->holelist , arg1->numberofholes   * 2, double);
              }
          }
      if (arg1->numberofregions && arg1->regionlist) {
          Newx(arg2->regionlist, arg1->numberofregions * (2 + arg1->numberoftriangleattributes), double);
          Copy(arg1->regionlist, arg2->regionlist , arg1->numberofregions * (2 + arg1->numberoftriangleattributes), double);
          if (doVoronoi) {
              Newx(arg3->regionlist, arg1->numberofregions * (2 + arg1->numberoftriangleattributes), double);
              Copy(arg1->regionlist, arg3->regionlist , arg1->numberofregions * (2 + arg1->numberoftriangleattributes), double);
              }
          }

MODULE=Math::Geometry::Delaunay      PACKAGE=Math::Geometry::Delaunay::Triangulateio      PREFIX=triio_
PROTOTYPES: DISABLE

Math::Geometry::Delaunay::Triangulateio
triio_new(char * CLASS)
    CODE:
        Zero((void*)&RETVAL, sizeof(RETVAL), char);
        /* initialize everything in the struct */
        RETVAL.pointlist                  = (double *) NULL;     /* In / out */
        RETVAL.pointattributelist         = (double *) NULL;     /* In / out */
        RETVAL.pointmarkerlist            = (int *)    NULL;     /* In / out */
        RETVAL.numberofpoints             = (int)         0;     /* In / out */
        RETVAL.numberofpointattributes    = (int)         0;     /* In / out */

        RETVAL.trianglelist               = (int *)    NULL;     /* In / out */
        RETVAL.triangleattributelist      = (double *) NULL;     /* In / out */
        RETVAL.trianglearealist           = (double *) NULL;     /* In only */
        RETVAL.neighborlist               = (int *)    NULL;     /* Out only */
        RETVAL.numberoftriangles          = (int)         0;     /* In / out */
        RETVAL.numberofcorners            = (int)         3;     /* In / out */
        RETVAL.numberoftriangleattributes = (int)         0;     /* In / out */

        RETVAL.segmentlist                = (int *)    NULL;     /* In / out */
        RETVAL.segmentmarkerlist          = (int *)    NULL;     /* In / out */
        RETVAL.numberofsegments           = (int)         0;     /* In / out */

        RETVAL.holelist                   = (double *) NULL;     /* In / pointer to array copied out */
        RETVAL.numberofholes              = (int)         0;     /* In / copied out */

        RETVAL.regionlist                 = (double *) NULL;     /* In / pointer to array copied out */
        RETVAL.numberofregions            = (int)         0;     /* In / copied out */

        RETVAL.edgelist                   = (int *)    NULL;     /* Out only */
        RETVAL.edgemarkerlist             = (int *)    NULL;     /* Not used with Voronoi diagram; out only */
        RETVAL.normlist                   = (double *) NULL;     /* Used only with Voronoi diagram; out only */
        RETVAL.numberofedges              = (int)         0;     /* Out only */
    OUTPUT: 
        RETVAL

# getter/setter for all the "number of" members of the triangulateio struct

int
numberof(THIS, newval = 0)
    Math::Geometry::Delaunay::Triangulateio THIS
    int newval
    PROTOTYPE: DISABLE
    ALIAS:
        numberofpoints             = 1
        numberofpointattributes    = 2
        numberoftriangles          = 3
        numberofcorners            = 4
        numberoftriangleattributes = 5
        numberofsegments           = 6
        numberofholes              = 7
        numberofregions            = 8
        numberofedges              = 9
    PREINIT:
        struct triangulateio * p;
        STRLEN len;
        char *s;
    CODE:
        if (!sv_derived_from(ST(0), "Math::Geometry::Delaunay::Triangulateio")) {croak("Wrong type to numberof()");} 
        s = SvPV((SV*)SvRV(ST(0)), len);
        if (len != sizeof(THIS)) {croak("Size %d of packed data != expected %d", len, sizeof(THIS));}
        p = (struct triangulateio *) s;
        switch (ix) {
            case 1 : if ( items > 1 ) {p->numberofpoints             = newval;} RETVAL = p->numberofpoints; break;
            case 2 : if ( items > 1 ) {p->numberofpointattributes    = newval;} RETVAL = p->numberofpointattributes; break;
            case 3 : if ( items > 1 ) {p->numberoftriangles          = newval;} RETVAL = p->numberoftriangles; break;
            case 4 : if ( items > 1 ) {p->numberofcorners            = newval;} RETVAL = p->numberofcorners; break;
            case 5 : if ( items > 1 ) {p->numberoftriangleattributes = newval;} RETVAL = p->numberoftriangleattributes; break;
            case 6 : if ( items > 1 ) {p->numberofsegments           = newval;} RETVAL = p->numberofsegments; break;
            case 7 : if ( items > 1 ) {p->numberofholes              = newval;} RETVAL = p->numberofholes; break;
            case 8 : if ( items > 1 ) {p->numberofregions            = newval;} RETVAL = p->numberofregions; break;
            case 9 : if ( items > 1 ) {p->numberofedges              = newval;} RETVAL = p->numberofedges; break;
            default : RETVAL = 0;
            }
    OUTPUT: 
        RETVAL

# getters/setters for all the arrays in the triangulateio struct

# This uses T_ARRAY typmap, with it's INPUT code overridden to accomodate
# always getting an object reference as the first, and sometimes the only, 
# argument.
#
# Note that the use of T_ARRAY requires certain typedefs and that the modified 
# INPUT code relies on a couple defines, that follow the typedefs above.
#
# Also, the modified INPUT code imposes a constraint on the maximum value 
# integers that can be stored in Triangulateio's lists. Integer input is 
# interpreted as a double, then cast back to an integer. If your integers have 
# more than 53 bits they will lose the least significant bits that don't fit in 
# a double.
#
# 53 bit integers let you have up to 9,007,199,254,740,992 indexable items of 
# any kind in your triangulations.

# getter/setter for all the arrays of doubles in the triangulateio struct

doubleArray *
doubleList(doubleArray * array, ... )
    ALIAS:
        pointlist             = 1
        pointattributelist    = 2
        triangleattributelist = 3
        trianglearealist      = 4
        holelist              = 5
        regionlist            = 6
        normlist              = 7
    PREINIT:
        U32 size_RETVAL;
        struct triangulateio * p;
        STRLEN len;
        char *s;
        int orig_items_cnt = (int) items;
    CODE:
        s = SvPV((SV*)SvRV(ST(0)), len);
        p = (struct triangulateio *) s;
        /* setter */
        if (orig_items_cnt > 1) {
          switch (ix) {
              case 1 :  if (p->pointlist)             {trifree(p->pointlist);}             p->pointlist             = array ; p->numberofpoints  = (int) (ix_array)/2; break;
              case 2 :  if (p->pointattributelist)    {trifree(p->pointattributelist);}    p->pointattributelist    = array ; break;
              case 3 :  if (p->triangleattributelist) {trifree(p->triangleattributelist);} p->triangleattributelist = array ; break;
              case 4 :  if (p->trianglearealist)      {trifree(p->trianglearealist);}      p->trianglearealist      = array ; break;
              case 5 :  if (p->holelist)              {trifree(p->holelist);}              p->holelist              = array ; p->numberofholes   = (int) (ix_array)/2; break;
              case 6 :  if (p->regionlist)            {trifree(p->regionlist);}            p->regionlist            = array ; p->numberofregions = (int) (ix_array)/2; break;
              case 7 :  if (p->normlist)              {trifree(p->normlist);}              p->normlist              = array ; break;
              }
          /* return count of how many items added */
          ST(0) = sv_newmortal();
          sv_setnv(ST(0), (double)(ix_array));
          XSRETURN(1);
          }
        /* getter */
        else {
          switch (ix) {
              case 1 : RETVAL = p->pointlist;             size_RETVAL = (U32) p->numberofpoints*2; break;
              case 2 : RETVAL = p->pointattributelist;    size_RETVAL = (U32) p->numberofpoints*p->numberofpointattributes; break;
              case 3 : RETVAL = p->triangleattributelist; size_RETVAL = (U32) p->numberoftriangles*p->numberoftriangleattributes; break;
              case 4 : RETVAL = p->trianglearealist;      size_RETVAL = (U32) p->numberoftriangles; break;
              case 5 : RETVAL = p->holelist;              size_RETVAL = (U32) p->numberofholes*2; break;
              case 6 : RETVAL = p->regionlist;            size_RETVAL = (U32) p->numberofregions*p->numberoftriangleattributes; break;
              case 7 : RETVAL = p->normlist;              size_RETVAL = (U32) p->numberofedges*2; break;
              default: RETVAL = (double *) NULL;          size_RETVAL = (U32) 0;
              }
          }
    OUTPUT: 
        RETVAL
    CLEANUP:
        XSRETURN(size_RETVAL);

# getter/setter for all the arrays of ints in the triangulateio struct

intArray *
intList( intArray * array, ... )
    ALIAS:
        pointmarkerlist     = 1
        trianglelist        = 2
        neighborlist        = 3
        segmentlist         = 4
        segmentmarkerlist   = 5
        edgelist            = 6
        edgemarkerlist      = 7
    PREINIT:
        U32 size_RETVAL;
        struct triangulateio * p;
        STRLEN len;
        char *s;
        int orig_items_cnt = (int) items;
    CODE:
        s = SvPV((SV*)SvRV(ST(0)), len);
        p = (struct triangulateio *) s;
        /* setter */
        if (orig_items_cnt > 1) {
            switch (ix) {
                case 1 :  if (p->pointmarkerlist)   {trifree(p->pointmarkerlist);}   p->pointmarkerlist = array;   break;
                case 2 :  if (p->trianglelist)      {trifree(p->trianglelist);}      p->trianglelist = array;      p->numberoftriangles = (ix_array)/3; break;
                case 3 :  if (p->neighborlist)      {trifree(p->neighborlist);}      p->neighborlist = array;      break;
                case 4 :  if (p->segmentlist)       {trifree(p->segmentlist);}       p->segmentlist = array;       p->numberofsegments  = (ix_array)/2; break;
                case 5 :  if (p->segmentmarkerlist) {trifree(p->segmentmarkerlist);} p->segmentmarkerlist = array; break;
                case 6 :  if (p->edgelist)          {trifree(p->edgelist);}          p->edgelist = array;          p->numberofedges     = (ix_array)/2; break;
                case 7 :  if (p->edgemarkerlist)    {trifree(p->edgemarkerlist);}    p->edgemarkerlist = array;    break;
                }
            /* return count of how many items added */
            ST(0) = sv_newmortal();
            sv_setiv(ST(0), (int)ix_array);
            XSRETURN(1);
            }
        /* getter */
        else {
            switch (ix) {
                case 1 : RETVAL = p->pointmarkerlist;   size_RETVAL = p->pointmarkerlist   ? p->numberofpoints      : 0; break;
                case 2 : RETVAL = p->trianglelist;      size_RETVAL = p->trianglelist      ? p->numberoftriangles * p->numberofcorners : 0; break;
                case 3 : RETVAL = p->neighborlist;      size_RETVAL = p->neighborlist      ? p->numberoftriangles*3 : 0; break;
                case 4 : RETVAL = p->segmentlist;       size_RETVAL = p->segmentlist       ? p->numberofsegments*2  : 0; break;
                case 5 : RETVAL = p->segmentmarkerlist; size_RETVAL = p->segmentmarkerlist ? p->numberofsegments    : 0; break;
                case 6 : RETVAL = p->edgelist;          size_RETVAL = p->edgelist          ? p->numberofedges*2     : 0; break;
                case 7 : RETVAL = p->edgemarkerlist;    size_RETVAL = p->edgemarkerlist    ? p->numberofedges       : 0; break;
                default: RETVAL = (int *) NULL;         size_RETVAL = 0;
                }
            }
    OUTPUT:
        RETVAL
    CLEANUP:
        XSRETURN(size_RETVAL);

Math::Geometry::Delaunay::TriangulateioPtr
to_ptr(Math::Geometry::Delaunay::Triangulateio THIS = NO_INIT)
    PREINIT:
        STRLEN len;
        char *s;
    CODE:
        if (sv_derived_from(ST(0), "Math::Geometry::Delaunay::Triangulateio")) {
            s = SvPV((SV*)SvRV(ST(0)), len);
            if (len != sizeof(THIS)) { croak("Size %d of packed data != expected %d", len, sizeof(THIS)); }
            RETVAL = (struct triangulateio *)s;
            }
        else { croak("THIS is not of type Math::Geometry::Delaunay::Triangulateio"); }
    OUTPUT:
        RETVAL

void *
triio_DESTROY(Math::Geometry::Delaunay::Triangulateio THIS = NO_INIT)
    PREINIT:
        struct triangulateio * p;
        STRLEN len;
        char *s;
    CODE:
        s = SvPV((SV*)SvRV(ST(0)), len);
        if (len != sizeof(THIS)) {croak("In DESTROY, size %d of packed data != expected %d", len, sizeof(THIS));}
        p = (struct triangulateio *) s;

        if (sv_derived_from(ST(0), "Math::Geometry::Delaunay::Triangulateio")) {
            if (len != sizeof(THIS)) {
                croak("Size %d of packed data != expected %d", len, sizeof(THIS));
                }
            if (p->pointlist)             {trifree(p->pointlist);}
            if (p->pointattributelist)    {trifree(p->pointattributelist);}
            if (p->pointmarkerlist)       {trifree(p->pointmarkerlist);}
            if (p->trianglelist)          {trifree(p->trianglelist);}
            if (p->triangleattributelist) {trifree(p->triangleattributelist);}
            if (p->trianglearealist)      {trifree(p->trianglearealist);}
            if (p->neighborlist)          {trifree(p->neighborlist);}
            if (p->segmentlist)           {trifree(p->segmentlist);}
            if (p->segmentmarkerlist)     {trifree(p->segmentmarkerlist);}
            if (p->holelist)              {trifree(p->holelist);}
            if (p->regionlist)            {trifree(p->regionlist);}
            if (p->edgelist)              {trifree(p->edgelist);}
            if (p->edgemarkerlist)        {trifree(p->edgemarkerlist);}
            if (p->normlist)              {trifree(p->normlist);}
            }
        else { croak("THIS is not of type Math::Geometry::Delaunay::Triangulateio"); }

