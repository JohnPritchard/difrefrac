/* "@(#) $Id: libdifrefrac.h 3021 2017-02-01 08:23:08Z jpritcha $" */
/*
  This is libdifrefrac.h

  Copyright (C) 2000  J.D.Pritchard <j.pritchard@eso.org>

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <libskycalc.h>

#define MAX(A,B) (((A)>(B)) ? (A) : (B))
#define MIN(A,B) (((A)<(B)) ? (A) : (B))

#ifdef __cplusplus
extern "C" {
#endif
void calcSafty( const char *telescope, double Dec, double safezone[2][2] );
void calcsssr( double mjd, double sssr[2], double lat, double longit, double elev );
double get_coorddes( const char *str );
void myput_date( struct date_time date, short p);
void myput_coords(double deci, short precision);
void sfmoon( double jds, double jdf, double moon[11], double lat, double longit, double elev );
void calcBaEofNight( double jdmid, double sssr[2], double moon[11], double lat, double longit, double elev, double altBEoN );
void calcminus13( double mjd, double sssr[2], double moon[11], double lat, double longit, double elev );
double hainm12top12 ( double jd, double longit, double RA);
double compPhase ( double jd, double T0, double P );
void printint( const char *text, double jdmid, double jdi, double jdf, double phi, double phf, double RA, double Dec, double lat, double longit, int ed, int preoh, int postoh, const char *label, double moon[11], short bPrintOpts );
void setTelescope( const char *telescope, short bptrec, short btrec, short *bCS, short *insz, double *longit, double *lat, double *elev, const char *EWtext[]);
#ifdef __cplusplus
}
#endif
