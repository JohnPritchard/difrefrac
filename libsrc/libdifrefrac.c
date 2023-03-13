/* "@(#) $Id: libdifrefrac.c 3021 2017-02-01 08:23:08Z jpritcha $" */
/*
  This is liblsosc.c

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

#include "libdifrefrac.h"

void
calcSafty( const char *telescope, double Dec, double safezone[2][2] )
{
  /*
    This subroutine returns safezone a 2x2 array which gives the safe HA ranges
    for the given Dec.

    The [i][0] refer to E
    The [i][1] refer to W
  */

  double Edata[53][3],Wdata[53][3];
  int i,iE,iW;

  iE=iW=0;

  if ( ! strcmp(telescope,"DK154")) {
    Edata[0][0]=-90.000;Edata[0][1]=6.000;Edata[0][2]=-0.450;
    Edata[1][0]=-87.565;Edata[1][1]=6.000;Edata[1][2]=-0.450;
    Edata[2][0]=-85.130;Edata[2][1]=6.000;Edata[2][2]=-0.450;
    Edata[3][0]=-82.695;Edata[3][1]=6.000;Edata[3][2]=-0.450;
    Edata[4][0]=-80.260;Edata[4][1]=6.000;Edata[4][2]=-0.450;
    Edata[5][0]=-77.825;Edata[5][1]=6.000;Edata[5][2]=-0.450;
    Edata[6][0]=-75.390;Edata[6][1]=6.000;Edata[6][2]=-0.450;
    Edata[7][0]=-72.955;Edata[7][1]=6.000;Edata[7][2]=-0.450;
    Edata[8][0]=-70.520;Edata[8][1]=6.000;Edata[8][2]=-0.450;
    Edata[9][0]=-68.085;Edata[9][1]=6.000;Edata[9][2]=-0.450;
    Edata[10][0]=-65.650;Edata[10][1]=6.000;Edata[10][2]=-0.450;
    Edata[11][0]=-63.215;Edata[11][1]=6.000;Edata[11][2]=-0.450;
    Edata[12][0]=-60.780;Edata[12][1]=6.000;Edata[12][2]=-0.450;
    Edata[13][0]=-58.345;Edata[13][1]=6.000;Edata[13][2]=-0.450;
    Edata[14][0]=-55.910;Edata[14][1]=6.000;Edata[14][2]=-0.450;
    Edata[15][0]=-53.475;Edata[15][1]=5.925;Edata[15][2]=-0.450;
    Edata[16][0]=-51.040;Edata[16][1]=5.925;Edata[16][2]=-0.525;
    Edata[17][0]=-48.605;Edata[17][1]=5.850;Edata[17][2]=-0.600;
    Edata[18][0]=-46.170;Edata[18][1]=5.850;Edata[18][2]=-0.750;
    Edata[19][0]=-43.735;Edata[19][1]=5.850;Edata[19][2]=-1.500;
    Edata[20][0]=-41.300;Edata[20][1]=5.775;Edata[20][2]=-3.750;
    Edata[21][0]=-38.865;Edata[21][1]=5.700;Edata[21][2]=-3.900;
    Edata[22][0]=-36.430;Edata[22][1]=5.700;Edata[22][2]=-4.050;
    Edata[23][0]=-33.995;Edata[23][1]=5.625;Edata[23][2]=-4.125;
    Edata[24][0]=-31.560;Edata[24][1]=5.550;Edata[24][2]=-4.125;
    Edata[25][0]=-29.125;Edata[25][1]=5.550;Edata[25][2]=-4.050;
    Edata[26][0]=-26.690;Edata[26][1]=5.475;Edata[26][2]=-3.975;
    Edata[27][0]=-24.255;Edata[27][1]=5.400;Edata[27][2]=-3.900;
    Edata[28][0]=-21.820;Edata[28][1]=5.250;Edata[28][2]=-3.900;
    Edata[29][0]=-19.385;Edata[29][1]=5.175;Edata[29][2]=-3.750;
    Edata[30][0]=-16.950;Edata[30][1]=5.100;Edata[30][2]=-3.750;
    Edata[31][0]=-14.515;Edata[31][1]=5.025;Edata[31][2]=-3.675;
    Edata[32][0]=-12.080;Edata[32][1]=4.875;Edata[32][2]=-3.600;
    Edata[33][0]=-9.645;Edata[33][1]=4.800;Edata[33][2]=-3.525;
    Edata[34][0]=-7.210;Edata[34][1]=4.650;Edata[34][2]=-3.450;
    Edata[35][0]=-4.775;Edata[35][1]=4.500;Edata[35][2]=-3.300;
    Edata[36][0]=-2.340;Edata[36][1]=4.350;Edata[36][2]=-3.225;
    Edata[37][0]=0.095;Edata[37][1]=4.200;Edata[37][2]=-3.150;
    Edata[38][0]=2.530;Edata[38][1]=4.050;Edata[38][2]=-3.075;
    Edata[39][0]=4.965;Edata[39][1]=3.975;Edata[39][2]=-3.000;
    Edata[40][0]=7.400;Edata[40][1]=3.900;Edata[40][2]=-2.925;
    Edata[41][0]=9.835;Edata[41][1]=3.750;Edata[41][2]=-2.850;
    Edata[42][0]=12.270;Edata[42][1]=3.600;Edata[42][2]=-2.775;
    Edata[43][0]=14.705;Edata[43][1]=3.450;Edata[43][2]=-2.700;
    Edata[44][0]=17.140;Edata[44][1]=3.300;Edata[44][2]=-2.625;
    // New limit imposed by Cable-Protection-System
    Edata[45][0]=18.000;Edata[45][1]=3.250;Edata[45][2]=-2.600;
    /* The old limits...
    Edata[45][0]=19.575;Edata[45][1]=3.150;Edata[45][2]=-2.550;
    Edata[46][0]=22.010;Edata[46][1]=2.850;Edata[46][2]=-2.400;
    Edata[47][0]=24.445;Edata[47][1]=2.700;Edata[47][2]=-2.250;
    Edata[48][0]=26.880;Edata[48][1]=2.250;Edata[48][2]=-2.100;
    Edata[49][0]=29.315;Edata[49][1]=2.100;Edata[49][2]=-1.800;
    Edata[50][0]=31.750;Edata[50][1]=1.500;Edata[50][2]=-1.200;
    Edata[51][0]=34.185;Edata[51][1]=1.050;Edata[51][2]=-0.750;
    */
    iE=45;

    Wdata[0][0]=-90.000;Wdata[0][1]=0.400;Wdata[0][2]=-6.050;
    Wdata[1][0]=-87.600;Wdata[1][1]=0.400;Wdata[1][2]=-6.050;
    Wdata[2][0]=-85.200;Wdata[2][1]=0.400;Wdata[2][2]=-6.050;
    Wdata[3][0]=-82.800;Wdata[3][1]=0.400;Wdata[3][2]=-6.050;
    Wdata[4][0]=-80.400;Wdata[4][1]=0.400;Wdata[4][2]=-6.050;
    Wdata[5][0]=-78.000;Wdata[5][1]=0.400;Wdata[5][2]=-6.050;
    Wdata[6][0]=-75.600;Wdata[6][1]=0.400;Wdata[6][2]=-6.050;
    Wdata[7][0]=-73.200;Wdata[7][1]=0.400;Wdata[7][2]=-6.050;
    Wdata[8][0]=-70.800;Wdata[8][1]=0.400;Wdata[8][2]=-6.050;
    Wdata[9][0]=-68.400;Wdata[9][1]=0.400;Wdata[9][2]=-6.050;
    Wdata[10][0]=-66.000;Wdata[10][1]=0.400;Wdata[10][2]=-6.050;
    Wdata[11][0]=-63.600;Wdata[11][1]=0.475;Wdata[11][2]=-6.050;
    Wdata[12][0]=-61.200;Wdata[12][1]=0.520;Wdata[12][2]=-6.050;
    Wdata[13][0]=-58.800;Wdata[13][1]=0.550;Wdata[13][2]=-6.050;
    Wdata[14][0]=-56.400;Wdata[14][1]=0.625;Wdata[14][2]=-6.050;
    Wdata[15][0]=-54.000;Wdata[15][1]=0.700;Wdata[15][2]=-6.050;
    Wdata[16][0]=-51.600;Wdata[16][1]=0.850;Wdata[16][2]=-5.975;
    Wdata[17][0]=-49.200;Wdata[17][1]=3.700;Wdata[17][2]=-5.975;
    Wdata[18][0]=-46.800;Wdata[18][1]=4.300;Wdata[18][2]=-5.900;
    Wdata[19][0]=-44.400;Wdata[19][1]=4.450;Wdata[19][2]=-5.825;
    Wdata[20][0]=-42.000;Wdata[20][1]=4.555;Wdata[20][2]=-5.750;
    Wdata[21][0]=-39.600;Wdata[21][1]=4.600;Wdata[21][2]=-5.675;
    Wdata[22][0]=-37.200;Wdata[22][1]=4.600;Wdata[22][2]=-5.525;
    Wdata[23][0]=-34.800;Wdata[23][1]=4.600;Wdata[23][2]=-5.450;
    Wdata[24][0]=-32.400;Wdata[24][1]=4.600;Wdata[24][2]=-5.375;
    Wdata[25][0]=-30.000;Wdata[25][1]=4.450;Wdata[25][2]=-5.300;
    Wdata[26][0]=-27.600;Wdata[26][1]=4.450;Wdata[26][2]=-5.225;
    Wdata[27][0]=-25.200;Wdata[27][1]=4.375;Wdata[27][2]=-5.150;
    Wdata[28][0]=-22.800;Wdata[28][1]=4.150;Wdata[28][2]=-5.075;
    Wdata[29][0]=-20.400;Wdata[29][1]=4.075;Wdata[29][2]=-5.000;
    Wdata[30][0]=-18.000;Wdata[30][1]=4.000;Wdata[30][2]=-4.925;
    Wdata[31][0]=-15.600;Wdata[31][1]=3.850;Wdata[31][2]=-4.850;
    Wdata[32][0]=-13.200;Wdata[32][1]=3.700;Wdata[32][2]=-4.775;
    Wdata[33][0]=-10.800;Wdata[33][1]=3.550;Wdata[33][2]=-4.700;
    Wdata[34][0]=-8.400;Wdata[34][1]=3.400;Wdata[34][2]=-4.625;
    Wdata[35][0]=-6.000;Wdata[35][1]=3.250;Wdata[35][2]=-4.475;
    Wdata[36][0]=-3.600;Wdata[36][1]=3.025;Wdata[36][2]=-4.400;
    Wdata[37][0]=-1.200;Wdata[37][1]=2.800;Wdata[37][2]=-4.250;
    Wdata[38][0]=1.200;Wdata[38][1]=2.650;Wdata[38][2]=-4.100;
    Wdata[39][0]=3.600;Wdata[39][1]=2.500;Wdata[39][2]=-4.025;
    Wdata[40][0]=6.000;Wdata[40][1]=2.350;Wdata[40][2]=-3.950;
    Wdata[41][0]=8.400;Wdata[41][1]=2.275;Wdata[41][2]=-3.800;
    Wdata[42][0]=10.800;Wdata[42][1]=2.200;Wdata[42][2]=-3.650;
    Wdata[43][0]=13.200;Wdata[43][1]=2.050;Wdata[43][2]=-3.500;
    Wdata[44][0]=15.600;Wdata[44][1]=1.975;Wdata[44][2]=-3.425;
    Wdata[45][0]=18.000;Wdata[45][1]=1.900;Wdata[45][2]=-3.200;
    Wdata[46][0]=20.400;Wdata[46][1]=1.825;Wdata[46][2]=-2.975;
    Wdata[47][0]=22.800;Wdata[47][1]=1.750;Wdata[47][2]=-2.750;
    Wdata[48][0]=25.200;Wdata[48][1]=1.525;Wdata[48][2]=-2.450;
    Wdata[49][0]=27.600;Wdata[49][1]=1.300;Wdata[49][2]=-2.000;
    // New limit imposed by Cable-Protection-System
    Wdata[50][0]=28.000;Wdata[50][1]=1.250;Wdata[50][2]=-1.950;
    /* The old limits...
    Wdata[50][0]=30.000;Wdata[50][1]=1.000;Wdata[50][2]=-1.700;
    Wdata[51][0]=32.400;Wdata[51][1]=0.850;Wdata[51][2]=-1.250;
    Wdata[52][0]=34.800;Wdata[52][1]=0.250;Wdata[52][2]=-0.650;
    */
    iW=50;
  }

  if ( ! strcmp(telescope,"ESO152")) {
    Edata[0][0]=-90.000;Edata[0][1]=6.000;Edata[0][2]=-0.230;
    Edata[1][0]=-46.590;Edata[1][1]=6.000;Edata[1][2]=-0.430;
    Edata[2][0]=-40.900;Edata[2][1]=6.000;Edata[2][2]=-0.930;
    Edata[3][0]=-40.900;Edata[3][1]=6.000;Edata[3][2]=-1.500;
    Edata[4][0]=-40.000;Edata[4][1]=6.000;Edata[4][2]=-1.580;
    Edata[5][0]=-40.000;Edata[5][1]=6.000;Edata[5][2]=-5.380;
    Edata[6][0]=-05.680;Edata[6][1]=6.000;Edata[6][2]=-4.450;
    Edata[7][0]=-04.540;Edata[7][1]=6.000;Edata[7][2]=-4.370;
    Edata[8][0]=-04.540;Edata[8][1]=5.640;Edata[8][2]=-4.370;
    Edata[9][0]=+23.860;Edata[9][1]=4.660;Edata[9][2]=-3.090;
    Edata[10][0]=+23.860;Edata[10][1]=4.660;Edata[10][2]=-0.830;
    Edata[11][0]=+46.020;Edata[11][1]=3.000;Edata[11][2]=-0.510;
    Edata[12][0]=+46.020;Edata[12][1]=3.000;Edata[12][2]=+0.000;
    Edata[13][0]=+50.000;Edata[13][1]=1.660;Edata[13][2]=+0.390;
    iE=13;

    Wdata[0][0]=-90.000;Wdata[0][1]=0.250;Wdata[0][2]=-6.000;
    Wdata[1][0]=-42.040;Wdata[1][1]=0.450;Wdata[1][2]=-6.000;
    Wdata[2][0]=-42.040;Wdata[2][1]=5.260;Wdata[2][2]=-6.000;
    Wdata[3][0]=-10.000;Wdata[3][1]=4.180;Wdata[3][2]=-6.000;
    Wdata[4][0]=-02.270;Wdata[4][1]=3.720;Wdata[4][2]=-6.000;
    Wdata[5][0]=+00.000;Wdata[5][1]=3.500;Wdata[5][2]=-5.640;
    Wdata[6][0]=+10.000;Wdata[6][1]=3.000;Wdata[6][2]=-5.000;
    Wdata[7][0]=+20.000;Wdata[7][1]=2.510;Wdata[7][2]=-4.560;
    Wdata[8][0]=+20.000;Wdata[8][1]=0.690;Wdata[8][2]=-4.560;
    Wdata[9][0]=+23.860;Wdata[9][1]=0.390;Wdata[9][2]=-4.460;
    Wdata[10][0]=+26.130;Wdata[10][1]=0.390;Wdata[10][2]=-4.400;
    Wdata[11][0]=+44.310;Wdata[11][1]=0.310;Wdata[11][2]=-2.250;
    Wdata[12][0]=+44.310;Wdata[12][1]=0.000;Wdata[12][2]=-2.250;
    Wdata[13][0]=+50.000;Wdata[13][1]=-0.470;Wdata[13][2]=-1.530;
    iW=13;
  }

  safezone[0][0]=-99.;
  safezone[1][0]=-99.;
  safezone[0][1]=-99.;
  safezone[1][1]=-99.;

  i=0;
  while (( Edata[i][0] < Dec ) && ( i < (iE+1) )) i++;
  if (( i != 0 ) && ( i != (iE+1) )) {
    safezone[0][0]=Edata[i-1][2]+(Dec-Edata[i-1][0])/(Edata[i][0]-Edata[i-1][0])*(Edata[i][2]-Edata[i-1][2]);
    safezone[1][0]=Edata[i-1][1]+(Dec-Edata[i-1][0])/(Edata[i][0]-Edata[i-1][0])*(Edata[i][1]-Edata[i-1][1]);
  }
  i=0;
  while (( Wdata[i][0] < Dec ) && ( i < (iW+1) )) i++;
  if (( i != 0 ) && ( i != (iW+1) )) {
    safezone[0][1]=Wdata[i-1][2]+(Dec-Wdata[i-1][0])/(Wdata[i][0]-Wdata[i-1][0])*(Wdata[i][2]-Wdata[i-1][2]);
    safezone[1][1]=Wdata[i-1][1]+(Dec-Wdata[i-1][0])/(Wdata[i][0]-Wdata[i-1][0])*(Wdata[i][1]-Wdata[i-1][1]);
  }
  return;

  /*
    From DSZC, by J. Brewer, April 1997, provided in file DSZC.tar.gz by...

    An email from Heath Jones, 2p2 Team member...
    =============================================
    Date: Fri, 13 Oct 2000 15:59:37 -0400
    From: "2.2 telescope team" <2p2team@eso.org>
    Organization: European Southern Observatory
    To: John Pritchard <j.pritchard@astro.ku.dk>
    CC: Patrick FRANCOIS <fpatrick@eso.org>, James Brewer <jbrewer@astro.ubc.ca>
    Subject: Re D1.54 Safety Zone map

    **) E.data
#     Dec      Min       Max
#
  -90.000    6.000   -0.450
  -87.565    6.000   -0.450
  -85.130    6.000   -0.450
  -82.695    6.000   -0.450
  -80.260    6.000   -0.450
  -77.825    6.000   -0.450
  -75.390    6.000   -0.450
  -72.955    6.000   -0.450
  -70.520    6.000   -0.450
  -68.085    6.000   -0.450
  -65.650    6.000   -0.450
  -63.215    6.000   -0.450
  -60.780    6.000   -0.450
  -58.345    6.000   -0.450
  -55.910    6.000   -0.450
  -53.475    5.925   -0.450
  -51.040    5.925   -0.525
  -48.605    5.850   -0.600
  -46.170    5.850   -0.750
  -43.735    5.850   -1.500
  -41.300    5.775   -3.750
  -38.865    5.700   -3.900
  -36.430    5.700   -4.050
  -33.995    5.625   -4.125
  -31.560    5.550   -4.125
  -29.125    5.550   -4.050
  -26.690    5.475   -3.975
  -24.255    5.400   -3.900
  -21.820    5.250   -3.900
  -19.385    5.175   -3.750
  -16.950    5.100   -3.750
  -14.515    5.025   -3.675
  -12.080    4.875   -3.600
   -9.645    4.800   -3.525
   -7.210    4.650   -3.450
   -4.775    4.500   -3.300
   -2.340    4.350   -3.225
    0.095    4.200   -3.150
    2.530    4.050   -3.075
    4.965    3.975   -3.000
    7.400    3.900   -2.925
    9.835    3.750   -2.850
   12.270    3.600   -2.775
   14.705    3.450   -2.700
   17.140    3.300   -2.625
   19.575    3.150   -2.550
   22.010    2.850   -2.400
   24.445    2.700   -2.250
   26.880    2.250   -2.100
   29.315    2.100   -1.800
   31.750    1.500   -1.200
   34.185    1.050   -0.750

    **) W.data
#     Dec      Min       Max
#
  -90.000    0.400   -6.050
  -87.600    0.400   -6.050
  -85.200    0.400   -6.050
  -82.800    0.400   -6.050
  -80.400    0.400   -6.050
  -78.000    0.400   -6.050
  -75.600    0.400   -6.050
  -73.200    0.400   -6.050
  -70.800    0.400   -6.050
  -68.400    0.400   -6.050
  -66.000    0.400   -6.050
  -63.600    0.475   -6.050
  -61.200    0.520   -6.050
  -58.800    0.550   -6.050
  -56.400    0.625   -6.050
  -54.000    0.700   -6.050
  -51.600    0.850   -5.975
  -49.200    3.700   -5.975
  -46.800    4.300   -5.900
  -44.400    4.450   -5.825
  -42.000    4.555   -5.750
  -39.600    4.600   -5.675
  -37.200    4.600   -5.525
  -34.800    4.600   -5.450
  -32.400    4.600   -5.375
  -30.000    4.450   -5.300
  -27.600    4.450   -5.225
  -25.200    4.375   -5.150
  -22.800    4.150   -5.075
  -20.400    4.075   -5.000
  -18.000    4.000   -4.925
  -15.600    3.850   -4.850
  -13.200    3.700   -4.775
  -10.800    3.550   -4.700
   -8.400    3.400   -4.625
   -6.000    3.250   -4.475
   -3.600    3.025   -4.400
   -1.200    2.800   -4.250
    1.200    2.650   -4.100
    3.600    2.500   -4.025
    6.000    2.350   -3.950
    8.400    2.275   -3.800
   10.800    2.200   -3.650
   13.200    2.050   -3.500
   15.600    1.975   -3.425
   18.000    1.900   -3.200
   20.400    1.825   -2.975
   22.800    1.750   -2.750
   25.200    1.525   -2.450
   27.600    1.300   -2.000
   30.000    1.000   -1.700
   32.400    0.850   -1.250
   34.800    0.250   -0.650

  */
}

void calcsssr( double mjd, double sssr[2], double lat, double longit, double elev )
{
  int y;

  /* *****************************************************
     From skycalendar.c...
			John Thorstensen
			Dept. of Physics and Astronomy
			Dartmouth College
			Hanover, NH 03755
			John.Thorstensen@dartmouth.edu

  */
  double jd, jdmid, jdeve, jdmorn, jdsunset, jdsunrise;
  double jdmoonrise,jdmoonset,jdetwilight,jdmtwilight,jd1,jd2,jdout;
  double rasun, decsun;
  double jdbdst, jdedst;  /* jd at beg and end of dst this year */
  double stmid;
  double hamoonset, hasunset, hatwilight, tmoonset, tmoonrise;

  short year;

  double stdz = 4.;
  short use_dst = -1;
  double moon_pr = 7.5;
  double horiz;

  struct date_time date;
  short  dow;
  caldat((mjd+2400000.5),&date,&dow);
  year=(short)date.y;


  /* *****************************************************
     Hacked from skycalendar.c...
  */
  horiz = pow((2. * elev / 6378140.),0.5) * DEG_IN_RADIAN;
  if(use_dst != 0) find_dst_bounds(year,stdz,use_dst,&jdbdst,&jdedst);
  jd=mjd;
  jdmid = 2400000.5 + mjd + zone(use_dst,stdz,jd,jdbdst,jdedst) / 24.;  /* corresponding ut */
  stmid = lst(jdmid,longit);

  lpsun(jdmid,&rasun,&decsun);
  hasunset = ha_alt(decsun,lat,-(0.83+horiz));
/* there follows some awful flow of control.  This is an artifact of the
   history of the program - it started as a program for mid-latitudes, and
   then I changed it to handle polar latitudes reasonably well.  This led
   to a lot of really weird branching.  Sorry.  It appears to work, though. */
  if(hasunset > 900.) {  /* sun isn't going to set; no twilight */
    moon_pr = 10.;   /* print moonrise/set if within 10 hr of midn. */
    goto MIDNIGHT_SUN; /* horrible flow of control */
  }
  if(hasunset < -900.) {  /* sun ain't gonna rise, but may be twilight */
    moon_pr = 10.;
    goto TWILIGHT;    /* more horrible flow of control */
  }
  jdsunset = jdmid + adj_time(rasun+hasunset-stmid)/24.; /* initial guess */
  jdsunset = jd_sun_alt(-(0.83+horiz),
			jdsunset,lat,longit); /* refinement */
  if(jdsunset > 0.) {/* convergence */
    sssr[0]=jdsunset;
  } else printf(" .....");
 TWILIGHT:
 SUNRISE:
  if(hasunset < -900.) printf(" .....");
  else {
    jdsunrise = jdmid + adj_time(rasun-hasunset-stmid)/24.;
    jdsunrise = jd_sun_alt(-(0.83+horiz),jdsunrise,lat,longit);
    if(jdsunrise > 0.) {
      sssr[1]=jdsunrise;
    } else printf(" .....");
  }
 MIDNIGHT_SUN:
  return;

}

double get_coorddes( const char *str )

/* Based on get_coord()
 */

/* Reads a string from the terminal and converts it into
   a double-precision coordinate.  This is trickier than
   it appeared at first, since a -00 tests as non-negative;
   the sign has to be picked out and handled explicitly. */
/* Prompt for input in the calling routine.*/
{
   short sign;
   double hrs, mins, secs;
   char hh_string[20];  /* string with the first coord (hh) */
   char hh1[1];
   short i = 0;

   if ( str == NULL ) {
     printf("str is NULL\n");
     return( -99. );
   }

   /* read and handle the hour (or degree) part with sign */

   sscanf(str,"%s",hh_string);
   while ( i < 20 ) {
     if ( hh_string[i] == ':' ) hh_string[i]=' ';
     i++;
   }

   i=0;
   hh1[0] = hh_string[i];
   while(hh1[0] == ' ') {
       /* discard leading blanks */
       i++;
       hh1[0] = hh_string[i];
   }
   sign=1;
   if(hh1[0] == '-') sign = -1;

   /* read in the hours, minutes and seconds normally */
   sscanf(hh_string, "%lf %lf %lf",&hrs,&mins,&secs);
   if(sign == -1) hrs = -1. * hrs;

   return( (double) sign * (hrs + mins / 60. + secs / 3600.));
}


void myput_date( struct date_time date, short p)
{
/* prints out a struct coord in a nice format; precision
   is a code for how accurate you want it.  The options are:
     precision = 0;   minutes rounded to the nearest minute
     precision = 1;   minutes rounded to the nearest tenth.
     precision = 2;   seconds rounded to the nearest second
     precision = 3;   seconds given to the tenth
     precision = 4;   seconds given to the hundredth
   The program assumes that the line is ready for the coord
   to be printed and does NOT deliver a new line at the end
   of the output. */

  int ihh,imm,iss;
  double dhh,dmm,dss;
  double mht;

  ihh=date.h;imm=date.mn,dss=date.s;
  dhh=(double)ihh;iss=(int)dss;

  mht=0.5;
  if ( p == 1 ) mht=0.95;
  if ( p <= 1 ) {
    dmm=(double)imm + dss/60.;
    if ( (dmm - (double)imm) > mht ) {
      dss=0.;
      imm++;
      dmm=(double)imm + dss/60.;
    }
  }
  printf("%2.2d:",ihh);
  if ( p == 1 ) {
    printf("%04.1f",dmm);
  } else {
    printf("%2.2d",imm);
  }

  if ( p > 1 ) {
    if ( p == 2 ) printf(":%2.2d",((int)(dss+0.5)));
    if ( p == 3 ) printf(":%04.1f",dss);
    if ( p == 4 ) printf(":%05.2f",dss);
  }

  return;

}

void myput_coords(double deci, short precision)


/* prints out a struct coord in a nice format; precision
   is a code for how accurate you want it.  The options are:
     precision = 0;   minutes rounded to the nearest minute
     precision = 1;   minutes rounded to the nearest tenth.
     precision = 2;   seconds rounded to the nearest second
     precision = 3;   seconds given to the tenth
     precision = 4;   seconds given to the hundredth
   The program assumes that the line is ready for the coord
   to be printed and does NOT deliver a new line at the end
   of the output. */

{

   double minutes;  /* for rounding off if necess. */
   struct coord out_coord, coords;
   char out_string[20];  /* for checking for nasty 60's */

   dec_to_bab(deci,&coords);  /* internally convert to coords*/

   if(precision == 0) {   /* round to nearest minute */
      minutes = coords.mm + coords.ss / 60.;
	   /* check to be sure minutes aren't 60 */
      sprintf(out_string,"%.0f %02.0f",coords.hh,minutes);
      sscanf(out_string,"%lf %lf",&out_coord.hh,&out_coord.mm);
      if(fabs(out_coord.mm - 60.) < 1.0e-7) {
	 out_coord.mm = 0.;
	 out_coord.hh = out_coord.hh + 1.;
      }
      /* if(out_coord.hh < 100.) oprntf(" "); */  /* put in leading blanks explicitly
         for 'h' option below. */
      //if(out_coord.hh < 10.) oprntf(" ");
      if(coords.sign == -1) oprntf("-");
	else oprntf(" ");   /* preserves alignment */
      oprntf("%02.0f:%02.0f",out_coord.hh,out_coord.mm);
   }

   else if(precision == 1) {    /* keep nearest tenth of a minute */
      minutes = coords.mm + coords.ss / 60.;
	   /* check to be sure minutes are not 60 */
      sprintf(out_string,"%.0f %04.1f",coords.hh,minutes);
      sscanf(out_string,"%lf %lf",&out_coord.hh, &out_coord.mm);
      if(fabs(out_coord.mm - 60.) < 1.0e-7) {
	 out_coord.mm = 0.;
	 out_coord.hh = out_coord.hh + 1.;
      }
      //if(out_coord.hh < 10.) oprntf(" ");
      if(coords.sign == -1) oprntf("-");
	else oprntf(" ");   /* preserves alignment */
      oprntf("%02.0f:%04.1f", out_coord.hh, out_coord.mm);
   }
   else if(precision == 2) {
	  /* check to be sure seconds are not 60 */
      sprintf(out_string,"%.0f %02.0f %02.0f",coords.hh,coords.mm,coords.ss);
      sscanf(out_string,"%lf %lf %lf",&out_coord.hh,&out_coord.mm,
	   &out_coord.ss);
      if(fabs(out_coord.ss - 60.) < 1.0e-7) {
	  out_coord.mm = out_coord.mm + 1.;
	  out_coord.ss = 0.;
	  if(fabs(out_coord.mm - 60.) < 1.0e-7) {
	      out_coord.hh = out_coord.hh + 1.;
	      out_coord.mm = 0.;
	  }
      }
      //if(out_coord.hh < 10.) oprntf(" ");
      if(coords.sign == -1) oprntf("-");
	 else oprntf(" ");   /* preserves alignment */
      oprntf("%02.0f:%02.0f:%02.0f",out_coord.hh,out_coord.mm,out_coord.ss);
   }
   else if(precision == 3) {
	  /* the usual shuffle to check for 60's */
      sprintf(out_string,"%.0f %02.0f %04.1f",coords.hh, coords.mm, coords.ss);
      sscanf(out_string,"%lf %lf %lf",&out_coord.hh,&out_coord.mm,
	   &out_coord.ss);
      if(fabs(out_coord.ss - 60.) < 1.0e-7) {
	  out_coord.mm = out_coord.mm + 1.;
	  out_coord.ss = 0.;
	  if(fabs(out_coord.mm - 60.) < 1.0e-7) {
	     out_coord.hh = out_coord.hh + 1.;
	     out_coord.mm = 0.;
	  }
      }
      //if(out_coord.hh < 10.) oprntf(" ");
      if(coords.sign == -1) oprntf("-");
	 else oprntf(" ");   /* preserves alignment */
      oprntf("%02.0f:%02.0f:%04.1f",out_coord.hh,out_coord.mm,out_coord.ss);
   }
   else {
      sprintf(out_string,"%.0f %02.0f %05.2f",coords.hh,coords.mm,coords.ss);
      sscanf(out_string,"%lf %lf %lf",&out_coord.hh,&out_coord.mm,
	   &out_coord.ss);
      if(fabs(out_coord.ss - 60.) < 1.0e-6) {
	 out_coord.mm = out_coord.mm + 1.;
	 out_coord.ss = 0.;
	 if(fabs(out_coord.mm - 60.) < 1.0e-6) {
	    out_coord.hh = out_coord.hh + 1.;
	    out_coord.mm = 0.;
	 }
      }
      //if(out_coord.hh < 10.) oprntf(" ");
      if(coords.sign == -1) oprntf("-");
	 else oprntf(" ");   /* preserves alignment */
      oprntf("%02.0f:%02.0f:%05.2f",out_coord.hh, out_coord.mm, out_coord.ss);
   }
}

void sfmoon( double jds, double jdf, double moon[11], double lat, double longit, double elev )
{
  double geora, geodec, geodist;  /* geocent for moon, not used here.*/
  double rasun, decsun;
  double stmid, ramoon, decmoon, distmoon;

  stmid = lst((((int) jds) + 0.5 + longit/24. ),longit);

  accumoon(jds,lat,stmid,elev,
           &geora,&geodec,&geodist,&ramoon,&decmoon,&distmoon);
  lpsun(jds,&rasun,&decsun);
  moon[6]=ramoon;
  moon[7]=decmoon;
  moon[8]=0.5*(1.-cos(subtend(ramoon,decmoon,rasun,decsun)));
  if ( jdf > 0. ) {
    accumoon(jdf,lat,stmid,elev,
             &geora,&geodec,&geodist,&ramoon,&decmoon,&distmoon);
    lpsun(jdf,&rasun,&decsun);
    moon[9]=ramoon;
    moon[10]=decmoon;
    moon[11]=0.5*(1.-cos(subtend(ramoon,decmoon,rasun,decsun)));
  }

}

void calcBaEofNight( double jdmid, double sssr[2], double moon[11], double lat, double longit, double elev, double altBEoN )
{

  /* *****************************************************
     From skycalendar.c...
			John Thorstensen
			Dept. of Physics and Astronomy
			Dartmouth College
			Hanover, NH 03755
			John.Thorstensen@dartmouth.edu

  */
  double jd, jdeve, jdmorn, jdsunset, jdsunrise;
  double jdmoonrise,jdmoonset,jdetwilight,jdmtwilight,jd1,jd2,jdout;
  double geora, geodec, geodist;  /* geocent for moon, not used here.*/
  double rasun, decsun;
  double jdbdst, jdedst;  /* jd at beg and end of dst this year */
  double hamoonset, hasunset, hatwilight, tmoonset, tmoonrise;
  double stmid, ramoon, decmoon, distmoon;

  short year;

  double stdz = 4.;
  short use_dst = -1;
  double moon_pr = 7.5;
  double horiz;

  struct date_time date;
  short  dow;
  caldat(jdmid,&date,&dow);
  year=(short)date.y;

  jdsunset=jdsunrise=-1.;

  /* *****************************************************
     Hacked from skycalendar.c...
  */
  //horiz = pow((2. * elev / 6378140.),0.5) * DEG_IN_RADIAN;
  horiz = sqrt(2. * elev / 6378140.) * DEG_IN_RADIAN;
  if(use_dst != 0) find_dst_bounds(year,stdz,use_dst,&jdbdst,&jdedst);
  //jdmid = 2400000.5 + mjd + zone(use_dst,stdz,jd,jdbdst,jdedst) / 24.;  /* corresponding ut */
  stmid = lst(jdmid,longit);

  accumoon(jdmid,lat,stmid,elev,
           &geora,&geodec,&geodist,&ramoon,&decmoon,&distmoon);
  lpsun(jdmid,&rasun,&decsun);
  hasunset = ha_alt(decsun,lat,(altBEoN-horiz));
/* there follows some awful flow of control.  This is an artifact of the
   history of the program - it started as a program for mid-latitudes, and
   then I changed it to handle polar latitudes reasonably well.  This led
   to a lot of really weird branching.  Sorry.  It appears to work, though. */
  if(hasunset > 900.) {  /* sun isn't going to set; no twilight */
    moon_pr = 10.;   /* print moonrise/set if within 10 hr of midn. */
    goto MIDNIGHT_SUN; /* horrible flow of control */
  }
  if(hasunset < -900.) {  /* sun ain't gonna rise, but may be twilight */
    moon_pr = 10.;
    goto TWILIGHT;    /* more horrible flow of control */
  }
  jdsunset = jdmid + adj_time(rasun+hasunset-stmid)/24.; /* initial guess */
  jdsunset = jd_sun_alt((altBEoN-horiz),
			jdsunset,lat,longit); /* refinement */
  sssr[0]=jdsunset;
 TWILIGHT:
 SUNRISE:
  if(hasunset < -900.) {
    printf(" .....");
  } else {
    jdsunrise = jdmid + adj_time(rasun-hasunset-stmid)/24.;
    jdsunrise = jd_sun_alt((altBEoN-horiz),jdsunrise,lat,longit);
    sssr[1]=jdsunrise;
  }
 MIDNIGHT_SUN:
  moon[1]=moon[2]=moon[3]=moon[4]=moon[5]=moon[6]=moon[7]=moon[8]=moon[9]=moon[10]=moon[11]=-1.;
  hamoonset = ha_alt(decmoon,lat,(0.-horiz)); /* rough approx. */
  if(fabs(hamoonset) < 100.) { /* if it's really gonna set */
    tmoonrise = adj_time(ramoon-hamoonset-stmid);
    tmoonset = adj_time(ramoon+hamoonset-stmid);
    jdmoonrise = jdmid + tmoonrise / 24.;
    jdmoonrise = jd_moon_alt((0.-horiz),jdmoonrise,lat,longit,elev);
    moon[1]=jdmoonrise;
    jdmoonset = jdmid + tmoonset / 24.;
    if ( jdmoonset < jdmoonrise ) jdmoonset=jdmoonset+1.;
    jdmoonset = jd_moon_alt((0.-horiz),jdmoonset,lat,longit,elev);
    moon[2]=jdmoonset;
  }
  moon[3]=ramoon;
  moon[4]=decmoon;
  moon[5]=0.5*(1.-cos(subtend(ramoon,decmoon,rasun,decsun)));
  accumoon(jdsunset,lat,stmid,elev,
           &geora,&geodec,&geodist,&ramoon,&decmoon,&distmoon);
  lpsun(jdsunset,&rasun,&decsun);
  moon[6]=ramoon;
  moon[7]=decmoon;
  moon[8]=0.5*(1.-cos(subtend(ramoon,decmoon,rasun,decsun)));
  accumoon(jdsunrise,lat,stmid,elev,
           &geora,&geodec,&geodist,&ramoon,&decmoon,&distmoon);
  lpsun(jdsunrise,&rasun,&decsun);
  moon[9]=ramoon;
  moon[10]=decmoon;
  moon[11]=0.5*(1.-cos(subtend(ramoon,decmoon,rasun,decsun)));

  return;

}

void calcminus13( double mjd, double sssr[2], double moon[11], double lat, double longit, double elev )
{
  calcBaEofNight(  mjd, sssr, moon, lat, longit, elev, -13.0 );
  return;
}

double hainm12top12 ( double jd, double longit, double RA)
{
  double ha=lst(jd,longit)-RA;
  if ( ha > 12. ) ha=ha-24.;
  if ( ha < -12. ) ha=ha+24.;
  return( ha );
}

double compPhase ( double jd, double T0, double P )
{
  double ph=-1.;
  if ( P >= 0. ) {
    ph=((jd-T0)/P - ((double) ((int) ((jd-T0)/P))));
    if ( ph <= 0. ) ph=1.+ph;
  }
  return ( ph );
}

void printint( const char *text, double jdmid, double jdi, double jdf, double phi, double phf, double RA, double Dec, double lat, double longit, int ed, int preoh, int postoh, const char *label, double moon[11], short bPrintOpts )
{
  short  dow;
  struct date_time date;
  double az,mjd;
  int i,hr,min,emin,esec;

  if ( text != NULL ) {
    if (( jdf < 0. ) || ( jdf > jdi )) {
      printf("%s ",text);

      // The Local date of the night
      caldat( (jdmid-0.5), &date, &dow);
      printf("%4.4d-%2.2d-%2.2d/",date.y,date.mo,date.d);
      caldat( (jdmid+0.5), &date, &dow);
      printf("%2.2d %9.1f",date.d,jdmid);

      // The Moon
      if ( moon[0] > 0. ) printf("  %4.2f ",moon[5]);
      if ( moon[0] > 1. ) {
        for ( i=1;i<3;i++) {
          if ( ( moon[i] > jdi ) && ( moon[i] < jdf ) ) {
            caldat( moon[i], &date, &dow);
            myput_date(date,0); printf(" ");
          } else {
            printf("--:-- ");
          }
        }
      }

      caldat( jdi, &date, &dow);
      printf(" : JD[s]=%13.5f %4.4d-%2.2d-%2.2d %2.2d:%2.2d:%02.0f [%2.2d:%2.2d",jdi,date.y,date.mo,date.d,date.h,date.mn,date.s,(int)lst(jdi,longit),((int)((lst(jdi,longit)-(int)lst(jdi,longit))*60.+0.5)));
      if ( RA >=0. ) {
	printf(",");
	myput_coords(hainm12top12(jdi,longit,RA),0);
	printf(",%6.2f",secant_z(altit(Dec,hainm12top12(jdi,longit,RA),lat,&az)));
        if ( bPrintOpts%2 > 0 ) printf(",%6.1f",parang((lst(jdi,longit)-RA),Dec,lat));
        if ( bPrintOpts%4 > 1 ) printf(",%6.1f,%6.1f",altit(Dec,(lst(jdi,longit)-RA),lat,&az),az);
      }
      printf("] ");
      if ( phi >= 0. ) printf("%6.4f ",phi);
      if ( moon[0] > 0. ) printf(" %4.2f  ",moon[8]);
      if (( moon[0] > 0. ) && ( RA >=0. )) printf("%5.1f  ",(180.*subtend(RA,Dec,moon[6],moon[7])/PI));
      // The finish date:time
      if ( jdf > 0. ) {
	caldat( jdf, &date, &dow );
	printf("-- JD[f]=%13.5f %4.4d-%2.2d-%2.2d %2.2d:%2.2d:%02.0f [%2.2d:%2.2d",jdf,date.y,date.mo,date.d,date.h,date.mn,date.s,(int)lst(jdf,longit),((int)((lst(jdf,longit)-(int)lst(jdf,longit))*60.+0.5)));
	if ( RA >=0. ) {
	  printf(",");
	  myput_coords(hainm12top12(jdf,longit,RA),0);
	  printf(",%6.2f",secant_z(altit(Dec,hainm12top12(jdf,longit,RA),lat,&az)));
          if ( bPrintOpts%2 > 0 ) printf(",%6.1f",parang((lst(jdf,longit)-RA),Dec,lat));
          if ( bPrintOpts%4 > 1 ) printf(",%6.1f,%6.1f",altit(Dec,(lst(jdf,longit)-RA),lat,&az),az);
	}
	printf("] ");
	if ( phf >= 0. ) printf("%6.4f ",phf);
	if ( moon[0] > 0. ) printf(" %4.2f  ",moon[11]);
	if (( moon[0] > 0. ) && ( RA >=0. )) printf("%5.1f  ",(180.*subtend(RA,Dec,moon[9],moon[10])/PI));
	// The window duration, exposure and overhead
	hr=(jdf-jdi)*24.;
	min=((jdf-jdi)*24.-hr)*60.;
	printf("  %2.2d:%2.2d  ",hr,min);
	if ( (ed+preoh+postoh) > 0 ) {
	  emin=(ed+preoh+postoh)/60;
	  esec=(ed+preoh+postoh)%60;
	  printf("%4d %4d  %2.2d:%2.2d ",ed,preoh,emin,esec);
	}
	//printf("  Observable\n",hr,min,ed,preoh,emin,esec);
	if ( label != 0 ) printf(" %s",label);
      }
      printf("\n");
    }
  } else {
    /*
"    |             |         |      |||--------------------Start of Condition-------------------------------------||||---------------------End of Condition--------------------------------------|        |    |Pre-|     |       "
"Cond|  Civil-Date | JD@Midn | Moon |||         JD        | UT-Date  |UT-Time | LST  |  HA  |   X  | phase | Moon ||||         JD        | UT-Date  |UT-Time |  LST |  HA  |   X  | phase | Moon |Duration| Exp| OH | Exec|  Label"
"    |             |         |      |||                   |yyyy-mm-dd|hh:mm:ss|      | hh:mm|      |       |      ||||                   |yyyy-mm-dd|hh:mm:ss|      | hh:mm|      |       |      |  hh:mm | (s)| (s)|mm:ss|       "
"EAST 2000-12-01/02 2451880.7  0.00  : JD[s]=2451880.68793 2000-12-02 04:30:37 [04:33,-00:27,  1.24] 0.8472  0.00  -- JD[f]=2451880.84895 2000-12-02 08:22:29 [08:26, 03:26,  1.49] 0.8774  0.00    03:51   100  120 05:50   osmc"
"EAST 2000-11-07/08 2451856.7  0.00  : JD[s]=2451856.58865 2000-11-08 02:07:39 [00:35,-00:27,  1.38] 0.4764  0.00  -- JD[f]=2451856.85667 2000-11-08 08:33:36 [07:02, 06:00,  2.15] 0.5061  0.00    06:25  osmc09-117107"
"EAST 2000-11-07/08 2451856.7  0.00  : JD[s]=2451856.51262 2000-11-08 00:18:10 [22:45, 02:42,  1.58]  0.00  -- JD[f]=2451856.56738 2000-11-08 01:37:02 [00:04, 04:02,  2.49]  0.00   01:18"
"SSSR 2000-11-07/08 2451856.7  0.00  : JD[s]=2451856.46712 2000-11-07 23:12:39 [21:40]  0.00  -- JD[f]=2451856.90308 2000-11-08 09:40:26 [08:09]  0.00  10:27"
    */
    printf("    |             |         ");
    if ( moon[0] > 0. ) {
      printf("| ");
      if ( moon[0] > 1. ) printf("      ");
      printf("     ");
      if ( moon[0] > 1. ) printf("      ");
    }
    printf("|||");
    // Start...
    if ( RA >=0. ) {
      printf("-------");
      if ( bPrintOpts%2 > 0 ) printf("----");
      if ( bPrintOpts%4 > 1 ) printf("-------");
    }
    if ( phi >= 0. ) printf("----");
    if ( moon[0] > 0. ) printf("---");
    if (( moon[0] > 0. ) && ( RA >=0. )) printf("---");
    printf("--------------------Start of Condition--------");
    if ( RA >=0. ) {
      printf("-------");
      if ( bPrintOpts%2 > 0 ) printf("---");
      if ( bPrintOpts%4 > 1 ) printf("-------");
    }
    printf("-");
    if ( phi >= 0. ) printf("---");
    if ( moon[0] > 0. ) printf("----");
    if (( moon[0] > 0. ) && ( RA >=0. )) printf("----");
    // End...
    printf("-||");
    if ( phi >= 0. ) printf("----");
    if ( RA >=0. ) {
      printf("-------");
      if ( bPrintOpts%2 > 0 ) printf("----");
      if ( bPrintOpts%4 > 1 ) printf("-------");
    }
    printf("----------------------End of Condition---------");
    if ( moon[0] > 0. ) printf("----");
    if (( moon[0] > 0. ) && ( RA >=0. )) printf("---");
    if ( RA >=0. ) {
      printf("-------");
      if ( bPrintOpts%2 > 0 ) printf("---");
      if ( bPrintOpts%4 > 1 ) printf("-------");
    }
    printf("-");
    if ( phi >= 0. ) printf("---");
    if ( moon[0] > 0. ) printf("---");
    if (( moon[0] > 0. ) && ( RA >=0. )) printf("----");
    // Duration...
    printf("|        |");
    if ( ed+preoh+postoh > 0  ) printf("    |Pre-|      |");
    if ( label != 0 ) printf("     ");
    printf("\n");

    printf("Cond|  Civil-Date | JD@Midn ");
    if ( moon[0] > 0. ) {
      printf("| ");
      if ( moon[0] > 1. ) printf("      ");
      printf("Moon ");
      if ( moon[0] > 1. ) printf("      ");
    }
    printf("|||         JD        | UT-Date  |UT-Time |  LST ");
    if ( RA >=0. ) {
      printf("|  HA  |   X  ");
      if ( bPrintOpts%2 > 0 ) printf("|  PA  ");
      if ( bPrintOpts%4 > 1 ) printf("|  Alt |  Az  ");
    }
    printf(" ");
    if ( phi >= 0. ) printf("|phase ");
    if ( moon[0] > 0. ) {
      if ( RA >= 0. ) {
        printf("|    Moon     ");
      } else {
        printf("| Moon ");
      }
    }
    printf(" ||          JD        | UT-Date  |UT-Time |  LST ");
    if ( RA >=0. ) {
      printf("|  HA  |   X  ");
      if ( bPrintOpts%2 > 0 ) printf("|  PA  ");
      if ( bPrintOpts%4 > 1 ) printf("|  Alt |  Az  ");
    }
    printf(" ");
    if ( phi >= 0. ) printf("|phase ");
    if ( moon[0] > 0. ) {
      if ( RA >= 0. ) {
        printf("|    Moon     ");
      } else {
        printf("| Moon ");
      }
    }
    printf("|Duration|");
    if ( ed+preoh+postoh > 0  ) printf(" Exp| OH | Exec |");
    if ( label != 0 ) printf(" Label");
    printf("\n");

    printf("    |             |         ");
    if ( moon[0] > 0. ) {
      printf("| ");
      printf("frac ");
      if ( moon[0] > 1. ) printf(" rise ");
      if ( moon[0] > 1. ) printf(" set  ");
    }
    printf("|||                   |yyyy-mm-dd|hh:mm:ss|      ");
    if ( RA >=0. ) {
      printf("| hh:mm|      ");
      if ( bPrintOpts%2 > 0 ) printf("| degr ");
      if ( bPrintOpts%4 > 1 ) printf("| degr | degr ");
    }
    printf(" ");
    if ( phi >= 0. ) printf("|      ");
    if ( moon[0] > 0. ) printf("| frac ");
    if (( moon[0] > 0. ) && ( RA >=0. )) printf("  dist ");
    printf(" ||                    |yyyy-mm-dd|hh:mm:ss|      ");
    if ( RA >=0. ) {
      printf("| hh:mm|      ");
      if ( bPrintOpts%2 > 0 ) printf("| degr ");
      if ( bPrintOpts%4 > 1 ) printf("| degr | degr ");
    }
    printf(" ");
    if ( phi >= 0. ) printf("|      ");
    if ( moon[0] > 0. ) printf("| frac ");
    if (( moon[0] > 0. ) && ( RA >=0. )) printf("  dist ");
    printf("|  hh:mm |");
    if ( ed+preoh+postoh > 0  ) printf(" (s)| (s)| mm:ss|");
    if ( label != 0 ) printf("     ");
    printf("\n");

    printf("----|-------------|---------");
    if ( moon[0] > 0. ) {
      printf("|-");
      if ( moon[0] > 1. ) printf("------");
      printf("-----");
      if ( moon[0] > 1. ) printf("------");
    }
    printf("|||-------------------|----------|--------|------");
    if ( RA >=0. ) {
      printf("|------|------");
      if ( bPrintOpts%2 > 0 ) printf("|------");
      if ( bPrintOpts%4 > 1 ) printf("|------|------");
    }
    printf("-");
    if ( phi >= 0. ) printf("|------");
    if ( moon[0] > 0. ) printf("|------");
    if (( moon[0] > 0. ) && ( RA >=0. )) printf("-------");
    printf("-||--------------------|----------|--------|------");
    if ( RA >=0. ) {
      printf("|------|------");
      if ( bPrintOpts%2 > 0 ) printf("|------");
      if ( bPrintOpts%4 > 1 ) printf("|------|------");
    }
    printf("-");
    if ( phi >= 0. ) printf("|------");
    if ( moon[0] > 0. ) printf("|------");
    if (( moon[0] > 0. ) && ( RA >=0. )) printf("-------");
    printf("|--------|");
    if ( ed+preoh+postoh > 0  ) printf("----|----|------|");
    if ( label != 0 ) printf("------>");
    printf("\n");
  }
  return;
}

void setTelescope( const char *telescope, short bptrec, short btrec, short *bCS, short *insz, double *longit, double *lat, double *elev, const char *EWtext[])
{
      while ( btrec == 0 ) {
        if ( ! strcmp(telescope,"LSO")) {
          btrec=1;
          *bCS=0;
          // From ...
          *longit = ((70.+44./60.+04.540/3600.)/180.*12.);
          *lat = (-29.-15./60.-15.430/3600.);
          *elev = 2335.00;
          *insz=1;
          EWtext[0]="AML ";
        }

        if ( ! strcmp(telescope,"2p2")) {
          btrec=1;
          *bCS=0;
          // From ...
          *longit = ((70.+44./60.+04.540/3600.)/180.*12.);
          *lat = (-29.-15./60.-15.430/3600.);
          *elev = 2335.00;
          *insz=1;
          EWtext[0]="AML ";
        }

        if ( ! strcmp(telescope,"NTT")) {
          btrec=1;
          *bCS=0;
          // From ...
          *longit = ((70.+44./60.+04.540/3600.)/180.*12.);
          *lat = (-29.-15./60.-15.430/3600.);
          *elev = 2335.00;
          *insz=1;
          EWtext[0]="AML ";
        }

        if ( ! strcmp(telescope,"3p6")) {
          btrec=1;
          *bCS=0;
          // From ...
          *longit = ((70.+44./60.+04.540/3600.)/180.*12.);
          *lat = (-29.-15./60.-15.430/3600.);
          *elev = 2335.00;
          *insz=1;
          EWtext[0]="AML ";
        }

        if ( ! strcmp(telescope,"DK154")) {
          btrec=1;
          *bCS=1;
          *longit = ((70.+44./60.+07.662/3600.)/180.*12.);
          *lat = (-29.-15./60.-14.235/3600.);
          *elev = 2340.00;
          *insz=2;
        }

        if ( ! strcmp(telescope,"ESO152")) {
          btrec=1;
          // http://www.ls.eso.org/lasilla/Telescopes/2p2T/E1p5M/
          // The telescope is at a geographical location of 70°44'12"865W 29°15'7"422S and an altitude of 2335 metres.
          *longit = ((70.+44./60.+12.865/3600.)/180.*12.);
          *lat = (-29.-15./60.-07.422/3600.);
          *elev = 2335.00;
        }

        if (( ! strcmp(telescope,"LasCampanas") ) || ( ! strcmp(telescope,"OGLE") )) {
          btrec=1;
          *bCS=0;
          // From XEphem...
          *longit = ((70.+42./60.+00.000/3600.)/180.*12.);
          *lat = (-29.-00./60.-00.000/3600.);
          *elev = 2300.00;
          *insz=1;
          EWtext[0]="AML ";
        }

        if ( ! strcmp(telescope,"MJUO") ) {
          btrec=1;
          *bCS=0;
          // Latitude 43 deg 59.2 min S; *Longitude 170 deg 27.9 min Altitude 1029 metres
          *longit = -1.*((170.+27.9/60.+00.000/3600.)/180.*12.);
          *lat = -1.*(43.+59.2/60.+00.000/3600.);
          *elev = 1029.00;
          *insz=1;
          EWtext[0]="AML ";
        }

        if ( (! strcmp(telescope,"VLT")) || (! strcmp(telescope,"UT1")) ) {
          /*
            VLT, ESO, Paranal
            From http://www.hq.eso.org/paranal/site/paranal.html:
            =====================================================
            The Paranal Observatory is located on the top of Cerro Paranal in
            the Atacama Desert in the northern part of Chile and what is
            believed to be the driest area on Earth. Cerro Paranal is a 2,635-m
            high mountain, about 120 km south of the town of Antofagasta and 12
            km inland from the Pacific Coast. The geographical coordinates are
            24o 40' S, 70o 25' W.


            But at http://www.hq.eso.org/paranal/site/paranal.html#GeoInfo
            ==============================================================
            Altitude above sea level        2635.43 m
            Distance from coast             12 km
            Highest neighboring peak        Cerro la Chira, 2569 m, 11 km to the NNE
            Road distances                  130 km to Antofagasta
                                            Appr. 1200 km to Santiago
                                            Appr. 600 km to La Silla

            Telescope number                Latitude                       *Longitude
            UT1                             -24deg 37 min 33.117 sec       70 deg 24 min 11.642 sec
            UT2                             -24deg 37 min 31.465 sec       70 deg 24 min 10.855 sec
            UT3                             -24deg 37 min 30.300 sec       70 deg 24 min 9.896 sec
            UT4                             -24deg 37 min 31.000 sec       70 deg 24 min 8.000 sec

            // UT2::
          */
          btrec=1;
          *bCS=0;
          *longit = ((70.+24./60.+11.642/3600.)/180.*12.);
          *lat = (-24.-37./60.-33.117/3600.);
          *elev = 2635.43;
          *insz=1;
          EWtext[0]="AML ";
        }

        if ( ! strcmp(telescope,"UT2") ) {
          btrec=1;
          *bCS=0;
          *longit = ((70.+24./60.+10.855/3600.)/180.*12.);
          *lat = (-24.-37./60.-31.465/3600.);
          *elev = 2635.43;
          *insz=1;
          EWtext[0]="AML ";
        }

        if ( ! strcmp(telescope,"UT3") ) {
          btrec=1;
          *bCS=0;
          *longit = ((70.+24./60.+9.896/3600.)/180.*12.);
          *lat = (-24.-37./60.-30.300/3600.);
          *elev = 2635.43;
          *insz=1;
          EWtext[0]="AML ";
        }

        if ( ! strcmp(telescope,"UT4") ) {
          btrec=1;
          *bCS=0;
          *longit = ((70.+24./60.+8.000/3600.)/180.*12.);
          *lat = (-24.-37./60.-31.000/3600.);
          *elev = 2635.43;
          *insz=1;
          EWtext[0]="AML ";
        }

        if ( btrec == 0 ) {
          printf("Warning!!! Telescope %s NOT recognised. ",telescope);
          telescope="LSO";
          printf("Using DEFAULT telescope %s instead.\n",telescope);
        } else {
          if ( bptrec == 1 ) printf("Using telescope %s instead.\n",telescope);
        }
      }
      return;
}
