/* "@(#) $Id: difrefrac.cc 3021 2017-02-01 08:23:08Z jpritcha $" */
/*
  This is difrefrac.c

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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <getopt.h>
#include <libgen.h>
#include <cpgplot.h>

#include "libdifrefrac.h"

void usage( char ** argv )
{
  /*
    How this programme (should) work...
    1) You give it the <T0> <P> <RA> <Dec> <phase-i> <phase-f> <jd-i> <jd-f> of the MCEB
    2) For each day between jd-i & jd-f it...
      (i)  Calculates sunrise and sunset...
     (ii)  Calculates the time of day at which the object is within the safty zones (E & W of pier)
    (iii)  Calculates the interval within that determined in (ii) during which the object is within the specified phase range.
  */

  printf ("Usage : %s [-h|--help|-v]\n\
             [-PGXWZ]\n\
             [-p <UpDate-Period> [integer sec] (3600)>]\n\
             [-a <HA[i],HA[f],HA[s] [decimal hrs] (-6.0,6.0,1/3600.)>]\n\
             [-d <Dec[i],Dec[f],Dec[s] [decimal degrees] (-80.0,+20.0,10.0)>]\n\
             [-w <WL[i],WL[f],WL[s] [nm] (350.0,350.0,200.0>]\n\
             [-x <x[i],x[f] [Zenith Distance range] (0.0,70.0)>]\n\
             [-y <y[i],y[f] [arcsec] (0.0,5.0)>]\n\
             [-z <minZ>]\n\
             [-r <reference-WL [nm] (550.454 i.e. Johnson V)>]\n\
             [-A <T,P [C,mb] (10.0,775.0)>]\n\
             [-T <telescope (2p2)>]\n\
             [-D] (differential refraction rather than UDP separation)\n\
             [-F]\n\
\n\
Note:\n\
 *)  -G, Make output PNG files, default file names.\n\
 *)  -P, Make output PS files, default file names.\n\
 *)  -X, No X display.\n\
 *)  -Z, Plot with respect to Zenith Distance (deafult is wrt Hour Angle).\n\
 *)  -Y, Plot with respect to Airmass (deafult is wrt Hour Angle).\n\
 *)  -W, Wait for <enter> after making plots\n\
 *)  <RA,Dec[,Epoch]> should be specified in hh:mm:ss.s,dd:mm:ss.s[,yyyy.y]\n\
format. If no coordinates are specified Sunset and Sunrise [SSSR] and Begining\n\
and End of Night [BEoN] are calcualted. If co-ordinates, but no Epoch are\n\
specified 2000.0 is assumed.\n\
 *)  Valid <telescope> are LSO, 2p2, NTT, 3p6, DK154, ESO152, VLT, UT1, UT2,\n\
UT3, UT4, MJUO, and LasCampanas==OGLE. The default is LSO.\n\
\n\
",basename(argv[0]));

}

void
versioninfo( char **argv )
{
  // $Revision: 3021 $ and $Date: 2017-02-01 09:23:08 +0100 (Wed, 01 Feb 2017) $ are CVS keywords
  /*
  char *rev={"$Revision: 3021 $"};
  char *date={"$Date: 2017-02-01 09:23:08 +0100 (Wed, 01 Feb 2017) $"};
  char *rstr,*dstr;
  rstr=strtok(rev,"$Revision: ");
  printf("%s\n",rstr);
  rstr=strtok(NULL," ");
  printf("%s\n",rstr);
  dstr=strtok(date,"$Date: "); dstr=strtok(NULL," ");
  printf("%s - Version %s, %s
  */
  printf("%s - Version , $Revision: 3021 $, $Date: 2017-02-01 09:23:08 +0100 (Wed, 01 Feb 2017) $\n\
\n\
By J.D.Pritchard <j.pritchard@eso.org>\n\
Based on the skycalc.c and skycalendar.c codes by John Thorstensen\n\
",basename(argv[0]));
}

double nmo ( double WL, double P, double T )
{
  /* See:
     http://www.ls.eso.org/lasilla/sciops/2p2/E2p2M/FEROS/Projects/ADC/references/refraction/index.html
  */
  double nsmo, wlmusq;
  double Ps, Ts, h0;

  Ps=1013.25E+2; /* Pa */
  Ts=288.15;     /* K  */
  h0=7.0;        /* km */
  wlmusq=WL*WL/1000./1000.;

  nsmo=(64.328+(29498.1/(146.0-(1.0/wlmusq)))+(255.4/(41.0-(1.0/wlmusq))))/1.0E+6;
  /*
  printf("WL=%f : P=%f : T=%f : nsmo=%f :: n=%f\n",WL,P,T,nsmo,(1.0+nsmo*P/T*Ts/Ps));
  */
  return(nsmo*P/T*Ts/Ps);
}

double n ( double WL, double P, double T )
{
  /* See:
     http://www.ls.eso.org/lasilla/sciops/2p2/E2p2M/FEROS/Projects/ADC/references/refraction/index.html
  */
  return(1.0+nmo(WL,P,T));
}

int main ( int argc, char **argv )
{
  double Pa, Ta, Zs, Zm, Ze, difrefrac;
  double refWL, WL, WLi, WLf, WLs;
  double rdec, rdeci, rdecf, rdecs;
  double rha, rhai, rhaf, rhas;
  double minZ;
  double PrismAngle;

  float rse, rsm, rme;
  float drf, ddrf, *vx, *vy;
  int nhapts;
  float xs, xf, ys, yf, yt;

  bool bXdisp=1, bPwrtX=0, bPwrtZ=0, bPwrtD=0, bMkPS=0, bAutoPSFileNames=0, bMkpng=0, bWaitAP=0;
  short psDev, XDev, pngDev;
  char ofile[200], title[300], ystr[100];

  int udp=0,preoh=0,postoh=0,y,m,d,j,icn,i,iew,hr,min,sec,emin,esec,icloser;
  double jdmid,mjds,fd;
  double T0=0.,T0shift=0.,P=-1.0,phi=0.,phishift=0.,phf=1.,phfshift=0.,jdi,jdf,jdph[2][2],ha[2][3],jdha[2],sssr[2][2],ssssr[2];

  double maxX=100.00,alt,az,minAlt=0.,haAtMinAlt[2],tha;
  double altBEoN=-13.0;
  double RA,Dec,tRA,tDec,Epoch,Equinox;
  const char *str;
  short  dow;
  struct date_time date,di,df;
  extern char *_tzname[2];

  double jdUTC;

  double moon[11],minMoonDist,minMoonDistTonight,maxMoonDistTonight,moond[3];
  short  bdmoon=0;

  double jd,ph[2][2],pht[2][3],safezone[2][2];
  short bfpr=0,bCalcNow=0,bEphem=0,bError=0,bPrintOpts=0,bPhase=0,bQuiet=0;
  short btrec=0;
  short bptrec=0;

  const char *telescope={"DK154"};
  const char *label={0};
  const char *EWtext[]={"EAST","WEST",0};

  double geora, geodec, geodist;  /* geocent for moon, not used here.*/
  double rasun, decsun;
  double stmid, ramoon, decmoon, distmoon;

  /*
    Danish-1.54, ESO, La Silla

    An email from Heath Jones, 2p2 Team member...
    =============================================
    Date: Fri, 13 Oct 2000 15:59:37 -0400
    From: "2.2 telescope team" <2p2team@eso.org>
    Organization: European Southern Observatory
    To: John Pritchard <j.pritchard@astro.ku.dk>
    CC: Patrick FRANCOIS <fpatrick@eso.org>, James Brewer <jbrewer@astro.ubc.ca>
    Subject: Re D1.54 Safety Zone map

    <<<SNIP>>>

    The exact coordinates of the D154 are 
        70 deg 44' 07"662 W
        29 deg 15' 14"235 S
    The telescope is at an altitude of 2340 m.

    <<<SNIP>>>

  */

// Danish-1.54m::
  double longit = ((70.+44./60.+07.662/3600.)/180.*12.), *plongit;
  double lat = (-29.-15./60.-14.235/3600.), *plat;
  double elev = 2340.00, *pelev;
  short bCS, *pbCS;
  short insz=2, *pinsz;

  /*
    From the gettimeofday man ctime:
SYNOPSIS
       #include <time.h>

       char *asctime(const struct tm *timeptr);

       char *ctime(const time_t *timep);

       struct tm *gmtime(const time_t *timep);

       struct tm *localtime(const time_t *timep);

       time_t mktime(struct tm *timeptr);

       extern char *_tzname[2];
       long int timezone;
       extern int daylight;

    From man time
SYNOPSIS
       #include <time.h>

       time_t time(time_t *t);

DESCRIPTION
       time returns the time since the Epoch (00:00:00 UTC, January 1, 1970), measured in seconds.

       If t is non-NULL, the return value is also stored in the memory pointed to by t.

RETURN VALUE
       On success, the value of time in seconds since the Epoch is returned.  On error, ((time_t)-1) is returned, and errno is set appropriately.
  */
  long int timezone;
  extern int daylight;
  struct tm gmt,lt;
  time_t ct,tzoff;

  int c;
  int digit_optind = 0;
  int this_option_optind = optind ? optind : 1;
  int option_index = 0;
  struct option longopts[] =
  {
    /* { name  has_arg  *flag  val } */
    {"BEoNAlt",       1, 0, 's'}, /* Altitude for Begining/End of night calculation */
    {"airmass",       1, 0, 'X'}, /* Maximum Airmass */
    {"altitude",      1, 0, 'A'}, /* Minimum Altitude */
    {"PwrtD",         1, 0, 'D'}, /* Print Differential Refraction not UDP seperation */
    {"coordinates",   1, 0, 'c'}, /* Coordinates : RA,Dec */
    {"disttomoon",    1, 0, 'd'}, /* Distance to the moon */
    {"ephemeris",     1, 0, 'E'}, /* Ephemeris */
    {"expduration",   1, 0, 'e'}, /* Exposure Duration */
    {"help",          0, 0, 'h'}, /* help */
    {"jd",            1, 0, 'j'}, /* Julian Date interval */
    {"label",         1, 0, 'l'}, /* Label */
    {"moon",          0, 0, 'm'}, /* Display moon parameters */
    {"overhead",      1, 0, 'o'}, /* Overhead */
    {"phaseinterval", 1, 0, 'p'}, /* Phase interval */
    {"PParAng",       1, 0, 'P'}, /* Print Paralactic Angle */
    {"PwrtX",         1, 0, 'Y'}, /* Print with respect to Airmass */
    {"PwrtZ",         1, 0, 'Z'}, /* Print with respect to ZD */
    {"now",           0, 0, 'N'}, /* Calculate for "Now" */
    {"UTC",           1, 0, 'U'}, /* Calculate for "UTC */
    { 0, 0, 0, 0 }
  };

  ct=time(NULL)+12*3600;
  lt=*localtime( &ct );
  lt.tm_hour=0;
  lt.tm_min=0;
  lt.tm_sec=0;
  ct=mktime( &lt );
  gmt=*gmtime( &ct );
  tzoff=(gmt.tm_hour)*3600+(gmt.tm_min)*60+(gmt.tm_sec);
  //printf("Local Date at nearest MidNight is %4d-%2.2d-%2.2d\n",(1900+lt.tm_year),(1+lt.tm_mon),lt.tm_mday);
  //printf("GMT Date:Time at nearest MidNight is %4d-%2.2d-%2.2d %2.2d:%2.2d:%2.2d\n",(1900+gmt.tm_year),(1+gmt.tm_mon),gmt.tm_mday,gmt.tm_hour,gmt.tm_min,gmt.tm_sec);
  di.y=(1900+gmt.tm_year);di.mo=(gmt.tm_mon+1);di.d=gmt.tm_mday;di.h=gmt.tm_hour,di.mn=gmt.tm_min;di.s=gmt.tm_sec;
  jdi=jdf=date_to_jd( di );
  //printf("Corresponding JD is : %9.1f\n",jdi);
  RA=-1.0;Dec=0.;
  moon[0]=0.;
  Equinox=2000.0;
  minMoonDist=-200.;
  maxMoonDistTonight=-199.;
  jdUTC=-1.;

  /*
    Check for how the program was call...
    Recognised options...
      dk154sc
      eso152sc
      2p2sc
      vltsc
  */

  telescope="LSO";
  if ( ! strcmp(basename(argv[0]),"dk154sc" ) ) telescope="DK154";
  if ( ! strcmp(basename(argv[0]),"eso152sc") ) telescope="ESO152";
  if ( ! strcmp(basename(argv[0]),"2p2sc"   ) ) telescope="2p2";
  if ( ! strcmp(basename(argv[0]),"nttsc"   ) ) telescope="NTT";
  if ( ! strcmp(basename(argv[0]),"3p6sc"   ) ) telescope="3p6";
  setTelescope( telescope, bptrec, btrec, &bCS, &insz, &longit, &lat, &elev, EWtext );

  /* Set default values */

  refWL=550.454;  /* nm */
  Pa=775.0;       /* mb */
  Ta=10.0;        /* C  */
  WL=350.00;      /* nm */
  udp=20;          /* sec */
  WLi=350.0;WLf=350.0;WLs=200.0;
  rdeci=-80.0;rdecf=+20.0;rdecs=10.0;
  rhai=-6.0;rhaf=(6.0-udp/3600.0);rhas=1.0/3600.0;
  minZ=70.0;
  xs=0.0; xf=+70.0; ys=0.0; yf=0.5; yt=0.2;
  PrismAngle=1.5/60.0/60.0;

  while (1) {
    c = getopt_long (argc, argv, "A:DGNPT:X:U:WXZa:c:d:e:E:hj:l:mo:p:qr:s:t:vw:x:y:z:",
                     longopts, &option_index);
    if (c == -1)
      break;

    switch (c) {
    case 0:
      printf ("option %s", longopts[option_index].name);
      if (optarg)
        printf (" with arg %s", optarg);
      printf ("\n");
      break;


    case 'A':
      sscanf(optarg,"%lf,%lf",&Ta,&Pa);
      break;

    case 'D':
      bPwrtD=1;
      break;

    case 'E':
      sscanf(optarg,"%lf,%lf",&T0,&P);
      //printf ("option E with value `%s' :: Ephemeris T0=%20.6f, P=%20.5f\n", optarg,T0,P);
      bEphem=1;
      break;

    case 'G':
      bMkpng=1;
      break;

    case 'N':
      bCalcNow=1;
      break;

    case 'P':
      bMkPS=1;
      bAutoPSFileNames=1;
      break;

    case 'T':
      telescope=optarg;
      bptrec=1;
      setTelescope( telescope, bptrec, btrec, &bCS, &insz, &longit, &lat, &elev, EWtext );
      break;

    case 'U':
      str=strtok(optarg,":");
      if ( str != NULL ) {
        di.h=atoi( str );di.mn=atoi( strtok(NULL,":"));di.s=atof( strtok(NULL,","));
        str=strtok(NULL,"-");
        if ( str != NULL ) {
          //if ( index( str,'-' ) != NULL ) {
          di.y=atoi( str );di.mo=atoi( strtok(NULL,"-"));di.d=atoi( strtok(NULL,","));
          //} else if ( index( str,'/' ) != NULL ) {
          //   di.y=atoi( strtok(NULL,"/") );di.mo=atoi( strtok(NULL,"/"));di.d=atoi( strtok(NULL,","));
          //}
        } else {
          ct=time(NULL);
          lt=*localtime( &ct );
          ct=mktime( &lt );
          gmt=*gmtime( &ct );
          di.y=(1900+gmt.tm_year);di.mo=(gmt.tm_mon+1);di.d=gmt.tm_mday;
          //di.h=gmt.tm_hour,di.mn=gmt.tm_min;di.s=gmt.tm_sec;
        }
        jdUTC = date_to_jd( di );
      } else {
        jdUTC = -1.;
      }
      break;

    case 'X':
      bXdisp=0;
      //sscanf(optarg,"%lf",&maxX);
      //printf ("option X with value `%s'  :: Maximum Airmass'\n", optarg);
      break;

    case 'W':
      bWaitAP=1;
      break;

    case 'Y':
      printf('Plotting with respect to Airmass');
      bPwrtX=1;
      break;

    case 'Z':
      bPwrtZ=1;
      break;

    case 'a':
      sscanf(optarg,"%lf,%lf,%lf",&rhai,&rhaf,&rhas);
      break;

    case 'c':
      str=strtok(optarg,",");
      if ( str != NULL ) {
        RA=get_coorddes( str );
        str=strtok(NULL,",");
        if ( str != NULL ) {
          Dec=get_coorddes( str );
          str=strtok(NULL,",");
          if ( str != NULL ) Equinox=atof( str );
        } else {
          RA=-1.;
        }
      }
      break;

    case 'd':
      sscanf(optarg,"%lf,%lf,%lf",&rdeci,&rdecf,&rdecs);
      break;

    case 'h':
      usage( argv );
      exit(0);
      break;

    case 'j':
      if ( index( optarg,'-' ) != NULL ) {
        di.y=atoi( strtok(optarg,"-") );di.mo=atoi( strtok(NULL,"-"));di.d=atoi( strtok(NULL,","));
        df.y=atoi( strtok(NULL,"-"));df.mo=atoi( strtok(NULL,"-"));df.d=atoi( strtok(NULL,"-"));
        di.h=(int) (12.+longit);di.mn=0;di.s=0;df.h=(int) (12.+longit);df.mn=0;df.s=0;
        jdi=date_to_jd( di );jdf=date_to_jd( df );
      } else if ( index( optarg,'/' ) != NULL ) {
        di.y=atoi( strtok(optarg,"/") );di.mo=atoi( strtok(NULL,"/"));di.d=atoi( strtok(NULL,","));
        df.y=atoi( strtok(NULL,"/"));df.mo=atoi( strtok(NULL,"/"));df.d=atoi( strtok(NULL,"/"));
        di.h=(int) (12.+longit);di.mn=0;di.s=0;df.h=(int) (12.+longit);df.mn=0;df.s=0;
        jdi=date_to_jd( di );jdf=date_to_jd( df );
      } else {
        sscanf(optarg,"%lf,%lf",&jdi,&jdf);
      }
      break;

    case 'l':
      label=optarg;
      break;

    case 'm':
      bdmoon=1;
      moon[0]++;
      break;

    case 'o':
      sprintf(ofile,"%s/CPS",optarg);
      bMkPS=1;
      bAutoPSFileNames=0;
      break;

    case 'p':
      sscanf(optarg,"%d",&udp);
      break;

    case 'q':
      bQuiet=1;
      break;

    case 'r':
      sscanf(optarg,"%lf",&refWL);
      break;

    case 's':
      sscanf(optarg,"%lf",&altBEoN);
      break;

    case 't':
      sscanf(optarg,"%f",&yt);
      break;

    case 'v':
      versioninfo( argv );
      exit(0);
      break;

    case 'w':
      sscanf(optarg,"%lf,%lf,%lf",&WLi,&WLf,&WLs);
      break;

    case 'x':
      sscanf(optarg,"%f,%f",&xs,&xf);
      break;

    case 'y':
      sscanf(optarg,"%f,%f",&ys,&yf);
      break;

    case 'z':
      sscanf(optarg,"%lf",&minZ);
      break;

    case '?':
      bError=1;
      break;

    default:
      printf ("?? getopt returned character code 0%o ??\n", c);
    }
  }

  /*
  if (( (argc-optind) != 2 ) || ( bError )) {
    usage( argv );
    exit (1);
  }
  */

  Ta+=273.15;  /* C --> K  */
  Pa=Pa*100.0; /* mb --> Pa */

  if (( bPhase ) && ( ! bEphem )) {
    printf("ERROR!!! You MUST specify -T <T0,P> if you specify a phase range.");
    exit (1);
  }

  if ( phi > phf ) {
    phishift=-1.*phf;
    phfshift=1.-1.*phf;
    T0shift=-1.*P*phishift;
  }

  WL=WLi;

  nhapts=(int)((rhaf-rhai)/rhas+1);
  vx=new float[nhapts];
  vy=new float[nhapts];

  while ( WL < (WLf+WLs) ) {
    j=0;
    printf("drf-%04d-%05.1f-%05.1f-%05.1f-%05.1f\n",udp,WL,refWL,Ta,Pa/100.);
    if ( bMkpng ) {
      sprintf(ofile,"drf-%04d-%05.1f-%05.1f-%05.1f-%05.1f.png/PNG",udp,WL,refWL,Ta,Pa/100.);
      pngDev=cpgopen(ofile);
      if(pngDev < 1) return EXIT_FAILURE;
      cpgsci(1);
    }
    if ( bMkPS ) {
      if ( bAutoPSFileNames ) {
        sprintf(ofile,"drf-%04d-%05.1f-%05.1f-%05.1f-%05.1f.ps/CPS",udp,WL,refWL,Ta,Pa/100.);
      }
      psDev=cpgopen(ofile);
      if(psDev < 1) return EXIT_FAILURE;
      cpgsci(1);
    }
    /* if(cpgbeg(0, "/XWINDOW", 1, 1) != 1) return EXIT_FAILURE; */
    if ( bXdisp ) {
      XDev=cpgopen("/XWINDOW");
      if(XDev < 1) return EXIT_FAILURE;
      cpgsci(1);
      cpgask(bWaitAP);
    }
    /*
     * Call PGENV to specify the range of the axes and to draw a box, and
     * PGLAB to label it. The x-axis runs from 0 to 10, and y from 0 to 20.
     */
    if ( bMkpng ) {
      cpgslct(pngDev);
      if ( bPwrtX ) {
        cpgenv(1./cos(xs*DEG_IN_RADIAN), 1./cos(xf*DEG_IN_RADIAN), ys, yf, 0, 1);
      } else {
        cpgenv(xs, xf, ys, yf, 0, 1);
      }
    }
    if ( bMkPS ) {
      cpgslct(psDev);
      cpgenv(xs, xf, ys, yf, 0, 1);
    }
    if ( bXdisp ) {
      cpgslct(XDev);
      cpgenv(xs, xf, ys, yf, 0, 1);
    }
    str="HA [hrs]";
    printf('bPwrtX');
    printf(bPwrtX);
    if ( bPwrtX ) str="Airmass";
    if ( bPwrtZ ) str="Zenith Distance (degrees)";
    print(str);
    sprintf(ystr,"Separation (r) [arcsec]");
    sprintf(title,"UDP=%d[sec], WL=%5.1f[nm], refWL=%5.1f[nm], Temp=%5.1f[K], Pres=%5.1f[mB]",udp,WL,refWL,Ta,Pa/100.);
    if ( bPwrtD ) {
      sprintf(ystr,"Differential Refraction (d) [arcsec]");
      sprintf(title,"WL=%5.1f[nm], refWL=%5.1f[nm], Temp=%5.1f[K], Pres=%5.1f[mB]",WL,refWL,Ta,Pa/100.);
    }
    if ( bMkpng ) {
      cpgslct(pngDev);
      cpglab(str, ystr, title);
      if (( ys < yt ) && ( yf > yt )) {
        cpgmove(xs,yt);
        cpgdraw(xf,yt);
      }
      if ( bPwrtX && (( xs < 60. ) && ( xf > 060 ))) {
        cpgmove(1./cos(60.*DEG_IN_RADIAN),ys);
        cpgdraw(1./cos(60.*DEG_IN_RADIAN),yf);
      }
      if ( bPwrtZ && (( xs < 60. ) && ( xf > 060 ))) {
        cpgmove(60.,ys);
        cpgdraw(60.,yf);
      }
    }
    if ( bMkPS ) {
      cpgslct(psDev);
      cpglab(str, ystr, title);
      if (( ys < 0.2 ) && ( yf > 0.2 )) {
        cpgmove(xs,0.2);
        cpgdraw(xf,0.2);
      }
      if ( bPwrtX && (( xs < 60. ) && ( xf > 060 ))) {
        cpgmove(1./cos(60.*DEG_IN_RADIAN),ys);
        cpgdraw(1./cos(60.*DEG_IN_RADIAN),yf);
      }
      if ( bPwrtZ && (( xs < 60. ) && ( xf > 060 ))) {
        cpgmove(60.,ys);
        cpgdraw(60.,yf);
      }
    }
    if ( bXdisp ) {
      cpgslct(XDev);
      cpglab(str, ystr, title);
      if (( ys < 0.2 ) && ( yf > 0.2 )) {
        cpgmove(xs,0.2);
        cpgdraw(xf,0.2);
      }
      if ( bPwrtX && (( xs < 60. ) && ( xf > 060 ))) {
        cpgmove(1./cos(60.*DEG_IN_RADIAN),ys);
        cpgdraw(1./cos(60.*DEG_IN_RADIAN),yf);
      }
      if ( bPwrtZ && (( xs < 60. ) && ( xf > 060 ))) {
        cpgmove(60.,ys);
        cpgdraw(60.,yf);
      }
    }
    
    rdec=rdeci;
    while ( rdec < (rdecf+rdecs) ) {
      i=0;
      difrefrac=206264.80625*(nmo(WL,Pa,Ta)-nmo(refWL,Pa,Ta));
      printf("WL[%02d]=%f, n[%5.1f nm]=%f :: n[ref=%5.1f nm]=%f\nC=DifReFrac=%f arcsec=%e rad :: fC=%f\nR[%5.1fnm]=%f :: R[ref=%5.1fnm]=%f\n",
             j,WL,WL,n(WL,Pa,Ta),refWL,n(refWL,Pa,Ta),fabs(difrefrac),fabs(difrefrac)/60/60/DEG_IN_RADIAN,
             fabs(difrefrac)/60/60/DEG_IN_RADIAN*28914,
             WL,(206264.80625*nmo(WL,Pa,Ta)),refWL,(206264.80625*nmo(refWL,Pa,Ta))
             );
      printf("C/A=%f\n",(fabs(difrefrac)/2.0/(PrismAngle*60.0*60.0)));
      rha=rhai;
      while ( rha < (rhaf+rhas) ) {
        /*
From http://www.ls.eso.org/lasilla/sciops/2p2/E2p2M/FEROS/Projects/ADC/index.html
  dR = R(lambda)-R(ref) = 206264.80625 * [ n(lambda) - n(ref) ] * tan (Z)
      */
        /* Zenith Distance at Start of exposure */
        Zs=(90.0-altit(rdec,rha,lat,&az))/DEG_IN_RADIAN;
        /* Zenith Distance at midpoint of exposure */
        Zm=(90.0-altit(rdec,(rha+((double)udp)/2.0/3600.0),lat,&az))/DEG_IN_RADIAN;
        /* Zenith Distance at End of exposure */
        Ze=(90.0-altit(rdec,(rha+((double)udp)/3600.0),lat,&az))/DEG_IN_RADIAN;
        if (( Zs*DEG_IN_RADIAN < minZ ) && ( Ze*DEG_IN_RADIAN < minZ )) {
          rse=fabs(difrefrac*sqrt(tan(Zs)*tan(Zs)+tan(Ze)*tan(Ze)-2*tan(Zs)*tan(Ze)*cos((parang(rha,rdec,lat)-parang((rha+((double)udp)/3600.0),rdec,lat))/DEG_IN_RADIAN)));
          rsm=fabs(difrefrac*sqrt(tan(Zs)*tan(Zs)+tan(Zm)*tan(Zm)-2*tan(Zs)*tan(Zm)*cos((parang(rha,rdec,lat)-parang((rha+((double)udp)/2.0/3600.0),rdec,lat))/DEG_IN_RADIAN)));
          rme=fabs(difrefrac*sqrt(tan(Zm)*tan(Zm)+tan(Ze)*tan(Ze)-2*tan(Zm)*tan(Ze)*cos((parang((rha+((double)udp)/2.0/3600.0),rdec,lat)-parang((rha+((double)udp)/3600.0),rdec,lat))/DEG_IN_RADIAN)));
          drf=difrefrac*tan(Zm);
          ddrf=difrefrac*(tan(Zs)-tan(Ze));
          /*
            printf("%f %d %f %f %f %f %f :: WL=%f[nm] : Exp=%d[sec] : HA=",WL,udp,rha,rdec,Zs*DEG_IN_RADIAN,difrefrac*tan(Zs),difrefrac*(tan(Zs)-tan(Ze)),WL,udp);
            myput_coords(rha,2);
            printf(" : Dec=");
            myput_coords(rdec,2);
            printf(" : Z=%f[deg] : DRF=%f[arcsec] : DeltaDRF=%f[arcsec]\n",Zs*DEG_IN_RADIAN,difrefrac*tan(Zs),difrefrac*(tan(Zs)-tan(Ze)));
          */
          vx[i]=rha+1.0*udp/2.0/3600.0;
          if ( bPwrtX ) vx[i]=1./cos(Zm);
          if ( bPwrtZ ) vx[i]=Zm*DEG_IN_RADIAN;
          vy[i]=rse;
          if ( bPwrtD ) vy[i]=fabs(drf);
          i++;
        }
        rha+=rhas;
      }
      j++;
      if ( bMkpng ) {
        cpgslct(pngDev);
        cpgsci(j);
        cpgline(i,vx,vy);
      }
      if ( bMkPS ) {
        cpgslct(psDev);
        cpgsci(j);
        cpgline(i,vx,vy);
      }
      if ( bXdisp ) {
        cpgslct(XDev);
        cpgsci(j);
        cpgline(i,vx,vy);
      }
      rdec+=rdecs;
    }
    WL+=WLs;
    if ( bMkpng ) {
      cpgslct(pngDev);
      cpgclos();
    }
    if ( bMkPS ) {
      cpgslct(psDev);
      cpgclos();
    }
    if ( bXdisp ) {
      cpgslct(XDev);
      cpgclos();
    }
  }
  /*
   * Mark five points (coordinates in arrays XS and YS), using symbol
   * number 9.
   */
  /*  cpgpt(5, xs, ys, 9); */
  /*
   * Compute the function at 'n=60' points, and use PGLINE to draw it.
   */
  /*
  for(i=0; i<n; i++) {
    xr[i] = 0.1*i;
    yr[i] = xr[i]*xr[i];
  }
  cpgline(n, xr, yr);
  */
  /*
   * Finally, call PGEND to terminate things properly.
   */
  cpgend();
  return EXIT_SUCCESS;

  exit (0);
}
