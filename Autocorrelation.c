#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main (int argc, char **argv)
{
  /* check the format of input */
  if (argc != 4)
    {
      fprintf(stderr, "Usage: %s DATfile skip[sample] frame_length[sample]\n", argv[0]);
      exit(1);
    }
  FILE* fpDAT;
  int nskip;
  int framelen;

  /* check the validity of input */
  if( ( fpDAT = fopen( argv[1], "r" ) ) == NULL )  exit( 1 );
  if( ( nskip    = atoi( argv[2] ) ) < 0 )  exit( 1 );
  if( ( framelen = atoi( argv[3] ) ) < 0 )  exit( 1 );

  fprintf( stderr, "# DATfile = %s\n", argv[1] );
  fprintf( stderr, "# %d samples are skipped.\n", nskip );
  fprintf( stderr, "# 1 frame contains %d sampels.\n", framelen );

  short *sdata = (short*)calloc( framelen, sizeof(short) );
  double *autocorrelation = (double*)calloc(framelen,sizeof(double));
  double *xr = (double*)calloc(framelen,sizeof(double));
  if ( sdata == NULL) exit (1);

  fseek( fpDAT, nskip*sizeof(short), SEEK_SET );
  fread( sdata, sizeof(short), framelen, fpDAT);
  fclose( fpDAT );

  /* windowing */
  for (int i = 0 ; i < framelen; i++)
    {
      xr[i] = (0.54-0.46*cos(2*M_PI*i/(framelen-1))) *sdata[i];
    }


  /* cauculate autocorrelation */
  for (int r = 0; r < framelen; r++){
    for (int t = 0; t < framelen - 1 - r ;t++){
      //autocorrelation[r] += sdata[t]*sdata[t+r];
      autocorrelation[r] += xr[t]*xr[t+r];
    }
  }

  /* normalize */
  double r0=autocorrelation[0];
  for (int i = 0; i < framelen; i++){
    autocorrelation[i]=autocorrelation[i]/r0;
  }

  for (int i = 0; i < framelen; i++){
    printf("%lf\n", autocorrelation[i]);
  }

  return 0;
}
