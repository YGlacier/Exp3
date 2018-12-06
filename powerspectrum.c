/*
#include <math.h>
#include<stdio.h>
#include<stdlib.h>
#include<sys/stat.h>
#include<sys/types.h>
#include<unistd.h>
#include<fcntl.h>
#include<sys/stat.h>

void autocorrelation(double* sdata,double* autocorr,int start,int framelen);
void window(short* sdata,double* windowed_sdata,int start,int framelen);

void autocorrelation(double* sdata,double* autocorr,int start,int framelen){
  //fprintf(stderr, "in autocorrelation\n");
  for (int r = 0; r < framelen; r++){
    for (int t = 0; t < framelen; t++){
      autocorr[r+start] += sdata[t]*sdata[(t+r)%framelen];
    }
  }
    /* normalize
  double r0=autocorr[start];
  for (int i = 0; i < framelen; i++){
    autocorr[i+start]=autocorr[i+start]/r0;
  }
}

void window(short* sdata,double* windowed_sdata,int start,int framelen){
  //fprintf(stderr, "in window\n");
  for (int i = 0 ; i <framelen; i++)
    {
      windowed_sdata[i] = (0.54-0.46*cos(2*M_PI*i/(framelen-1))) *sdata[i+start];
    }
}

int main(int argc, char** argv){
  if(argc != 3 ){
    fprintf(stderr, "Usage: %s DATfile Framelenth \n", argv[0]);
    exit(1);
  }

  FILE* fin;
  int filesize;
  int fd,framelen,framenum;
  struct stat stbuf;
  if((fin=fopen(argv[1],"rb"))==NULL) exit(1);
  if( ( framelen = atoi( argv[2] ) ) < 0 )  exit( 1 );
  if( framelen%2 !=0 ){
    fprintf(stderr, "Please enter an even frame length\n");
  }

  /*ファイルのサイズを取得
  fd=open(argv[1],O_RDONLY);
  fstat(fd,&stbuf);
  filesize = stbuf.st_size/2;
  framenum=filesize/framelen;//the total number of frames
  filesize = filesize/framelen*framelen;//ignore the last part of file if its size is smaller than frame length

  short *sdata = (short*)calloc(filesize,sizeof(short));
  //int *f0 = (int*)calloc(framenum,sizeof(int));
  double *windowed_sdata = (double*)calloc(framelen,sizeof(double));
  double *autocorr = (double*)calloc(filesize,sizeof(double));

  fread( sdata, sizeof(short), filesize*2, fin);

  //fprintf(stderr,"%d\n",filesize);
  //fprintf(stderr,"%d\n",framenum);
  fprintf(stderr, "%d\n",framenum);
  int start=0;
  int i;
  for(i = 0; i < framenum; i++){
    //fprintf(stderr,"%d\n",start);
    window(sdata,windowed_sdata,start,framelen);
    autocorrelation(windowed_sdata,autocorr,start,framelen);
    //fprintf(stderr,"%d\n" ,f0[i]);
    start+=framelen;
  }
  //fprintf(stderr, "%d\n",i );
  //fprintf(stderr, "%d\n",framenum);
  for(i=0 ; i < filesize ; i++){
    printf("%f\n",autocorr[i]);
  }

  fclose(fin);
  return 0;
}
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void FFT (double *xr, double *xi, double *Xr, double *Xi, int N)
{
  int i,j,k,n,n2;
  double theta, wr, wi;

  static double *rbuf, *ibuf;
  static int bufsize = 0;

  /* memory allocation for buffers */
  if ( bufsize != N )
    {
      bufsize = N;
      rbuf = (double*)calloc(  bufsize, sizeof(double) );
      ibuf = (double*)calloc(  bufsize, sizeof(double)  );
    }

  /* bit reverse of xr[] & xi[] --> store to rbuf[] and ibuf[] */
  i = j = 0 ;
  rbuf[j] = xr[j];  ibuf[j] = xi[j];
  for( j = 1 ; j < N-1 ; j++ )
    {
      for( k = N/2 ; k <= i ; k /= 2 )  i -= k;
      i += k;
      rbuf[j] = xr[i];  ibuf[j] = xi[i];
    }
  rbuf[j] = xr[j];  ibuf[j] = xi[j];

   /* butterfly calculation */
  theta = -2.0*M_PI;
  for( n = 1 ; ( n2 = n*2 ) <= N ; n = n2 )
    {
      theta *= 0.5;
      for ( i = 0 ; i < n ; i++ )
	{
	  wr = cos(theta*i);  wi = sin(theta*i);
	  for ( j = i ; j < N ; j += n2 )
	    {
	      k = j + n;
	      /*
		X[j] = 1*buf[j] + W*buf[k];
		X[k] = 1*buf[j] - W*buf[k];
		Note : X[], buf[], and W are complex numbers.
		Re{ X[n] } = Xr[n], Im{ X[n] } = Xi[n];
		Re{ buf[n] } = rbuf[n], Im{ buf[n] } = ibuf[n];
		Re{ W } = wr, Im{ W } = wi;
	      */
	      //Xr[j] = rbuf[j] + 0;  /* ??????????, using wr, wi, rbuf, and ibuf */
        Xr[j] = rbuf[j] + wr*rbuf[k] - wi*ibuf[k];
	      //Xi[j] = ibuf[j] + 0;  /* ??????????, using wr, wi, rbuf, and ibuf */
        Xi[j] = ibuf[j] + wr*ibuf[k] + wi*rbuf[k];
	      //Xr[k] = rbuf[j] - 0;  /* ??????????, using wr, wi, rbuf, and ibuf */
        Xr[k] = rbuf[j] - wr*rbuf[k] + wi*ibuf[k];
	      //Xi[k] = ibuf[j] - 0;  /* ??????????, using wr, wi, rbuf, and ibuf */
        Xi[k] = ibuf[j] - wr*ibuf[k] - wi*rbuf[k];
	    }
	}
      for( i = 0 ; i < N ; i++ )
	{
	  rbuf[i] = Xr[i];
	  ibuf[i] = Xi[i];
	}
    }
  return;
}

void DFT (double *xr, double *xi, double *Xr, double *Xi, int N)
{
  //fprintf(stderr,"DFT\n");
  int k,n;
  for (k = 0 ; k < N ; k++ )
    {
      Xr[k] = Xi[k] = 0.0; // initialization
      for (n = 0 ; n < N ; n++)
	{
	  double wr = cos ( 2.0 * M_PI * n * k / N);
	  double wi = sin ( 2.0 * M_PI * n * k / N);
	  //Xr[k] += 0; /* should be modified. (hint using wr & wi) */
	  //Xi[k] += 0; /* should be modified. (hint using wr & wi) */
    Xr[k] += xr[n]*wr + xi[n]*wi;
    Xi[k] += (-1)*xr[n]*wi + xi[n]*wr;
	}
    }
  return;
}

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
  double *autocorrelationi = (double*)calloc(framelen,sizeof(double));
  double *Sr = (double*)calloc(framelen,sizeof(double));
  double *Si = (double*)calloc(framelen,sizeof(double));
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
    for (int t = 0; t < framelen ;t++){
      //autocorrelation[r] += sdata[t]*sdata[t+r];
      autocorrelation[r] += xr[t]*xr[(t+r)%framelen];
    }
    autocorrelationi[r] = 0;
  }

  /* normalize */
  /*
  double r0=autocorrelation[0];
  for (int i = 0 ; i< framelen ; i++){
    autocorrelation[i]/=r0;
  }
  */

  FFT(autocorrelation,autocorrelationi,Sr,Si,framelen);

  for (int i = 0; i < framelen; i++){
    printf("%lf\n",log10(Sr[i]));
    //printf("%lf\n",autocorrelation[i]);
  }

  return 0;
}
