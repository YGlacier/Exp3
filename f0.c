#include <math.h>
#include<stdio.h>
#include<stdlib.h>
#include<sys/stat.h>
#include<sys/types.h>
#include<unistd.h>
#include<fcntl.h>
#include<sys/stat.h>

void autocorrelation(double* sdata,double* autocorr,int framelen);
void window(short* sdata,double* windowed_sdata,int start,int framelen);
int calculatef0samples(double* autocorr, int framelen);

void autocorrelation(double* sdata,double* autocorr,int framelen){
  //fprintf(stderr, "in autocorrelation\n");
  for (int r = 0; r < framelen; r++){
    for (int t = 0; t < framelen - 1 - r ; t++){
      autocorr[r] += sdata[t]*sdata[(t+r)/*%framelen*/];
    }
  }
    /* normalize */
  double r0=autocorr[0];
  for (int i = 0; i < framelen; i++){
    autocorr[i]=autocorr[i]/r0;
  }
}

void window(short* sdata,double* windowed_sdata,int start,int framelen){
  //fprintf(stderr, "in window\n");
  for (int i = 0 ; i <framelen; i++)
    {
      windowed_sdata[i] = (0.54-0.46*cos(2*M_PI*i/(framelen-1))) *sdata[i+start];
    }
}

int calculatef0samples(double* autocorr,int framelen){
  //fprintf(stderr, "in calculatef0\n");
  int i=0;
  for(i = 0 ; ;i++){
    if(autocorr[i] * autocorr[i+1] < 0)
    break;
  }
  int maxindex=0;
  double max=-1;
  for(; i < framelen; i++){
    if(autocorr[i] > max){
      maxindex=i;
      max=autocorr[i];
    }
  }
  //fprintf(stderr, "%d\n",maxindex);
  return maxindex;
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

  /*ファイルのサイズを取得*/
  fd=open(argv[1],O_RDONLY);
  fstat(fd,&stbuf);
  filesize = stbuf.st_size/2;
  framenum=filesize/framelen*2-1;//the total number of frames
  filesize = filesize/framelen*framelen;//ignore the last part of file if its size is smaller than frame length

  short *sdata = (short*)calloc(filesize,sizeof(short));
  int *f0_samples = (int*)calloc(framenum,sizeof(int));
  double *windowed_sdata = (double*)calloc(framelen,sizeof(double));
  double *autocorr = (double*)calloc(framelen,sizeof(double));

  fread( sdata, sizeof(short), filesize*2, fin);

  //fprintf(stderr,"%d\n",filesize);
  //fprintf(stderr,"%d\n",framenum);
  fprintf(stderr, "%d\n",framenum);
  int start=0;
  int i;
  for(i = 0; i < framenum; i++){
    //fprintf(stderr,"%d\n",start);
    window(sdata,windowed_sdata,start,framelen);
    autocorrelation(windowed_sdata,autocorr,framelen);
    f0_samples[i] = calculatef0samples(autocorr,framelen);
    //fprintf(stderr,"%d\n" ,f0[i]);
    start+= (framelen/2);
  }

  double t0=1.0/16000;//サンプリングレートから、隣接する二つのshortの時間間隔を算出
  double *T0 = (double*)calloc(framenum,sizeof(double));

  for ( i = 0 ; i<framenum; i++){
    T0[i] = f0_samples[i] * t0;
  }

  double *logf0 = (double*)calloc(framenum,sizeof(double));

  for (i=0;i<framenum;i++){
    logf0[i] = log10(1/T0[i]);
  }
  for(i=0 ; i < framenum ; i++){
    printf("%lf\n",logf0[i]);
    //printf("%d\n",f0_samples[i]);
  }

  fclose(fin);
  return 0;
}
