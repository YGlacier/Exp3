#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void autocorrelation(double* sdata,double* autocorr);
void window(short* sdata,double* windowed_sdata,int start,int framelen);
double calculatef0(double* autocorr);

void autocorrelation(double* sdata,double* autocorr,int framelen){
    double *xr = (double*)calloc(framelen,sizeof(double));
    for (int i = 0 ; i < framelen; i++)
      {
        xr[i] = (0.54-0.46*cos(2*M_PI*i/(framelen-1))) *sdata[i];
      }


    /* cauculate autocorrelation */
    for (int r = 0; r < framelen; r++){
      for (int t = 0; t < framelen ;t++){
        //autocorrelation[r] += sdata[t]*sdata[t+r];
        autocorr[r] += xr[t]*xr[(t+r)%framelen];
      }
    }

    /* normalize */
    double r0=autocorr[0];
    for (int i = 0; i < framelen; i++){
      autocorr[i]=autocorr[i]/r0;
    }
}
void window(short* sdata,double* windowed_sdata,int start,int framelen){
  for (int i = start ; i < start + framelen; i++)
    {
      windowed_sdata[i] = (0.54-0.46*cos(2*M_PI*i/(framelen-1))) *sdata[i];
    }
}



void main(int argc, char** argv){
  if(argc != 3 ){
    fprintf(stderr, "Usage: %s DATfile Framelenth \n", argv[0]);
    exit(1);
  }

  FILE* fin;
  short* data;
  long filesize;
  int fd,framelen;
  struct stat stbuf;
  if((fin=fopen(argv[1],"rb"))==NULL) exit(1);
  if( ( framelen = atoi( argv[3] ) ) < 0 )  exit( 1 );
  if( framelen%2 !=0 ){
    fprintf(stderr, "Please enter an even frame length\n");
  }

  /*ファイルのサイズを取得*/
  fd=open(argv[1],O_RDONLY);
  fstat(fd,&stbuf);
  filesize = stbuf.st_size;
  framenum=filesize/framelen-1;//the total number of frames
  filesize = filesize/framelen*framelen;//ignore the last part of file if its size is smaller than frame length

  short *sdata = (short*)calloc(filesize,sizeof(short));
  double *f0 = (double*)calloc(framenum,sizeof(double));
  double *windowed_sdata = (double*)calloc(framelen,sizeof(double));
  double *autocorr = (double*)calloc(framelen,sizeof(double));


  for(int i = 0; i < framenum; i++){
    for(int r = 0; r < framelen; r++){
      for(int t = 0; j < framelen; j++){
        window(sdata,windowed_sdata,i,framelen);
        autocorrelation()
      }
    }
  }
  i+=(framelen/2);
}
