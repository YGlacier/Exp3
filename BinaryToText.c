#include<stdio.h>
#include<stdlib.h>
#include<sys/stat.h>
#include<sys/types.h>
#include<unistd.h>
#include<fcntl.h>
#include<sys/stat.h>

int main(int argc,char ** argv){
  if(argc!=2){
    fprintf(stderr,"Wrong Commandline\n");
    exit(1);
  }

  FILE* fin;
  short* data;
  long filesize;
  int fd;
  struct stat stbuf;
  if((fin=fopen(argv[1],"rb"))==NULL) exit(1);

  /*ファイルのサイズを取得（メモリ確保のため）*/
  fd=open(argv[1],O_RDONLY);
  fstat(fd,&stbuf);
  filesize = stbuf.st_size;

  data=(short*)calloc(filesize/2, sizeof(short));
  fread( data, sizeof(short), filesize, fin);
  for(int i = 0;i < filesize/2; i++){
    printf("%d\n",data[i]);
  }
  fclose(fin);
  return 0;
}
