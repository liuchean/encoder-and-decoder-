#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <sys/stat.h>
#include <sys/types.h>
#include<stdint.h>
/////////////////////這個header是給hoouses.bmp儲存他的header用的//////////
typedef struct Bmpheader {
	char identifier[2]; // 0x0000
	unsigned int filesize; // 0x0002
	unsigned short reserved; // 0x0006
	unsigned short reserved2;
	unsigned int bitmap_dataoffset; // 0x000A
	unsigned int bitmap_headersize; // 0x000E
	unsigned int width; // 0x0012
	unsigned int height; // 0x0016
	unsigned short planes; // 0x001A
	unsigned short bits_perpixel; // 0x001C
	unsigned int compression; // 0x001E
	unsigned int bitmap_datasize; // 0x0022
	unsigned int hresolution; // 0x0026
	unsigned int vresolution; // 0x002A
	unsigned int usedcolors; // 0x002E
	unsigned int importantcolors; // 0x0032
	unsigned int palette; // 0x0036
} Bitmap;

///////////////////儲存輸入檔案的data/////////////////////
typedef struct Gray{
	int R;
	int G;
	int B;
} ImgGray;

/////////////////用於儲存YCbCr///////////////////////
typedef struct _YCbCr{
	double Y;
	double Cb;
	double Cr;

  short QY;
  short QCb;
  short QCr;

  float eY;
	float eCb;
	float eCr;

  int rY;
  int rCb;
  int rCr;

  double rrY;
  double rrCb;
  double rrCr;
}YCbCr;

////////////fast/////////////////////////////////////////
YCbCr** mal_2D(YCbCr **array,int x,int y){
	array=(YCbCr**)malloc(sizeof(YCbCr*)*x);
	for(int i=0;i<x;i++){
		array[i]=(YCbCr*)malloc(sizeof(YCbCr)*y);
	}
	return array;
}
//創建一個函數用於free
void free_2D(YCbCr **array,int x){
	for(int i=0;i<x;i++)
		free(array[i]);
	free(array);
}
YCbCr* mal_1D(YCbCr *array,int x){
	array=(YCbCr*)malloc(sizeof(YCbCr)*x);
	return array;
}
////////////////////////////main function//////////////////////////////////
int main(int argc, char **argv){
////////////////////////////////////////////////////////////////方法0/////////////////////////////////////////////
  if(argc==7){
    char *qqq=argv[1];
    FILE *outfile;//輸入
	outfile = fopen(argv[2], "wb");//output輸出bmp file
   Bitmap header;
        header.identifier[0]='B';
        header.identifier[1]='M';
        header.filesize=36578358;
        header.reserved=0;
        header.reserved2=0;
        header.bitmap_dataoffset=54;
        header.bitmap_headersize=40;
        header.width=3024;
        header.height=4032;
        header.planes=1;
        header.bits_perpixel=24;
        header.compression=0; // 0x001E
	      header.bitmap_datasize=0; // 0x0022
        header.hresolution=3780;//x
        header.vresolution=3780;//y
        header.usedcolors=0;
        header.importantcolors=0;
/////////////////////////////////////////////////////////////////////////	
	ImgGray **Data_Gray = (ImgGray**)malloc(sizeof(ImgGray*)*header.height);

	for(int i=0; i<header.height ;i++)
        Data_Gray[i]=(ImgGray*)malloc(sizeof(ImgGray)*header.width);
        FILE *Rread=fopen(argv[3],"r");
        FILE *Gread=fopen(argv[4],"r");
        FILE *Bread=fopen(argv[5],"r");
        for(int i=0;i<header.height;i++){
		for(int j=0;j<header.width;j++){
			fscanf(Rread,"%d ",&Data_Gray[i][j].R);
			fscanf(Gread,"%d ",&Data_Gray[i][j].G);
			fscanf(Bread,"%d ",&Data_Gray[i][j].B);
		}
	}
  //用來提取H跟W的值
  FILE *dimread=fopen(argv[6],"r");
  int th,tw;
  fscanf(dimread,"H: %d\n",&th);
  fscanf(dimread,"W: %d",&tw);
  printf("%d %d",th,tw);

  fwrite(&header.identifier, sizeof(short), 1, outfile);
	fwrite(&header.filesize, sizeof(int), 1, outfile);
	fwrite(&header.reserved, sizeof(short), 1, outfile);
	fwrite(&header.reserved2, sizeof(short), 1, outfile);
	fwrite(&header.bitmap_dataoffset, sizeof(int), 1, outfile);
	fwrite(&header.bitmap_headersize, sizeof(int), 1, outfile);
	fwrite(&header.width, sizeof(int), 1, outfile);
	fwrite(&header.height, sizeof(int), 1, outfile);
	fwrite(&header.planes, sizeof(short), 1, outfile);
	fwrite(&header.bits_perpixel, sizeof(short), 1, outfile);
	fwrite(&header.compression, sizeof(int), 1, outfile);
	fwrite(&header.bitmap_datasize, sizeof(int), 1, outfile);
	fwrite(&header.hresolution, sizeof(int), 1, outfile);
	fwrite(&header.vresolution, sizeof(int), 1, outfile);
	fwrite(&header.usedcolors, sizeof(int), 1, outfile);
	fwrite(&header.importantcolors, sizeof(int), 1, outfile);

	for (int x = 0; x<header.height; x++){
		for (int y = 0; y<header.width; y++){
			fwrite(&Data_Gray[x][y].B, sizeof(char), 1, outfile);
			fwrite(&Data_Gray[x][y].G, sizeof(char), 1, outfile);
			fwrite(&Data_Gray[x][y].R, sizeof(char), 1, outfile);
		}
	}
	
	for(int i=0;i<th;i++)
		free(Data_Gray[i]);
	free(Data_Gray);
	fclose(outfile);

  printf("完成方法0\n");
	return 0;
  }
////////////////////////////////////////////////////////////////方法1/////////////////////////////////////////////
  else if(argc==10){
  char *qqq1=argv[1];
	FILE *outfile;//輸入
	outfile = fopen(argv[2], "wb");//output輸出bmp file
  
  //頭檔資料
   Bitmap header;
        header.identifier[0]='B';
        header.identifier[1]='M';
        header.filesize=36578358;
        header.reserved=0;
        header.reserved2=0;
        header.bitmap_dataoffset=54;
        header.bitmap_headersize=40;
        header.width=3024;
        header.height=4032;
        header.planes=1;
        header.bits_perpixel=24;
        header.compression=0; // 0x001E
	    header.bitmap_datasize=0; // 0x0022
        header.hresolution=3780;//x
        header.vresolution=3780;//y
        header.usedcolors=0;
        header.importantcolors=0;

  int H = header.height;
	int W = header.width;
  ImgGray **Data_Gray = (ImgGray**)malloc(sizeof(ImgGray*)*H);
	for(int i=0; i<H ;i++)
        Data_Gray[i]=(ImgGray*)malloc(sizeof(ImgGray)*W);
/////////////////////////////////////////////////////////////////////////
//generate basis vector for 2D-DCT
//By matlab formula
  double basic[8][8][8][8];
	for(int u=0; u<8; u++){
		for(int v=0; v<8; v++){

		uint8_t x[8][8];
			for(int r=0; r<8; r++){
				for(int c=0; c<8; c++){
					basic[u][v][r][c] = cos(u*M_PI*(2*r+1)/16)*cos(v*M_PI*(2*c+1)/16);
					x[r][c]= (uint8_t)(127.5+127.5*basic[u][v][r][c]);
				}
				
			}							
		}
	}

/////////////////////////////////////////////Iquantize////////////////////////////////////////////////
  YCbCr **F=mal_2D(F,H,W);
	YCbCr **quan=mal_2D(quan,H,W);
	double tempY=0;
	double tempCb=0;
	double tempCr=0;
	double beta[8];
  //根據matlab寫入beta資料
	for(int i=0; i<8; i++){//C[0] = (1/√2) , beat[1]~beat[7]為1
		if(i==0){ beta[i] = 1.0/sqrt(2); }
		else{ beta[i] = 1.0; }
	}
  
  double rQ[8][8];
  double rQ2[8][8];
  FILE *QYr=fopen(argv[3],"r");
	FILE *QCbr=fopen(argv[4],"r");
	FILE *QCrr=fopen(argv[5],"r");

  FILE *dimread=fopen(argv[6],"r");
  int th,tw;
  fscanf(dimread,"H: %d\n",&th);
  fscanf(dimread,"W: %d",&tw);
  printf("印出提取的H=%d W=%d\n",th,tw);

  H=th;
  W=tw;

	FILE *qY1=fopen(argv[7],"rb");
	FILE *qCb1=fopen(argv[8],"rb");
	FILE *qCr1=fopen(argv[9],"rb");

  FILE *eYr=fopen("eF_Y.raw","rb");
	FILE *eCbr=fopen("eF_Cb.raw","rb");
	FILE *eCrr=fopen("eF_Cr.raw","rb");

  printf("開始讀取量化表\n");
  for(int i=0;i<8;i++){
				for(int j=0;j<8;j++){
      fscanf(QYr,"%lf ",&rQ[i][j]);
			fscanf(QCbr,"%lf ",&rQ2[i][j]);
			fscanf(QCrr,"%lf ",&rQ2[i][j]);
        if(j==7){
      fscanf(QYr,"%lf \n",&rQ[i][j]);
			fscanf(QCbr,"%lf \n",&rQ2[i][j]);
			fscanf(QCrr,"%lf \n",&rQ2[i][j]);

            }
				}
			}
  printf("讀取完成量化表\n");
  printf("開始讀取資料\n");
	for(int i=0;i<H/8;i++){
		for(int j=0;j<W/8;j++){
			for(int u=0;u<8;u++){
				for(int v=0;v<8;v++){
  				fread(&quan[i*8+u][j*8+v].QY,sizeof(short),1,qY1);
				  fread(&quan[i*8+u][j*8+v].QCb,sizeof(short),1,qCb1);
				  fread(&quan[i*8+u][j*8+v].QCr,sizeof(short),1,qCr1);

				}
			}
		}
	}

          
	for(int i=0;i<H/8;i++){
		for(int j=0;j<W/8;j++){
			for(int u=0;u<8;u++){
				for(int v=0;v<8;v++){
					F[i*8+u][j*8+v].Y=quan[i*8+u][j*8+v].QY*rQ[u][v];
					F[i*8+u][j*8+v].Cb=quan[i*8+u][j*8+v].QCb*rQ2[u][v];
					F[i*8+u][j*8+v].Cr=quan[i*8+u][j*8+v].QCr*rQ2[u][v];
				}
			}
		}
	}
  double errorY=0,pY=0,p0Y=0;
  double errorCb=0,pCb=0,p0Cb=0;
  double errorCr=0,pCr=0,p0Cr=0;
  int count=0;
  for(int i=0;i<H/8;i++){
		for(int j=0;j<W/8;j++){
			for(int u=0;u<8;u++){
				for(int v=0;v<8;v++){
        errorY=quan[i*8+u][j*8+v].QY-F[i*8+u][j*8+v].Y;
        pY=pY+errorY*errorY;
        p0Y=p0Y+(quan[i*8+u][j*8+v].QY*quan[i*8+u][j*8+v].QY);

        errorCb=quan[i*8+u][j*8+v].QCb-F[i*8+u][j*8+v].Cb;
        pCb=pCb+errorCb*errorCb;
        p0Cb=p0Cb+(quan[i*8+u][j*8+v].QCb*quan[i*8+u][j*8+v].QCb);

        errorCr=quan[i*8+u][j*8+v].QCr-F[i*8+u][j*8+v].Cr;
        pCr=pCr+errorCr*errorCr;
        p0Cr=p0Cr+(quan[i*8+u][j*8+v].QCr*quan[i*8+u][j*8+v].QCr);
        count++;
				}
			}
		}
	}
  pY=pY/count;
  p0Y=p0Y/count;
  pCb=pCb/count;
  p0Cb=p0Cb/count;
  pCr=pCr/count;
  p0Cr=p0Cr/count;

  //sqnr值越大越好
  double sqnrY=0;
  sqnrY=10*log10(pY/p0Y);
  printf("%f db\n",sqnrY);
  double sqnrCb=0;
  sqnrCb=10*log10(pCb/p0Cb);
  printf("%f db\n",sqnrCb);
  double sqnrCr=0;
  sqnrCr=10*log10(pCr/p0Cr);
  printf("%f db\n",sqnrCr);


	free_2D(quan,H);
	fclose(qY1);
	fclose(qCb1);
	fclose(qCr1);
////////////////////////////////////IDCT//////////////////////////////////
//根據matlab進行inverse DCT
	tempY=0;
	tempCb=0;
	tempCr=0;
	YCbCr **f=mal_2D(f,H,W);
	for(int i=0;i<H/8;i++){
		for(int j=0;j<W/8;j++){
			for(int r=0;r<8;r++){
				for(int c=0;c<8;c++){
					for(int u=0;u<8;u++){
						for(int v=0;v<8;v++){
							tempY=tempY+0.25*beta[u]*beta[v]*F[i*8+u][j*8+v].Y*basic[u][v][r][c];
							tempCb=tempCb+0.25*beta[u]*beta[v]*F[i*8+u][j*8+v].Cb*basic[u][v][r][c];
							tempCr=tempCr+0.25*beta[u]*beta[v]*F[i*8+u][j*8+v].Cr*basic[u][v][r][c];
						}
					}
					f[i*8+r][j*8+c].Y=tempY+128;
					tempY=0;
					f[i*8+r][j*8+c].Cb=tempCb;
					tempCb=0;
					f[i*8+r][j*8+c].Cr=tempCr;
					tempCr=0;
				}
			}
		}
	}
	free_2D(F,H);//free
///////////////////////////////////////IYCbCr/////////////////////////////////
	Data_Gray = (ImgGray**)malloc(sizeof(ImgGray*)*H);
	for(int i=0; i<H ;i++)
        Data_Gray[i]=(ImgGray*)malloc(sizeof(ImgGray)*W);
             
	for(int i=0;i<H;i++){
		for(int j=0;j<W;j++){
			double y=f[i][j].Y;
			double cb=f[i][j].Cb;
			double cr=f[i][j].Cr;
			
			double tmpR=y+1.403*cr;
			double tmpG=y-0.344*cb-0.714*cr;
			double tmpB=y+1.773*cb;
			
			tmpR=round(tmpR);
			tmpG=round(tmpG);
			tmpB=round(tmpB);
			
			
      //把超過255或低於0的數字變成255或0
      //也就是避免overflow
      //再將資料寫回去             
                        //R
			if(tmpR>255)
				Data_Gray[i][j].R=255;
			else if(tmpR<0)
				Data_Gray[i][j].R=0;
			else
				Data_Gray[i][j].R=(unsigned char)tmpR;
			//G	
			if(tmpG>255)
				Data_Gray[i][j].G=255;
			else if(tmpG<0)
				Data_Gray[i][j].G=0;
			else
				Data_Gray[i][j].G=(unsigned char)tmpG;
			//B	
			if(tmpB>255)
				Data_Gray[i][j].B=255;
			else if(tmpB<0)
				Data_Gray[i][j].B=0;
			else
				Data_Gray[i][j].B=(unsigned char)tmpB;

		}
	}
	free_2D(f,H);//free
	//output_bmp(Data_Gray, fp_out, bmpheader, skip);//bmpdata,output file,header,skip

  fwrite(&header.identifier, sizeof(short), 1, outfile);
	fwrite(&header.filesize, sizeof(int), 1, outfile);
	fwrite(&header.reserved, sizeof(short), 1, outfile);
	fwrite(&header.reserved2, sizeof(short), 1, outfile);
	fwrite(&header.bitmap_dataoffset, sizeof(int), 1, outfile);
	fwrite(&header.bitmap_headersize, sizeof(int), 1, outfile);
	fwrite(&header.width, sizeof(int), 1, outfile);
	fwrite(&header.height, sizeof(int), 1, outfile);
	fwrite(&header.planes, sizeof(short), 1, outfile);
	fwrite(&header.bits_perpixel, sizeof(short), 1, outfile);
	fwrite(&header.compression, sizeof(int), 1, outfile);
	fwrite(&header.bitmap_datasize, sizeof(int), 1, outfile);
	fwrite(&header.hresolution, sizeof(int), 1, outfile);
	fwrite(&header.vresolution, sizeof(int), 1, outfile);
	fwrite(&header.usedcolors, sizeof(int), 1, outfile);
	fwrite(&header.importantcolors, sizeof(int), 1, outfile);

	for (int x = 0; x<header.height; x++){
		for (int y = 0; y<header.width; y++){
			fwrite(&Data_Gray[x][y].B, sizeof(char), 1, outfile);
			fwrite(&Data_Gray[x][y].G, sizeof(char), 1, outfile);
			fwrite(&Data_Gray[x][y].R, sizeof(char), 1, outfile);
		}
	}
	
	for(int i=0;i<4032;i++)
		free(Data_Gray[i]);
	free(Data_Gray);
	fclose(outfile);
  printf("完成方法1\n");
	return 0;
  }
  ////////////////////////////////////////////////////////////////方法1(b)/////////////////////////////////////////////
  else if(argc==13){
  
  char *qqq1=argv[1];
	FILE *outfile;//輸入
	outfile = fopen(argv[2], "wb");//output輸出bmp file

   Bitmap header;
        header.identifier[0]='B';
        header.identifier[1]='M';
        header.filesize=36578358;
        header.reserved=0;
        header.reserved2=0;
        header.bitmap_dataoffset=54;
        header.bitmap_headersize=40;
        header.width=3024;
        header.height=4032;
        header.planes=1;
        header.bits_perpixel=24;
        header.compression=0; // 0x001E
	      header.bitmap_datasize=0; // 0x0022
        header.hresolution=3780;//x
        header.vresolution=3780;//y
        header.usedcolors=0;
        header.importantcolors=0;

  int H = header.height;
	int W = header.width;
  ImgGray **Data_Gray = (ImgGray**)malloc(sizeof(ImgGray*)*H);
	for(int i=0; i<H ;i++)
        Data_Gray[i]=(ImgGray*)malloc(sizeof(ImgGray)*W);
/////////////////////////////////////////////////////////////////////////
//generate basis vector for 2D-DCT
//By matlab formula
  double basic[8][8][8][8];
	for(int u=0; u<8; u++){
		for(int v=0; v<8; v++){

		uint8_t x[8][8];
			for(int r=0; r<8; r++){
				for(int c=0; c<8; c++){
					basic[u][v][r][c] = cos(u*M_PI*(2*r+1)/16)*cos(v*M_PI*(2*c+1)/16);
					x[r][c]= (uint8_t)(127.5+127.5*basic[u][v][r][c]);
				}
				
			}							
		}
	}

/////////////////////////////////////////////Iquantize////////////////////////////////////////////////
  YCbCr **F=mal_2D(F,H,W);
	YCbCr **quan=mal_2D(quan,H,W);
	double tempY=0;
	double tempCb=0;
	double tempCr=0;
	double beta[8];
  //根據matlab寫入beta資料
	for(int i=0; i<8; i++){//C[0] = (1/√2) , beat[1]~beat[7]為1
		if(i==0){ beta[i] = 1.0/sqrt(2); }
		else{ beta[i] = 1.0; }
	}
  
  double rQ[8][8];
  double rQ2[8][8];
  FILE *QYr=fopen(argv[3],"r");
	FILE *QCbr=fopen(argv[4],"r");
	FILE *QCrr=fopen(argv[5],"r");

  FILE *dimread=fopen(argv[6],"r");
  int th,tw;
  fscanf(dimread,"H: %d\n",&th);
  fscanf(dimread,"W: %d",&tw);
  printf("印出提取的H=%d W=%d\n",th,tw);
  fclose(dimread);

  H=th;
  W=tw;

	FILE *qY1=fopen(argv[7],"rb");
	FILE *qCb1=fopen(argv[8],"rb");
	FILE *qCr1=fopen(argv[9],"rb");

  for(int i=0;i<8;i++){
				for(int j=0;j<8;j++){
      fscanf(QYr,"%lf ",&rQ[i][j]);
			fscanf(QCbr,"%lf ",&rQ2[i][j]);
			fscanf(QCrr,"%lf ",&rQ2[i][j]);
        if(j==7){
      fscanf(QYr,"%lf \n",&rQ[i][j]);
			fscanf(QCbr,"%lf \n",&rQ2[i][j]);
			fscanf(QCrr,"%lf \n",&rQ2[i][j]);

            }
				}
			}
  printf("完成陣列的寫入");

 
	for(int i=0;i<H/8;i++){
		for(int j=0;j<W/8;j++){
			for(int u=0;u<8;u++){
				for(int v=0;v<8;v++){
  				fread(&quan[i*8+u][j*8+v].QY,sizeof(short),1,qY1);
				  fread(&quan[i*8+u][j*8+v].QCb,sizeof(short),1,qCb1);
				  fread(&quan[i*8+u][j*8+v].QCr,sizeof(short),1,qCr1);

				}
			}
		}
	}
  FILE *eYr=fopen(argv[10],"rb");
	FILE *eCbr=fopen(argv[11],"rb");
	FILE *eCrr=fopen(argv[12],"rb");
  printf("完成quan的寫入\n");
  for(int i=0;i<H/8;i++){
		for(int j=0;j<W/8;j++){
			for(int u=0;u<8;u++){
				for(int v=0;v<8;v++){
  				fread(&quan[i*8+u][j*8+v].eY,sizeof(float),1,eYr);
				  fread(&quan[i*8+u][j*8+v].eCb,sizeof(float),1,eCbr);
				  fread(&quan[i*8+u][j*8+v].eCr,sizeof(float),1,eCrr);

				}
			}
		}
	}
  printf("完成e的寫入\n");


          
	for(int i=0;i<H/8;i++){
		for(int j=0;j<W/8;j++){
			for(int u=0;u<8;u++){
				for(int v=0;v<8;v++){
					F[i*8+u][j*8+v].Y=quan[i*8+u][j*8+v].eY+quan[i*8+u][j*8+v].QY*rQ[u][v];
					F[i*8+u][j*8+v].Cb=quan[i*8+u][j*8+v].eCb+quan[i*8+u][j*8+v].QCb*rQ2[u][v];
					F[i*8+u][j*8+v].Cr=quan[i*8+u][j*8+v].eCr+quan[i*8+u][j*8+v].QCr*rQ2[u][v];
				}
			}
		}
	}

  printf("完成1b的圖形復原\n");
  


	free_2D(quan,H);

  fclose(QYr);
  fclose(QCbr);
  fclose(QCrr);
  fclose(eYr);
  fclose(eCbr);
  fclose(eCrr);
	fclose(qY1);
	fclose(qCb1);
	fclose(qCr1);
////////////////////////////////////IDCT//////////////////////////////////
//根據matlab進行inverse DCT
	tempY=0;
	tempCb=0;
	tempCr=0;
	YCbCr **f=mal_2D(f,H,W);
	for(int i=0;i<H/8;i++){
		for(int j=0;j<W/8;j++){
			for(int r=0;r<8;r++){
				for(int c=0;c<8;c++){
					for(int u=0;u<8;u++){
						for(int v=0;v<8;v++){
							tempY=tempY+0.25*beta[u]*beta[v]*F[i*8+u][j*8+v].Y*basic[u][v][r][c];
							tempCb=tempCb+0.25*beta[u]*beta[v]*F[i*8+u][j*8+v].Cb*basic[u][v][r][c];
							tempCr=tempCr+0.25*beta[u]*beta[v]*F[i*8+u][j*8+v].Cr*basic[u][v][r][c];
						}
					}
					f[i*8+r][j*8+c].Y=tempY+128;
					tempY=0;
					f[i*8+r][j*8+c].Cb=tempCb;
					tempCb=0;
					f[i*8+r][j*8+c].Cr=tempCr;
					tempCr=0;
				}
			}
		}
	}
  printf("完成IDCT");
	free_2D(F,H);//free
///////////////////////////////////////IYCbCr/////////////////////////////////
	Data_Gray = (ImgGray**)malloc(sizeof(ImgGray*)*H);
	for(int i=0; i<H ;i++)
        Data_Gray[i]=(ImgGray*)malloc(sizeof(ImgGray)*W);
             
	for(int i=0;i<H;i++){
		for(int j=0;j<W;j++){
			double y=f[i][j].Y;
			double cb=f[i][j].Cb;
			double cr=f[i][j].Cr;
			
			double tmpR=y+1.403*cr;
			double tmpG=y-0.344*cb-0.714*cr;
			double tmpB=y+1.773*cb;
			
			tmpR=round(tmpR);
			tmpG=round(tmpG);
			tmpB=round(tmpB);
			
			
      //把超過255或低於0的數字變成255或0
      //也就是避免overflow
      //再將資料寫回去             
                        //R
			if(tmpR>255)
				Data_Gray[i][j].R=255;
			else if(tmpR<0)
				Data_Gray[i][j].R=0;
			else
				Data_Gray[i][j].R=(unsigned char)tmpR;
			//G	
			if(tmpG>255)
				Data_Gray[i][j].G=255;
			else if(tmpG<0)
				Data_Gray[i][j].G=0;
			else
				Data_Gray[i][j].G=(unsigned char)tmpG;
			//B	
			if(tmpB>255)
				Data_Gray[i][j].B=255;
			else if(tmpB<0)
				Data_Gray[i][j].B=0;
			else
				Data_Gray[i][j].B=(unsigned char)tmpB;

		}
	}
  printf("完成IYCbCr");
	free_2D(f,H);//free
	//output_bmp(Data_Gray, fp_out, bmpheader, skip);//bmpdata,output file,header,skip
  printf("進行檔案寫入");
  fwrite(&header.identifier, sizeof(short), 1, outfile);
	fwrite(&header.filesize, sizeof(int), 1, outfile);
	fwrite(&header.reserved, sizeof(short), 1, outfile);
	fwrite(&header.reserved2, sizeof(short), 1, outfile);
	fwrite(&header.bitmap_dataoffset, sizeof(int), 1, outfile);
	fwrite(&header.bitmap_headersize, sizeof(int), 1, outfile);
	fwrite(&header.width, sizeof(int), 1, outfile);
	fwrite(&header.height, sizeof(int), 1, outfile);
	fwrite(&header.planes, sizeof(short), 1, outfile);
	fwrite(&header.bits_perpixel, sizeof(short), 1, outfile);
	fwrite(&header.compression, sizeof(int), 1, outfile);
	fwrite(&header.bitmap_datasize, sizeof(int), 1, outfile);
	fwrite(&header.hresolution, sizeof(int), 1, outfile);
	fwrite(&header.vresolution, sizeof(int), 1, outfile);
	fwrite(&header.usedcolors, sizeof(int), 1, outfile);
	fwrite(&header.importantcolors, sizeof(int), 1, outfile);

	for (int x = 0; x<header.height; x++){
		for (int y = 0; y<header.width; y++){
			fwrite(&Data_Gray[x][y].B, sizeof(char), 1, outfile);
			fwrite(&Data_Gray[x][y].G, sizeof(char), 1, outfile);
			fwrite(&Data_Gray[x][y].R, sizeof(char), 1, outfile);
		}
	}
	
	for(int i=0;i<4032;i++)
		free(Data_Gray[i]);
	free(Data_Gray);
	fclose(outfile);
  printf("完成方法1(b)全部過程\n");
	return 0;
  }
////////////////////////////////////////////////////////////////////法2///////////////////////////////////////////////////////////
else if(argc==5){

  char *qqq1=argv[1];
	FILE *outfile;//輸入
	outfile = fopen(argv[2], "wb");//output輸出bmp file

   Bitmap header;
        header.identifier[0]='B';
        header.identifier[1]='M';
        header.filesize=36578358;
        header.reserved=0;
        header.reserved2=0;
        header.bitmap_dataoffset=54;
        header.bitmap_headersize=40;
        header.width=3024;
        header.height=4032;
        header.planes=1;
        header.bits_perpixel=24;
        header.compression=0; // 0x001E
	      header.bitmap_datasize=0; // 0x0022
        header.hresolution=3780;//x
        header.vresolution=3780;//y
        header.usedcolors=0;
        header.importantcolors=0;

  int H = header.height;
	int W = header.width;
  ImgGray **Data_Gray = (ImgGray**)malloc(sizeof(ImgGray*)*H);
	for(int i=0; i<H ;i++)
        Data_Gray[i]=(ImgGray*)malloc(sizeof(ImgGray)*W);
/////////////////////////////////////////////////////////////////////////
//generate basis vector for 2D-DCT
//By matlab formula
  double basic[8][8][8][8];
	for(int u=0; u<8; u++){
		for(int v=0; v<8; v++){

		uint8_t x[8][8];
			for(int r=0; r<8; r++){
				for(int c=0; c<8; c++){
					basic[u][v][r][c] = cos(u*M_PI*(2*r+1)/16)*cos(v*M_PI*(2*c+1)/16);
					x[r][c]= (uint8_t)(127.5+127.5*basic[u][v][r][c]);
				}
				
			}							
		}
	}
////////////////////////////////inverse RLE////////////////////////////////

  char *qqq2=argv[3];//binary or txt
  int bH=4032;
  int bW=3024;
  YCbCr *zig=mal_1D(zig,bH*bW);
  YCbCr *RLE=mal_1D(RLE,2*bH*bW);
  int big=10000;
  int checknum=2;
  if(strcmp(argv[3],"binary")==0) checknum=0;
  else if (strcmp(argv[3],"ascii")==0) checknum=1;
  
  int indexY=0;
  int indexCb=0;
  int indexCr=0;
  //txt
  if(checknum==1){
  printf("讀到acsii開始進行讀檔\n");
  FILE *ascii=fopen(argv[4], "r");
  for(int x=0;x<H/8;x++){
     for(int y=0;y<W/8;y++){
       ////////////Y//////////
       fscanf(ascii," %d %d Y",&x,&y);
       for(int indexY=0;indexY<64;indexY+=2){
        fscanf(ascii," %d %d",&RLE[indexY].rY,&RLE[indexY+1].rY);
       if(RLE[indexY].rY==0&&RLE[indexY+1].rY==0){
        fscanf(ascii," %d %d\n",&RLE[indexY].rY,&RLE[indexY+1].rY);
        break;
       }
      }
       /////////////////Cb//////////////
       fscanf(ascii," %d %d Cb",&x,&y);
       for(int indexCb=0;indexCb<64;indexCb+=2){
        fscanf(ascii," %d %d",&RLE[indexCb].rCb,&RLE[indexCb+1].rCb);
       if(RLE[indexCb].rCb==0&&RLE[indexCb+1].rCb==0){
        fscanf(ascii," %d %d\n",&RLE[indexCb].rCb,&RLE[indexCb+1].rCb);
        break;
       }
      }
       //////////////////////Cr/////////////////
       fscanf(ascii," %d %d Cr",&x,&y);
       for(int indexCr=0;indexCr<64;indexCr+=2){
        fscanf(ascii," %d %d",&RLE[indexCr].rCr,&RLE[indexCr+1].rCr);
       if(RLE[indexCr].rCr==0&&RLE[indexCr+1].rCr==0){
        fscanf(ascii," %d %d\n",&RLE[indexCr].rCr,&RLE[indexCr+1].rCr);
        break;
       }

       }
     }
  }
  printf("結束讀取txt檔\n");
  fclose(ascii);
  }
  
  //bin
  else if(checknum==0){
  printf("讀取到binary開始進行讀取");
  FILE *binread=fopen(argv[4],"rb");
   for(indexY=0;indexY<big;indexY++){
  fread(&RLE[indexY].rY, sizeof(int), 1, binread);
  if(RLE[indexY].rY==0&&RLE[indexY+1].rY==0) break;
  }

   for(indexCb=0;indexCb<big;indexCb++){
  fread(&RLE[indexCb].rCb, sizeof(int), 1, binread);
     if(RLE[indexCb+1].rCb==0&&RLE[indexCb+1].rCb==0)
    break;
  }

 for(indexCr=0;indexCr<big;indexCr++){
  fread(&RLE[indexCr].rCr, sizeof(int), 1, binread);
  if(RLE[indexCr].rCr==0&&RLE[indexCr+1].rCr==0){
  break;
  }
  }
  printf("完成讀取檔案\n");
  fclose(binread);
  }
  

  int index=0;	
	for(int i=0;i<bH*bW;i+=64){
		for(int j=0;j<64;j++){
      ///首先分為兩個幾種情況
      //前不0後不0=>第一個if處理
      //前0後不0=>經第一個if處理完一定會到第二種情況
      //第三種情況就是encoder所設計的結尾為00,然後就break
			if(RLE[index].rY!=0){  
				RLE[index].rY--;
				zig[i+j].Y=0;
			}
			else if(RLE[index].rY==0&&RLE[index+1].rY!=0){  
				zig[i+j].Y=RLE[index+1].rY;
				index+=2;
			}
      if(RLE[index].rY==0&&RLE[index+1].rY==0) break;
		}
		
	}


	index=0;
	for(int i=0;i<bH*bW;i+=64){
		for(int j=0;j<64;j++){
			if(RLE[index].rCb!=0){
				RLE[index].rCb--;
				zig[i+j].Cb=0;
			}
			else if(RLE[index].rCb==0&&RLE[index+1].rCb!=0){
				zig[i+j].Cb=RLE[index+1].rCb;
				index+=2;
			}
      if(RLE[index].rCb==0&&RLE[index+1].rCb==0) break;
		}
		
	}

	index=0;
	for(int i=0;i<bH*bW;i+=64){
		for(int j=0;j<64;j++){
			if(RLE[index].rCr!=0){
				RLE[index].rCr--;
				zig[i+j].Cr=0;
			}
			else if(RLE[index].rCr==0&&RLE[index+1].rCr!=0){
				zig[i+j].Cr=RLE[index+1].rCr;
				index+=2;
			}
      if(RLE[index].rCr==0&&RLE[index+1].rCr==0) break;
		}
		
	}
	free(RLE);
////////////////////////////////zigzag table////////////////////////////////
	int vec[8][8]={{ 0, 1, 5, 6,14,15,27,28},
				   { 2, 4, 7,13,16,26,29,42},
 				   { 3, 8,12,17,25,30,41,43},
 				   { 9,11,18,24,31,40,44,53},
 				   {10,19,23,32,39,45,52,54},
 				   {20,22,33,38,46,51,55,60},
 				   {21,34,37,47,50,56,59,61},
 				   {35,36,48,49,57,58,62,63}};

////////////////////////////////Izigzag////////////////////////////////
	int count=0;
	YCbCr **DPCM=mal_2D(DPCM,H,W);
	for(int i=0;i<H/8;i++){
		for(int j=0;j<W/8;j++){
			for(int u=0;u<8;u++){
				for(int v=0;v<8;v++){
					DPCM[i*8+u][j*8+v].Y=zig[64*count+vec[u][v]].Y;
					DPCM[i*8+u][j*8+v].Cb=zig[64*count+vec[u][v]].Cb;
					DPCM[i*8+u][j*8+v].Cr=zig[64*count+vec[u][v]].Cr;
				}
			}
			count++;
		}
	}
	free(zig);
////////////////////////////////IDPCM////////////////////////////////
	YCbCr **quan=mal_2D(quan,H,W);

	for(int i=0; i<H; i++){
		for(int j=0; j<W; j++){
			if(i>0 && j==0){
				DPCM[i][j].Y= DPCM[i][j].Y+DPCM[i-1][j].Y;
				DPCM[i][j].Cb = DPCM[i][j].Cb +DPCM[i-1][j].Cb;
				DPCM[i][j].Cr = DPCM[i][j].Cr +DPCM[i-1][j].Cr;
			}
			else if(j>0){
				DPCM[i][j].Y=DPCM[i][j].Y + DPCM[i][j-1].Y;
				DPCM[i][j].Cb= DPCM[i][j].Cb +DPCM[i][j-1].Cb;
				DPCM[i][j].Cr=DPCM[i][j].Cr + DPCM[i][j-1].Cr;
			}
		}
	}
	for(int i=0; i<H; i++){
		for(int j=0; j<W; j++){
      quan[i][j].Y=DPCM[i][j].Y;
			quan[i][j].Cb=DPCM[i][j].Cb;
			quan[i][j].Cr=DPCM[i][j].Cr;

		}
	}
	free_2D(DPCM,H);
///////////////////////////////////////////quantize table////////////////////////////////////////
	double Q_Y[8][8]={{16, 11, 10, 16, 24, 40, 51, 61},
    			      {12, 12, 14, 19, 26, 58, 60, 55},
    			      {14, 13, 16, 24, 40, 57, 69, 56},
    			      {14, 17, 22, 29, 51, 87, 80, 62},
    			      {18, 22, 37, 56, 68,109,103, 77},
    			      {24, 35, 55, 64, 81,104,113, 92},
    			      {49, 64, 78, 87,103,121,120,101},
    			      {72, 92, 95, 98,112,100,103, 99}};


	double Q_CbCr[8][8]={{17,18,24,47,99,99,99,99},
    				    {18,21,26,66,99,99,99,99},
    				    {24,26,56,99,99,99,99,99},
    				    {47,66,99,99,99,99,99,99},
    				    {99,99,99,99,99,99,99,99},
    				    {99,99,99,99,99,99,99,99},
    				    {99,99,99,99,99,99,99,99},
    				    {99,99,99,99,99,99,99,99}};
/////////////////////////////////////////////Iquantize////////////////////////////////////////////////
	//根據matlab進行inverse q
  YCbCr **F=mal_2D(F,H,W);
	//H=512,W=512
	for(int i=0;i<H/8;i++){
		for(int j=0;j<W/8;j++){
			for(int u=0;u<8;u++){
				for(int v=0;v<8;v++){
					F[i*8+u][j*8+v].Y=quan[i*8+u][j*8+v].Y*Q_Y[u][v];
					F[i*8+u][j*8+v].Cb=quan[i*8+u][j*8+v].Cb*Q_CbCr[u][v];
					F[i*8+u][j*8+v].Cr=quan[i*8+u][j*8+v].Cr*Q_CbCr[u][v];

				}
			}
		}
	}
	free_2D(quan,H);
////////////////////////////////////IDCT//////////////////////////////////
//根據matlab進行inverse DCT
	int tempY=0;
	int tempCb=0;
	int tempCr=0;
	YCbCr **f=mal_2D(f,H,W);
  double beta[8];
  //根據matlab寫入beta資料
	for(int i=0; i<8; i++){//C[0] = (1/√2) , beat[1]~beat[7]為1
		if(i==0){ beta[i] = 1.0/sqrt(2); }
		else{ beta[i] = 1.0; }
	}

	for(int i=0;i<H/8;i++){
		for(int j=0;j<W/8;j++){
			for(int r=0;r<8;r++){
				for(int c=0;c<8;c++){
					for(int u=0;u<8;u++){
						for(int v=0;v<8;v++){
							tempY=tempY+0.25*beta[u]*beta[v]*F[i*8+u][j*8+v].Y*basic[u][v][r][c];
							tempCb=tempCb+0.25*beta[u]*beta[v]*F[i*8+u][j*8+v].Cb*basic[u][v][r][c];
							tempCr=tempCr+0.25*beta[u]*beta[v]*F[i*8+u][j*8+v].Cr*basic[u][v][r][c];
						}
					}
					f[i*8+r][j*8+c].Y=tempY+128;
					tempY=0;
					f[i*8+r][j*8+c].Cb=tempCb;
					tempCb=0;
					f[i*8+r][j*8+c].Cr=tempCr;
					tempCr=0;
				}
			}
		}
	}
	free_2D(F,H);//free
///////////////////////////////////////IYCbCr/////////////////////////////////
	Data_Gray = (ImgGray**)malloc(sizeof(ImgGray*)*H);
	for(int i=0; i<H ;i++)
        Data_Gray[i]=(ImgGray*)malloc(sizeof(ImgGray)*W);
             
	for(int i=0;i<H;i++){
		for(int j=0;j<W;j++){
			double y=f[i][j].Y;
			double cb=f[i][j].Cb;
			double cr=f[i][j].Cr;
			double tmpR=y+1.403*cr;
			double tmpG=y-0.344*cb-0.714*cr;
			double tmpB=y+1.773*cb;
			tmpR=round(tmpR);
			tmpG=round(tmpG);
			tmpB=round(tmpB);
	
      //把超過255或低於0的數字變成255或0
      //也就是避免overflow
      //再將資料寫回去

      //R
			if(tmpR>255)
				Data_Gray[i][j].R=255;
			else if(tmpR<0)
				Data_Gray[i][j].R=0;
			else
				Data_Gray[i][j].R=(unsigned char)tmpR;
			//G	
			if(tmpG>255)
				Data_Gray[i][j].G=255;
			else if(tmpG<0)
				Data_Gray[i][j].G=0;
			else
				Data_Gray[i][j].G=(unsigned char)tmpG;
			//B	
			if(tmpB>255)
				Data_Gray[i][j].B=255;
			else if(tmpB<0)
				Data_Gray[i][j].B=0;
			else
				Data_Gray[i][j].B=(unsigned char)tmpB;

		}
	}
	free_2D(f,H);//free
  	//output_bmp(Data_Gray, fp_out, bmpheader, skip);//bmpdata,output file,header,skip
  printf("進行檔案寫入");
  fwrite(&header.identifier, sizeof(short), 1, outfile);
	fwrite(&header.filesize, sizeof(int), 1, outfile);
	fwrite(&header.reserved, sizeof(short), 1, outfile);
	fwrite(&header.reserved2, sizeof(short), 1, outfile);
	fwrite(&header.bitmap_dataoffset, sizeof(int), 1, outfile);
	fwrite(&header.bitmap_headersize, sizeof(int), 1, outfile);
	fwrite(&header.width, sizeof(int), 1, outfile);
	fwrite(&header.height, sizeof(int), 1, outfile);
	fwrite(&header.planes, sizeof(short), 1, outfile);
	fwrite(&header.bits_perpixel, sizeof(short), 1, outfile);
	fwrite(&header.compression, sizeof(int), 1, outfile);
	fwrite(&header.bitmap_datasize, sizeof(int), 1, outfile);
	fwrite(&header.hresolution, sizeof(int), 1, outfile);
	fwrite(&header.vresolution, sizeof(int), 1, outfile);
	fwrite(&header.usedcolors, sizeof(int), 1, outfile);
	fwrite(&header.importantcolors, sizeof(int), 1, outfile);

	for (int x = 0; x<header.height; x++){
		for (int y = 0; y<header.width; y++){
			fwrite(&Data_Gray[x][y].B, sizeof(char), 1, outfile);
			fwrite(&Data_Gray[x][y].G, sizeof(char), 1, outfile);
			fwrite(&Data_Gray[x][y].R, sizeof(char), 1, outfile);
		}
	}
	
	for(int i=0;i<4032;i++)
		free(Data_Gray[i]);
	free(Data_Gray);
	fclose(outfile);
  printf("完成方法2(a)全部過程\n");
	return 0;

}
	
	}//底
