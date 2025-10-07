#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include <memory.h>
#include <sys/stat.h>
#include <sys/types.h>
#include<stdint.h>
///////////////////////table/////////////////////////////////
	double Q[8][8]={{16, 11, 10, 16, 24, 40, 51, 61},
    			      {12, 12, 14, 19, 26, 58, 60, 55},
    			      {14, 13, 16, 24, 40, 57, 69, 56},
    			      {14, 17, 22, 29, 51, 87, 80, 62},
    			      {18, 22, 37, 56, 68,109,103, 77},
    			      {24, 35, 55, 64, 81,104,113, 92},
    			      {49, 64, 78, 87,103,121,120,101},
    			      {72, 92, 95, 98,112,100,103, 99}};


	double Q2[8][8]={{17,18,24,47,99,99,99,99},
    				    {18,21,26,66,99,99,99,99},
    				    {24,26,56,99,99,99,99,99},
    				    {47,66,99,99,99,99,99,99},
    				    {99,99,99,99,99,99,99,99},
    				    {99,99,99,99,99,99,99,99},
    				    {99,99,99,99,99,99,99,99},
    				    {99,99,99,99,99,99,99,99}};
  int vec[8][8]={{ 0, 1, 5, 6,14,15,27,28},
				   { 2, 4, 7,13,16,26,29,42},
 				   { 3, 8,12,17,25,30,41,43},
 				   { 9,11,18,24,31,40,44,53},
 				   {10,19,23,32,39,45,52,54},
 				   {20,22,33,38,46,51,55,60},
 				   {21,34,37,47,50,56,59,61},
 				   {35,36,48,49,57,58,62,63}};

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
///////////////////////讀取houses.bmp的頭檔共 54byte//////////////////
void readheader(FILE* fp, Bitmap *x) {
	fread(&x->identifier, sizeof(x->identifier), 1, fp);
	fread(&x->filesize, sizeof(x->filesize), 1, fp);
	fread(&x->reserved, sizeof(x->reserved), 1, fp);
	fread(&x->reserved2, sizeof(x->reserved2), 1, fp);
	fread(&x->bitmap_dataoffset, sizeof(x->bitmap_dataoffset), 1, fp);
	fread(&x->bitmap_headersize, sizeof(x->bitmap_headersize), 1, fp);
	fread(&x->width, sizeof(x->width), 1, fp);
	fread(&x->height, sizeof(x->height), 1, fp);
	fread(&x->planes, sizeof(x->planes), 1, fp);
	fread(&x->bits_perpixel, sizeof(x->bits_perpixel), 1, fp);
	fread(&x->compression, sizeof(x->compression), 1, fp);
	fread(&x->bitmap_datasize, sizeof(x->bitmap_datasize), 1, fp);
	fread(&x->hresolution, sizeof(x->hresolution), 1, fp);
	fread(&x->vresolution, sizeof(x->vresolution), 1, fp);
	fread(&x->usedcolors, sizeof(x->usedcolors), 1, fp);
	fread(&x->importantcolors, sizeof(x->importantcolors), 1, fp);
}


///////////////////儲存輸入檔案的data/////////////////////
typedef struct Gray{
	int R;
	int G;
	int B;
} ImgGray;

//////////////////////////////////分割線///////////////////////////////////////////////////////
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
///////////////////////讀取RGB資料//////////////////////////////
void InputData2(FILE* fp, ImgGray **array, int H, int W, int skip){
	int temp;
	char skip_buf[3];
	int i, j;
	for (i = 0; i<H; i++){
		for (j = 0; j<W; j++){
			temp = fgetc(fp);
			array[i][j].B = temp;
			temp = fgetc(fp);
			array[i][j].G = temp;
			temp = fgetc(fp);
			array[i][j].R = temp;
		}
		if (skip != 0) fread(skip_buf, skip, 1, fp);
	}
}
////////////////////////////main function//////////////////////////////////
int main(int argc, char **argv){
////////////////////////////////////////////////////////////////////////////////////////////////////////////
if(argc==7){
  FILE *fp_in;//輸出
  char *qqq=argv[1];
	fp_in = fopen(argv[2], "rb");//input輸入 bmp file
	//read input bmp file header
	Bitmap bmpheader;//for input bmp file header
	readheader(fp_in, &bmpheader);
  // W:3024 x H:4032 for houses.bmp
	int H = bmpheader.height;
	int W = bmpheader.width;

  FILE *Rwrite=fopen(argv[3],"w");
  FILE *Gwrite=fopen(argv[4],"w");
  FILE *Bwrite=fopen(argv[5],"w");
        
	int temp;
	int i, j;
	for (i = 0; i<H; i++){
		for (j = 0; j<W; j++){
			temp = fgetc(fp_in);
			fprintf(Bwrite,"%d ",temp);
			
			temp = fgetc(fp_in);
			fprintf(Gwrite,"%d ",temp);
			
			temp = fgetc(fp_in);
			fprintf(Rwrite,"%d ",temp);
		}
	}
	fclose(Rwrite);
	fclose(Gwrite);
	fclose(Bwrite);
        FILE *dim=fopen(argv[6],"w");
        fprintf(dim,"H: %d\n",H);
        fprintf(dim,"W: %d\n",W);
        fclose(dim);
	// 8x8, if not multiples of 8, then bitmap padded, i.e. one more block
	////寬度或高度不是8的倍數，列或行的block個數+1
	int block_H = H / 8;
	if (H % 8 != 0) block_H++;
	int block_W = W / 8;
	if (W % 8 != 0) block_W++;
	//Height and width after blocking
	int bH=8*block_H;//處理過的H
	int bW=8*block_W;//處理過的W
	return 0;
}
  
////////////////////////////方法一/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  else if(argc==13){
  FILE *fp_in;//輸出
  char *qqq1=argv[1];
	fp_in = fopen(argv[2], "rb");//input輸入 bmp file

	//read input bmp file header
	Bitmap bmpheader;//for input bmp file header
	readheader(fp_in, &bmpheader);
  // W:512 x H:512 for houses.bmp
	int H = bmpheader.height;
	int W = bmpheader.width;
	// skip paddings at the end of each image line
	////bmp檔每行會對齊4bytes,所以若raw data不是4的倍數,讀檔時需要跳過skip個數/寫檔的時候需要補上skip個0.
	int skip = (4 - (bmpheader.width * 3) % 4);
	if (skip == 4) skip = 0;
	char skip_buf[3] = { 0, 0, 0 };
	// 8x8, if not multiples of 8, then bitmap padded, i.e. one more block
	////寬度或高度不是8的倍數，列或行的block個數+1
	int block_H = H / 8;
	if (H % 8 != 0) block_H++;
	int block_W = W / 8;
	if (W % 8 != 0) block_W++;
	//Height and width after blocking
	int bH=8*block_H;//處理過的H
	int bW=8*block_W;//處理過的W

	ImgGray **Data_Gray = (ImgGray**)malloc(sizeof(ImgGray*)*H);
	for(int i=0; i<H ;i++)
        Data_Gray[i]=(ImgGray*)malloc(sizeof(ImgGray)*W);
////////////////////////////////////read Gray////////////////////////////////////////////
  InputData2(fp_in, Data_Gray, bmpheader.height, bmpheader.width, skip);//read houses.bmp gray data
  fclose(fp_in);

	YCbCr **f=mal_2D(f,bH,bW);
	for(int i=0;i<H;i++){
		for(int j=0;j<W;j++){
			f[i][j].Y=0.299*(double)Data_Gray[i][j].R+0.587*(double)Data_Gray[i][j].G+0.114*(double)Data_Gray[i][j].B-128;
			f[i][j].Cb=-0.169*(double)Data_Gray[i][j].R-0.331*(double)Data_Gray[i][j].G+0.5*(double)Data_Gray[i][j].B;
			f[i][j].Cr=0.5*(double)Data_Gray[i][j].R-0.419*(double)Data_Gray[i][j].G-0.081*(double)Data_Gray[i][j].B;
			
		}
	}	
	//free
	for(int i=0;i<H;i++){
           free(Data_Gray[i]);
           }
	free(Data_Gray);
///////////////////////////////////output 64 data//////////////////////////////
//generate basis vector for 2D-DCT
//By matlab formula
  double basic[8][8][8][8];
  FILE *outfile;//輸出
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
	
/////////////////////////////////DCT//////////////////////////////////
	YCbCr **F=mal_2D(F,H,W);
	double tempY=0;
	double tempCb=0;
	double tempCr=0;
	double beta[8];
  //根據matlab寫入beta資料
	for(int i=0; i<8; i++){//C[0] = (1/√2) , beat[1]~beat[7]為1
		if(i==0){ beta[i] = 1.0/sqrt(2); }
		else{ beta[i] = 1.0; }
	}

  //根據matlab進行DCT
	for(int i=0;i<H/8;i++){
		for(int j=0;j<W/8;j++){
			for(int u=0;u<8;u++){
				for(int v=0;v<8;v++){
					for(int r=0;r<8;r++){
						for(int c=0;c<8;c++){
							tempY=tempY+f[i*8+r][j*8+c].Y*basic[u][v][r][c];
							tempCb=tempCb+f[i*8+r][j*8+c].Cb*basic[u][v][r][c];
							tempCr=tempCr+f[i*8+r][j*8+c].Cr*basic[u][v][r][c];
						}
					}
					F[i*8+u][j*8+v].Y=beta[u]*beta[v]*0.25*tempY; //0.25=(2/sqrt(8*8))
					tempY=0;
					F[i*8+u][j*8+v].Cb=beta[u]*beta[v]*0.25*tempCb; //0.25=(2/sqrt(8*8))
					tempCb=0;
					F[i*8+u][j*8+v].Cr=beta[u]*beta[v]*0.25*tempCr; //0.25=(2/sqrt(8*8))
					tempCr=0;
				}
			}
		}
	}
	free_2D(f,bH);

///////////////////////////////////////quantize//////////////////////////////////////////////////////	
	//根據matlab進行量化
	YCbCr **quan=mal_2D(quan,H,W);
   
  FILE *QY=fopen(argv[3],"w");
	FILE *QCb=fopen(argv[4],"w");
	FILE *QCr=fopen(argv[5],"w");
	
        for(int u=0;u<8;u++){
				for(int v=0;v<8;v++){
        fprintf(QY,"%f ",Q[u][v]);
        fprintf(QCb,"%f ",Q2[u][v]);
        fprintf(QCr,"%f ",Q2[u][v]);
        if(v==7){
        fprintf(QY,"%f \n",Q[u][v]);
        fprintf(QCb,"%f \n",Q2[u][v]);
        fprintf(QCr,"%f \n",Q2[u][v]);
         }
				}
			}
  FILE *dim=fopen(argv[6],"w");
      fprintf(dim,"H: %d\n",H);
      fprintf(dim,"W: %d\n",W);
  fclose(dim);
	
	FILE *qY=fopen(argv[7],"wb");
	FILE *qCb=fopen(argv[8],"wb");
	FILE *qCr=fopen(argv[9],"wb");
	
	FILE *eY=fopen(argv[10],"wb");
	FILE *eCb=fopen(argv[11],"wb");
	FILE *eCr=fopen(argv[12],"wb");
	

  for(int i=0;i<H/8;i++){
		for(int j=0;j<W/8;j++){
			for(int u=0;u<8;u++){
				for(int v=0;v<8;v++){
					quan[i*8+u][j*8+v].QY=(short)round(F[i*8+u][j*8+v].Y/Q[u][v]);
					quan[i*8+u][j*8+v].eY=F[i*8+u][j*8+v].Y-(quan[i*8+u][j*8+v].QY*Q[u][v]);
          
					fwrite(&quan[i*8+u][j*8+v].QY,sizeof(short),1,qY);
					fwrite(&quan[i*8+u][j*8+v].eY,sizeof(float),1,eY);
					
					quan[i*8+u][j*8+v].QCb=(short)round(F[i*8+u][j*8+v].Cb/Q2[u][v]);
					quan[i*8+u][j*8+v].eCb=F[i*8+u][j*8+v].Cb-(quan[i*8+u][j*8+v].QCb*Q2[u][v]);
          

					fwrite(&quan[i*8+u][j*8+v].QCb,sizeof(short),1,qCb);
					fwrite(&quan[i*8+u][j*8+v].eCb,sizeof(float),1,eCb);
					
					quan[i*8+u][j*8+v].QCr=(short)round(F[i*8+u][j*8+v].Cr/Q2[u][v]);
					quan[i*8+u][j*8+v].eCr=F[i*8+u][j*8+v].Cr-(quan[i*8+u][j*8+v].QCr*Q2[u][v]);
          
					fwrite(&quan[i*8+u][j*8+v].QCr,sizeof(short),1,qCr);
					fwrite(&quan[i*8+u][j*8+v].eCr,sizeof(float),1,eCr);
					
				}
			}
		}
	}
  /////////////計算 3*64 SQNR//////////////////
  double esqnrY[8][8]={0,0};
  double esqnrCb[8][8]={0,0};
  double esqnrCr[8][8]={0,0};

  double p0y[8][8]={0,0};
  double p0cb[8][8]={0,0};
  double p0cr[8][8]={0,0};
  
  //加起來
  for(int i=0;i<8;i++){
    for(int j=0;j<8;j++){
      for(int k=0;k<H/8;k++){
      esqnrY[i][j]=((double)quan[i+8*k][j].eY*(double)quan[i+8*k][j].eY)+esqnrY[i][j];
      p0y[i][j]=(F[i+8*k][j].Y*F[i+8*k][j].Y)+p0y[i][j];

      esqnrCb[i][j]=((double)quan[i+8*k][j].eCb*(double)quan[i+8*k][j].eCb)+esqnrCb[i][j];
      p0cb[i][j]=(F[i+8*k][j].Cb*F[i+8*k][j].Cb)+p0cb[i][j];

      esqnrCr[i][j]=((double)quan[i+8*k][j].eCr*(double)quan[i+8*k][j].eCr)+esqnrCr[i][j];
      p0cr[i][j]=(F[i+8*k][j].Cr*F[i+8*k][j].Cr)+p0cr[i][j];
      }
    }
  }
  //取平均
  for(int i=0;i<8;i++){
    for(int j=0;j<8;j++){
      esqnrY[i][j]=esqnrY[i][j]/504; //pY
      esqnrCb[i][j]=esqnrCb[i][j]/504;
      esqnrCr[i][j]=esqnrCr[i][j]/504;

      p0y[i][j]=p0y[i][j]/504;
      p0cb[i][j]=p0cb[i][j]/504;
      p0cr[i][j]=p0cr[i][j]/504;
    }
  }
  //算SQNR
  double sqnrY[8][8];
  double sqnrCb[8][8];
  double sqnrCr[8][8];
  for(int i=0;i<8;i++){
    for(int j=0;j<8;j++){
      sqnrY[i][j]=10*log10(esqnrY[i][j]/p0y[i][j]);
      sqnrCb[i][j]=10*log10(esqnrCb[i][j]/p0cb[i][j]);
      sqnrCr[i][j]=10*log10(esqnrCr[i][j]/p0cr[i][j]);
    } 
  }
  ///為了格式
  for(int i=0;i<8;i++){
    for(int u=0;u<8;u++){
      if(i==7&&u==7){
        printf("%.1f  \n",sqnrY[i][u]);
        break;
      }
  printf("%.1f  ",sqnrY[i][u]);
    }
  }
  for(int i=0;i<8;i++){
    for(int u=0;u<8;u++){
      if(i==7&&u==7){
        printf("%.1f  \n",sqnrCb[i][u]);
        break;
      }
  printf("%.1f  ",sqnrCb[i][u]);
    }
  }
  for(int i=0;i<8;i++){
    for(int u=0;u<8;u++){
      if(i==7&&u==7){
        printf("%.1f  \n",sqnrCr[i][u]);
        break;
      }
  printf("%.1f  ",sqnrCr[i][u]);
    }
  }
	free_2D(F,bH);
  free_2D(quan,H);

	fclose(QY);
	fclose(QCb);
	fclose(QCr);
	
	fclose(qY);
	fclose(qCb);
	fclose(qCr);

	fclose(eY);
	fclose(eCb);
	fclose(eCr);
}
////////////////////////////////////////////////////////////////////方法2//////////////////////////////////////////////////////////////////////
else if(argc==5){
  char *qqq1=argv[1];
  FILE *fp_in;//輸出
	fp_in = fopen(argv[2], "rb");//input輸入 bmp file
	//read input bmp file header
	Bitmap bmpheader;//for input bmp file header
	readheader(fp_in, &bmpheader);
  // W:512 x H:512 for houses.bmp
	int H = bmpheader.height;
	int W = bmpheader.width;
	// skip paddings at the end of each image line
	////bmp檔每行會對齊4bytes,所以若raw data不是4的倍數,讀檔時需要跳過skip個數/寫檔的時候需要補上skip個0.
	int skip = (4 - (bmpheader.width * 3) % 4);
	if (skip == 4) skip = 0;
	char skip_buf[3] = { 0, 0, 0 };
	// 8x8, if not multiples of 8, then bitmap padded, i.e. one more block
	////寬度或高度不是8的倍數，列或行的block個數+1
	int block_H = H / 8;
	if (H % 8 != 0) block_H++;
	int block_W = W / 8;
	if (W % 8 != 0) block_W++;
	//Height and width after blocking
	int bH=8*block_H;
	int bW=8*block_W;

	ImgGray **Data_Gray = (ImgGray**)malloc(sizeof(ImgGray*)*H);
	for(int i=0; i<H ;i++)
        Data_Gray[i]=(ImgGray*)malloc(sizeof(ImgGray)*W);
////////////////////////////////////read Gray////////////////////////////////////////////
  InputData2(fp_in, Data_Gray, bmpheader.height, bmpheader.width, skip);//read houses.bmp gray data
  fclose(fp_in);

	YCbCr **Data_FGray=mal_2D(Data_FGray,H,W);
	for(int i=0;i<H;i++){
		for(int j=0;j<W;j++){
			Data_FGray[i][j].Y=0.299*(double)Data_Gray[i][j].R+0.587*(double)Data_Gray[i][j].G+0.114*(double)Data_Gray[i][j].B-128;
			Data_FGray[i][j].Cb=-0.169*(double)Data_Gray[i][j].R-0.331*(double)Data_Gray[i][j].G+0.5*(double)Data_Gray[i][j].B;
			Data_FGray[i][j].Cr=0.5*(double)Data_Gray[i][j].R-0.419*(double)Data_Gray[i][j].G-0.081*(double)Data_Gray[i][j].B;
			
		}
	}	
	//free
	for(int i=0;i<H;i++){
           free(Data_Gray[i]);
           }
	free(Data_Gray);
	
	//in a temp
	YCbCr **f=mal_2D(f,bH,bW);
	for(int i=0;i<bH;i++){
		for(int j=0;j<bW;j++){
                        
				f[i][j].Y=Data_FGray[i][j].Y;
				f[i][j].Cb=Data_FGray[i][j].Cb;
				f[i][j].Cr=Data_FGray[i][j].Cr;
				
		}
	}
	free_2D(Data_FGray,H);//free
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
	
/////////////////////////////////DCT//////////////////////////////////
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

  //根據matlab進行DCT
	for(int i=0;i<H/8;i++){
		for(int j=0;j<W/8;j++){
			for(int u=0;u<8;u++){
				for(int v=0;v<8;v++){
					for(int r=0;r<8;r++){
						for(int c=0;c<8;c++){
							tempY=tempY+f[i*8+r][j*8+c].Y*basic[u][v][r][c];
							tempCb=tempCb+f[i*8+r][j*8+c].Cb*basic[u][v][r][c];
							tempCr=tempCr+f[i*8+r][j*8+c].Cr*basic[u][v][r][c];
						}
					}
					F[i*8+u][j*8+v].Y=beta[u]*beta[v]*0.25*tempY; //0.25=(2/sqrt(8*8))
					tempY=0;
					F[i*8+u][j*8+v].Cb=beta[u]*beta[v]*0.25*tempCb; //0.25=(2/sqrt(8*8))
					tempCb=0;
					F[i*8+u][j*8+v].Cr=beta[u]*beta[v]*0.25*tempCr; //0.25=(2/sqrt(8*8))
					tempCr=0;
				}
			}
		}
	}
	free_2D(f,bH);
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

///////////////////////////////////////quantize//////////////////////////////////////////////////////	
	//根據matlab進行量化
  for(int i=0;i<H/8;i++){
		for(int j=0;j<W/8;j++){
			for(int u=0;u<8;u++){
				for(int v=0;v<8;v++){
					quan[i*8+u][j*8+v].Y=round(F[i*8+u][j*8+v].Y/Q_Y[u][v]);
					quan[i*8+u][j*8+v].Cb=round(F[i*8+u][j*8+v].Cb/Q_CbCr[u][v]);
					quan[i*8+u][j*8+v].Cr=round(F[i*8+u][j*8+v].Cr/Q_CbCr[u][v]);
				}
			}
		}
	}
	free_2D(F,bH);
////////////////////////////////////////////DPCM///////////////////////////////////////////////////////////
	YCbCr **DPCM=mal_2D(DPCM,bH,bW);

	for(int i=bH-1;i>=0;i--){
		for(int j=bW-1;j>=0;j--){
			if(i>0 && j==0){
				quan[i][j].Y= quan[i][j].Y-quan[i-1][j].Y;
				quan[i][j].Cb = quan[i][j].Cb-quan[i-1][j].Cb;
				quan[i][j].Cr = quan[i][j].Cr-quan[i-1][j].Cr;
			}
			else if(j>0){
				quan[i][j].Y=quan[i][j].Y - quan[i][j-1].Y;
				quan[i][j].Cb= quan[i][j].Cb -quan[i][j-1].Cb;
				quan[i][j].Cr=quan[i][j].Cr - quan[i][j-1].Cr;
			}
		}
	}
	for(int i=0; i<bH; i++){
		for(int j=0; j<bW; j++){
      DPCM[i][j].Y=quan[i][j].Y;
			DPCM[i][j].Cb=quan[i][j].Cb;
			DPCM[i][j].Cr=quan[i][j].Cr;

		}
	}
	free_2D(quan,bH);
//////////////////////////////////////zigzag table/////////////////////////////
	int vec[8][8]={{ 0, 1, 5, 6,14,15,27,28},
				   { 2, 4, 7,13,16,26,29,42},
 				   { 3, 8,12,17,25,30,41,43},
 				   { 9,11,18,24,31,40,44,53},
 				   {10,19,23,32,39,45,52,54},
 				   {20,22,33,38,46,51,55,60},
 				   {21,34,37,47,50,56,59,61},
 				   {35,36,48,49,57,58,62,63}};

///////////////////////////////////////zigzag//////////////////////////////////////
	int count=0;
	YCbCr *zig=mal_1D(zig,bH*bW);

	for(int i=0;i<H/8;i++){
		for(int j=0;j<W/8;j++){
			for(int u=0;u<8;u++){
				for(int v=0;v<8;v++){
					//vec[][]=0~63
					//0~63 64~127.....
					zig[8*8*count+vec[u][v]].Y=DPCM[i*8+u][j*8+v].Y;
					zig[8*8*count+vec[u][v]].Cb=DPCM[i*8+u][j*8+v].Cb;
					zig[8*8*count+vec[u][v]].Cr=DPCM[i*8+u][j*8+v].Cr;
				}	
			}
			count++;
		}
	}
	free_2D(DPCM,bH);
//////////////////////////////////////RLE///////////////////////////////////////
	YCbCr *RLE=mal_1D(RLE,2*bH*bW);
	//依序讀經過zigzag過後的數據
  int total=0; 
  int nzeroy=0;
  int nzerocb=0;
  int nzerocr=0;
  int index=0;
  int check=0;
  char *qqq2=argv[3];//binary or ascii
  
  //判斷是ascii or binary
  int checknum=2;
  if(strcmp(argv[3],"binary")==0) checknum=0;
  else if (strcmp(argv[3],"ascii")==0) checknum=1;
  
  printf("%d",checknum);
  
//////////////////////////////////////////////txt version/////////////////////////////////////
  if(checknum==1){
  FILE *bin=fopen(argv[4], "wb");
  int average=(H*W)/64;//for more clear
	while(total<average){
		for(int buff=0;buff<64;buff++){
      //最後一個就為00
      //再來就是分兩種情況[0,數字][數字,數字]
      //第一格為0連續出現的數量
      //第二個為打斷0連續出現的數=>即可存入
      if(buff==63){
				RLE[index].rY=0;
				RLE[index+1].rY=0;
				break;

			}
			if(buff!=63&&zig[total*64+buff].Y==0) {
			nzeroy++;
			}
			else if(buff!=63&&zig[total*64+buff].Y!=0){
				RLE[index].rY=nzeroy;
        RLE[index].rY=zig[total*64+buff].Y;
        index=index+2;
				nzeroy=0;
			}
		}
    total=total+1;
  }
     
	index=0;
        total=0; 
	while(total<H*W){
		for(int buff=0;buff<64;buff++){
			if(buff==63){
				RLE[index].rCb=0;
				RLE[index+1].rCb=0;
        break;
			}			
			if(zig[total+buff].Cb==0){
			nzerocb++;
			}
			else if(zig[total+buff].Cb!=0){
				RLE[index].rCb=nzerocb;
        RLE[index+1].rCb=zig[total+buff].Cb;
				nzerocb=0;
				index+=2;
			}
		}
    total=total+64;
  }

	index=0;
        total=0;
	while(total<H*W){
		for(int buff=0;buff<64;buff++){
			if(buff==63){
				RLE[index].rCr=0;
				RLE[index+1].rCr=0;
				break;
			}
			if(zig[total+buff].Cr==0){
			 nzerocr++;
			 }
			else if(zig[total+buff].Cr!=0){
				RLE[index].rCr=nzerocr;
        RLE[index+1].rCr=zig[total+buff].Cr;
				nzerocr=0;
				index+=2;
			}
		}
    total=total+64;
  }
  
  free(zig);
  ////為了符合格式的寫法
  FILE *ascii=fopen(argv[4], "w");
  for(int x=0;x<H/8;x++){
     for(int y=0;y<W/8;y++){
       ////////////Y//////////
       fprintf(ascii," %d %d Y",x,y);
       for(int indexY=0;indexY<64;indexY+=2){
       if(RLE[indexY].rY==0&&RLE[indexY+1].rY==0){
       fprintf(ascii," %d %d\n",RLE[indexY].rY,RLE[indexY+1].rY);
       break;
       }
       else{
       fprintf(ascii," %d %d",RLE[indexY].rY,RLE[indexY+1].rY);
       }
       }
     /////////////////Cb//////////////
       fprintf(ascii," %d %d Cb",x,y);
       for(int indexCb=0;indexCb<64;indexCb+=2){
       if(RLE[indexCb].rCb==0&&RLE[indexCb+1].rCb==0){
       fprintf(ascii," %d %d\n",RLE[indexCb].rCb,RLE[indexCb+1].rCb);
       break;
       }
       else{
       fprintf(ascii," %d %d",RLE[indexCb].rCb,RLE[indexCb+1].rCb);
       
       }
       }
       //////////////////////Cr/////////////////
       fprintf(ascii," %d %d Cr",x,y);
       for(int indexCr=0;indexCr<64;indexCr+=2){
       if(RLE[indexCr].rCr==0&&RLE[indexCr+1].rCr==0){
       fprintf(ascii," %d %d\n",RLE[indexCr].rCr,RLE[indexCr+1].rCr);
       break;
       }
       else{
       fprintf(ascii," %d %d",RLE[indexCr].rCr,RLE[indexCr+1].rCr);
       }
       }
     
     }
  }
  
  fclose(ascii);
  free(RLE);
  }

////////////////////////////////////////ascii////////////////////////////////////////
  else if(checknum==0){
  printf("開始寫入bin\n");
  FILE *bin=fopen(argv[4], "wb");
  int average=(H*W)/64;

	while(total<average){
		for(int buff=0;buff<64;buff++){
       if(buff==63){
				RLE[index].rY=0;
				RLE[index+1].rY=0;
				fwrite(&RLE[index].rY,sizeof(int),1,bin);
        fwrite(&RLE[index+1].rY,sizeof(int),1,bin);
        break;
			}           
			if(buff!=63&&zig[total*64+buff].Y==0) {
			nzeroy++;
			}
			else if(buff!=63&&zig[total*64+buff].Y!=0){
				RLE[index].rY=nzeroy;
        RLE[index].rY=zig[total*64+buff].Y;
        fwrite(&RLE[index].rY,sizeof(int),1,bin);
        fwrite(&RLE[index+1].rY,sizeof(int),1,bin);
        index=index+2;
				nzeroy=0;
			}
			
		}
    total=total+1;
  }
     
	index=0;
        total=0; 
	while(total<H*W){
		for(int buff=0;buff<64;buff++){

			if(buff==63){ 
				RLE[index].rCb=0;
				RLE[index+1].rCb=0;
				fwrite(&RLE[index].rCb,sizeof(int),1,bin);
        fwrite(&RLE[index+1].rCb,sizeof(int),1,bin);
        break;
			}		
			if(zig[total+buff].Cb==0){
			nzerocb++;
			}
			if(zig[total+buff].Cb!=0){
    
				RLE[index].rCb=nzerocb;
        RLE[index+1].rCb=zig[total+buff].Cb;
        fwrite(&RLE[index].rCb,sizeof(int),1,bin);
        fwrite(&RLE[index+1].rCb,sizeof(int),1,bin);
        index=index+2;
				nzerocb=0;
			}

		}
    total=total+64;
  }

	index=0;
  total=0;
	while(total<H*W){
		for(int buff=0;buff<64;buff++){

      if(buff==63){ 
				RLE[index].rCr=0;
				RLE[index+1].rCr=0;
				fwrite(&RLE[index].rCb,sizeof(int),1,bin);
        fwrite(&RLE[index+1].rCb,sizeof(int),1,bin);
        break;
			}

			if(zig[total+buff].Cr==0){
			 nzerocr++;
			 }
			if(zig[total+buff].Cr!=0){
				RLE[index].rCr=nzerocr;
        RLE[index+1].rCr=zig[total+buff].Cr;
        fwrite(&RLE[index].rCb,sizeof(int),1,bin);
        fwrite(&RLE[index+1].rCb,sizeof(int),1,bin);
				nzerocr=0;
				index+=2;
			}
	
		}
    total=total+64;
  }
  
  free(zig);
  fclose(bin);
  printf("使用2b的壓縮率為17.09/34.88*100=百分之49\n");
  free(RLE);
}


}

        }//底

