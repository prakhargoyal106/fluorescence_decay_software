/*  EXP2.C
	Two exponential fit
	11/11/92;
	Simulate() added on 14/3/93.
	sd option added on 13/6/94.
	Batch mode added in Sep 98.
	LAB option added -- Oct 2000
	Invert data added -- Oct 2000
	25/2/03
	N.Periasamy */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
///#include <conio.h>
#include <math.h>
#include <malloc.h>

#define LAB 0U	    /*0 for TIFR;1 for Mely*/

#define FILES 9U	    /* Max number of exem files plus 1*/
#define PARMS 27U            /* Max number of parameters plus 1*/
char txt[]="2-exp fit";
#define square(x) ((x)*(x))
#define SWAP(a,b) {float temp=(a);(a)=(b);(b)=temp;}
int n=1;

void plot(char FileName[],int n){
  

  FILE *gnuplotPipe   ;
  
  gnuplotPipe = popen("gnuplot -persistent","w");
  
  fprintf(gnuplotPipe,"set multiplot \n");
  fprintf(gnuplotPipe, "set size 1.00,0.70 \n");
  fprintf(gnuplotPipe, "set origin 0.00,0.30\n");
  fprintf(gnuplotPipe, "set logscale y 10 \n");
  fprintf(gnuplotPipe, "set ylabel \"log counts \" \n");
  fprintf(gnuplotPipe, "set xlabel \"channel no.\" \n");
  
  fprintf(gnuplotPipe, "plot \"%s\" every:\"%d\" u 1:2 w lines  title \"ex vs channel\",\"%s\" every:\"%d\"  u 1:3 w  lines lt -1 title \"em vs channel\",\"%s\" every:\"%d\" u 1:4 w lines lt 1 title \"cal em vs channel\" \n ",FileName,n,FileName,n,FileName,n);
  
  fprintf(gnuplotPipe,"set size 1.00,0.30 \n");
  fprintf(gnuplotPipe,"set origin 0.00,0.00 \n");
  fprintf(gnuplotPipe,"unset logscale y \n");
  fprintf(gnuplotPipe,"set ylabel \"std deviation\" \n");
  
  fprintf(gnuplotPipe,"plot \"%s\" every:\"%d\" u 1:5 w lines title \"residual \" \n",FileName,n);
  
  fflush(gnuplotPipe);

}

void read_parameters(int *jx,int *npar,int *ma,int *nchan,int *iplt, int *cycle,float *xi,\
	float *pa,float *win,int *nb,int *ne,int *nx,int *ip, int *ib,float *bkgd,float *bkgl, char *fname);
void read_exem_file(int jx,char *fname,int nchan, float **il, float **ie,\
     char *zim);
void read_exem_melylab_file(int jx,char *fname,int nchan, float **il, float **ie, char *zim);
void read_ex_em_file(int jx,char *fname1,char *fname,int nchan,\
    float **il, float **ie, char *zim);
void covsrt(float **covar,int ma,int lista[],int mfit);
void gaussj(float **a,int n,float **b,int m);
void core_job(float *yex,float *yval,float *ycal,int jx,int ndata,int npar,int ma, int mfit, int *nb, int *ne, int *nx,int *lista,float *a,float *xi, float **alpha,float *beta,float *chisq,float *sigma);
void mrqcof(float *yex,float *yem,float *ycal,int jx,int ndata,int npar,int ma,	int mfit, int *nb, int *ne, int *nx,int *lista,float *a,float *xi, float **alpha,float *beta,float *chisq,float *sigma);
void mrqmin(float *yex,float *yem,float *ycal,int jx,int ndata,int npar,int ma,	int mfit, int *nb, int *ne, int *nx,int *lista,float *a,float *xi, float **covar,float **alpha,float *chisq,float *alamda,float *sigma);
void convolute(float *yex,int jx,int ndata,int npar,int ma,int mfit,\
	int *nb, int *ne, int *nx,int *lista,float *a,float *xi,float *ycal,\
	float **dyda);
void time_shift(float d,float *yll,float *yls,int nch);
void para_convert( int jx,int ma,float **b, float *a);
int find_max_chan(int jj,int nchan,float **ie);
void plot_res(int jx,int ntot, int *nb,int *ne,int *nx,float *yex, float\
	*yem,float *ycal,float *res, float chisq,int ma,float *a,\
	char *fname,char *txt,float dt);
void randomize_amps(int jx,int ma, float *a, float *pa);
void out_screen(int jx,int ma, float *a,float **alpha);
void out_file(FILE *fres,FILE *fres2,int jx,int cycle,int iter,float chisq,int ma,//
	float *a,float **alpha,float *yex,float *yem, float *ycal,float *res,float *acorr,//
	int ntot,int mfit, int *lista,int c_c);
void find_lifetime(int jx,int nchan,int *nb,int *ne,float *xi,\
    float **ie,float *pa);
void scale_amps(int jx,int nchan,int *nb,int *ne,int *nx,float *xi,float **il,//
	float **ie,float *pa);
void scale_amps_batch(int jx,int nchan,int *nb,int *ne,int *nx,float *xi,float **il,//
	float **ie,float *pa);
void display_exem_file(int jx,int nchan,float **il,float **ie,float *bkgl,//
	float *bkgd,int *nb,int *ne,int *nx);
void simulate(float *yex,float *ycal,int jx,int ntot,int npar,int mfit,//
	int *nb,int *ne,int *nx,int *lista,float *a,float *xi);
void invert_data(int jx,int nchan,float **il,float **ie);



float *vector(int start,int end){
	float *vec=malloc(sizeof(float)*(end-start+3));				///previous use of vector() and following functions not supported by linux 
	return vec;
} 
float **matrix(int a,int row,int b,int col){
	float**mat=malloc(sizeof(float*)*(row-a+3));
	int i;
	for(i=0;i<row-a+3;i++)
		mat[i]=malloc(sizeof(float)*(col-b+3));
	return mat;
}
int *ivector(int start,int end){
	int *vec=malloc(sizeof(int)*(end-start+3));
	return vec;
} 
void free_vector(float *vec,int start, int end){
	free(vec);
}
void free_ivector(int *vec,int start, int end){
	free(vec);
}
void free_matrix(float** mat,int a,int row, int b, int col){
/*	int i;
	for(i=0;i<row-a+3;i++)
		free(mat[i]);
	free(mat);*/
}

struct range{
	float low[PARMS];
	float high[PARMS];
	} par;
float pa[PARMS];

int main(int argc, char *argv[])
{
    int n =1;
	FILE *fd, *fres,*fres2;
	int i,j,k,ij,jj,kk;
	float chisq,alamda;
	int ntot,ma,*lista,mfit;
	float con_factor,chisq_old;
	int iter;
	int nch;
    int nb[FILES],ne[FILES],nx[FILES],ib[PARMS];
	float xi[FILES],a[PARMS],bkgl[FILES],bkgd[FILES],win[PARMS];
	char fname[10];
	float **covar,**alpha;
	int jx,npar,nchan,iplt,ip,cycle;
	float **ie,**il,*yex,*yem,*ycal,*res,*acorr,*sigma,*temp;
	float x;
	int file_open=0;
	int nmax,nch1;
	char ans;
	float chi_up;
	char zim[6];
	char sdfname[10];
	float dt=1.;
	char fname1[10]="ex";
	int jbat=0,jjbat=10;
	int lab=LAB;

	if(lab==0) printf("---TIFR DATA----\n");
	if(lab==1) printf("---MELY_LAB DATA---\n");
	//strcpy(zim,argv[1]);
	printf("Available choices:\n\t Press Any key: Analysis of Data;\
	\n\t batch: Batch files;\
	\n\t sim : simulate;\
	\n\t sd: s.d. data read from a file.\n\n\tYour choice:");
	gets(zim);
	puts(zim);
    if((strcmp(zim,"batch"))==0)
    {
	jx=1;npar=5;ma=5;nchan=512;iplt=0;cycle=5;
	pa[1]=1.;pa[2]=1.;pa[3]=0.5;pa[4]=0.5;pa[5]=0.;
	win[1]=0.;win[2]=0;win[3]=0.;win[4]=0.;win[5]=1.;
	nb[1]=1;ne[1]=512;nx[1]=-5;ip=0;bkgl[1]=0;bkgd[1]=0;
	mfit=ma;ntot=nchan;
	il = matrix(1,jx,1,nchan);
	ie = matrix(1,jx,1,nchan);
	covar=matrix(1,ma,1,ma);
	alpha=matrix(1,ma,1,ma);
	lista = ivector(1,ma);
	yex = vector(1,ntot);
	yem = vector(1,ntot);
	ycal = vector(1,ntot);
	res = vector(1,ntot);
	acorr = vector(1,ntot);
	sigma = vector(1,ntot);
	temp = vector(1,nchan);
	printf("\nEnter the name of Ex file:");
	scanf("%s",fname1);
	getchar();
	printf("\n Give time per channel [ps]:");
	scanf("%f",&xi[1]);
	getchar();
	printf("\n Time per channel = %f ps",xi[1]);
	xi[1]/=1000.;
    }
    else
    {
		read_parameters(&jx,&npar,&ma,&nchan,&iplt,&cycle,xi,pa,win,nb,//
			ne,nx,&ip,ib,bkgd,bkgl,fname);
		if (ma>PARMS){printf("  %d: Too many parameters...\n",ma);exit(0);}
    }
		if(file_open!=1)
		{
			fres=fopen("result","w");
			fres2=fopen("result2","w");
			fprintf(fres,"\n\n*****  T.I.F.R.  CHEMICAL DYNAMICS  ********\n");
			fprintf(fres2,"\n\n*****  T.I.F.R.  CHEMICAL DYNAMICS  ********\n");
			fprintf(fres2,"Two exp fit of File name:%s\n\n",fname);
			fprintf(fres,"%d %d %d %d %d %d\n",jx,npar,ma,nchan,iplt,cycle);
			for(i=1;i<=jx;i++) fprintf(fres,"%2.3f\n",xi[i]);
			for(i=1;i<=jx;i++) fprintf(fres,"%d %d %d\n",nb[i],ne[i],nx[i]);
			for(i=1;i<=ma;i++) fprintf(fres,"%2.3f ",pa[i]);
			fprintf(fres,"\n");
			for(i=1;i<=ma;i++)
			{
				if(win[i]==0.) fprintf(fres,"nw ");
				else fprintf(fres,"%2.3f ",win[i]);
			}
			fprintf(fres,"\n");
			fprintf(fres,"%d\n",ip);
			if (ip>0)
				for(i=1;i<=ip;i++) fprintf(fres,"%d ",ib[i]);
			for(i=1;i<=jx;i++) fprintf(fres,"%2.3f %2.3f \n",bkgl[i],bkgd[i]);
			fprintf(fres,"%s\n\n",fname);
			file_open=1;
		}
  do
  {
    if((strcmp(zim,"batch"))==0)
    {
	jbat++;
	cycle=5;
	pa[1]=1.;pa[2]=1.;pa[3]=0.5;pa[4]=0.5;pa[5]=0.;
	win[1]=0.;win[2]=0;win[3]=0.;win[4]=0.;win[5]=1.;
	nb[1]=1;ne[1]=512;nx[1]=-5;ip=0;bkgl[1]=0;bkgd[1]=0;
	do
	{
	    printf("\n-----\nEnter the name of Em file:");
	    scanf("%s",fname);
		getchar();
	    ans='y';
	    if((fd=fopen(fname,"rt")) == NULL)
	    {
		printf(" Emission file '%s' does not exist!!\n",fname);
		ans='n';
	    }
	}while(ans=='n');
	read_ex_em_file(jx,fname1,fname,nchan, il, ie,zim);
	invert_data(jx,nchan,il,ie);
///	display_exem_file(jx,nchan,il,ie,bkgl,bkgd,nb,ne,nx);
	find_lifetime(jx,nchan,nb,ne,xi,ie,pa);
	scale_amps_batch(jx,nchan,nb,ne,nx,xi,il,ie,pa); /* pa[] & nb,ne,nx modified; */
	for(k=1,i=1;i<=jx;i++)
		for(j=nb[i];j<=ne[i];j++)
		{
			yex[k] = il[i][j] - bkgl[i];
			yem[k] = ie[i][j] - bkgd[i];
			ycal[k] = 0.0;
			k++;
		}
	    ntot=k-1;  /* batch mode */
	//if(k != (ntot+1)){ printf(" Something wrong. Check!"); getchar();}
	for(i=1;i<=ntot;i++){
		if(yex[i]<0.) yex[i]=0.;
		if(yem[i]<0.) yem[i]=0.;
	}
    }
    else
    {
		il = matrix(1,jx,1,nchan);
		ie = matrix(1,jx,1,nchan);
		if(lab==0)read_exem_file(jx,fname,nchan, il, ie,zim);
		if(lab==1)read_exem_melylab_file(jx,fname,nchan, il, ie,zim);
		invert_data(jx,nchan,il,ie);
		//display_exem_file(jx,nchan,il,ie,bkgl,bkgd,nb,ne,nx);
		for(ntot=0,i=1;i<=jx;i++) ntot += ne[i]-nb[i]+1;
		mfit=ma-ip;
		covar=matrix(1,ma,1,ma);
		alpha=matrix(1,ma,1,ma);
		lista = ivector(1,ma);
		yex = vector(1,ntot);
		yem = vector(1,ntot);
		ycal = vector(1,ntot);
		res = vector(1,ntot);
		acorr = vector(1,ntot);
		sigma = vector(1,ntot);
		temp = vector(1,nchan);
		//find_lifetime(jx,nchan,nb,ne,xi,ie,pa);
		scale_amps(jx,nchan,nb,ne,nx,xi,il,ie,pa); /* pa[] & nx[] modified; */
		for(k=1,i=1;i<=jx;i++)
			for(j=nb[i];j<=ne[i];j++)
			{
				yex[k] = il[i][j] - bkgl[i];
				yem[k] = ie[i][j] - bkgd[i];
				ycal[k] = 0.0;
				k++;
			}
		if(k != (ntot+1)){ printf(" Something wrong. Check!"); getchar();}
		for(i=1;i<=ntot;i++){
			if(yex[i]<0.) yex[i]=0.;
			if(yem[i]<0.) yem[i]=0.;
		}
    }
	for(i=1;i<=ma;i++) a[i]=pa[i];
	for(i=1;i<=ma;i++)
	{
		if(win[i]==0.)
		{
			if(i<=2){
			par.low[i]=0.0;
			par.high[i]=100.;
			}
			else{
			par.low[i]=-10.;
			par.high[i]=10.;
			}
		}
		else
		{
			par.low[i]=pa[i]-win[i];
			par.high[i]=pa[i]+win[i];
		}
	}
		//getchar();
	for (i=1;i<=ip;i++) lista[ma-i+1] = ib[i];
	ij = 0;
	for(i=1;i<=ma;i++)
	{
		k=0;
		for (j=1;j<=ip;j++)
			if(i == ib[j]) k=1;
		if (k==0)
			lista[++ij] = i;
	}
	for(i=1;i<=ma;i++) printf("%d ",lista[i]);
	printf("\n pa[1]=%4.1f, pa[2]=%4.1f pa[3]=%4.1f pa[4]=%4.1f\n",\
	    pa[1],pa[2],pa[3],pa[4]);
	printf("\tOKAY?\n");

	if((strcmp(zim,"sim"))==0)
	{
		simulate(yex,ycal,jx,ntot,npar,mfit,nb,ne,nx,lista,a,xi);
		exit(1);
	}
	if((strcmp(zim,"sd"))==0)
	{
		printf("\nGive filename for s.d. data: ");
		scanf("%s",sdfname);
		getchar();
		printf("\n%s",sdfname);
		if((fd=fopen(sdfname,"rt"))!=NULL)
		{
			if(jx>1){
				printf("\n Only one exem file allowed for external s.d. file.");
				exit(0);
			}
			for(i=1;i<=64;i++) for(j=1;j<=8;j++)
				fscanf(fd,"%f",&temp[(i-1)*8+j]);
		}
		else{
			printf("\n No s.d. file.");
			exit(0);
		}
		for(i=1;i<=ntot;i++)
			sigma[i]=(temp[i+nb[1]-1]>1.) ? temp[i+nb[1]-1] : 1.;
		for(i=1;i<=ntot;i++){
			printf("\n%d %f %f",i,yem[i],sigma[i]);
			//getchar();
		}
		fclose(fd);
	}
	else
	{
		for(j=1;j<=ntot;j++)
			sigma[j]=(yem[j]>1.) ? (float) sqrt((double)yem[j]) : 1.;
	}

	printf("\n\n Give value of Upper chisq:");
	scanf("%f",&chi_up);
	getchar();
	if(chi_up<=0.5)
	{
		printf(" \n\n You are unreasonable with your chisq requirement!\
		Not likely!");
		getchar();
	}

	j=1;
	do{
		iter=0;
		alamda=-1.;
		mrqmin(yex,yem,ycal,jx,ntot,npar,ma,mfit,nb,ne,nx,lista,a,\
			xi,covar,alpha,&chisq,&alamda,sigma);
		chisq /= (ntot-mfit+1);
		//printf("Iter=%d, chisq=%f alamda=%f\n",iter,chisq,alamda);
repeat:
		//printf("\ntau1=%2.3f,tau2=%2.3f,amp1=%2.3f,amp2=%2.3f\n",a[1],a[2],a[3],a[4]);
		//printf("Iter=%d, chisq=%f alamda=%f\n",iter,chisq,alamda);
		printf(" Working...Wait! chisq=%f\r ",chisq);
		chisq_old=chisq;
		iter++;
		mrqmin(yex,yem,ycal,jx,ntot,npar,ma,mfit,nb,ne,nx,lista,a,\
			xi,covar,alpha,&chisq,&alamda,sigma);
		chisq /= (ntot-mfit+1);
		//getchar();
		if(chisq==chisq_old && alamda<10.) goto repeat;
		con_factor=(chisq_old-chisq)/chisq_old;
		if(con_factor>0.01&&iter<100&&alamda<.1) goto repeat;
		printf("Cycle %d Over...\n",j);
		printf("Iter=%d, chisq=%f alamda=%f\n",iter,chisq,alamda);
		if(chisq<chi_up)
		{
			ij = 0;
			alamda=0.;
			mrqmin(yex,yem,ycal,jx,ntot,npar,ma,mfit,nb,ne,nx,lista,a,\
				xi,covar,alpha,&chisq,&alamda,sigma);
			chisq /= (ntot-mfit+1);
			for(jj=1;jj<=jx;jj++)
			{
				nch1=ne[jj]-nb[jj]+1;
				for(i=ij+1;i<=ij+nch1;i++)
				{
					res[i] = 0.;
					if(i>=ij+nx[jj])
					{
						alamda = (yem[i]>1.0) ? (float)sqrt((double)yem[i]) : 1.0;
						res[i] = (ycal[i]-yem[i])/alamda;
					}
				}
				ij += nch1;
			}
			ij = nx[1]-1;
			for(i=1;i<=ntot;i++) acorr[i]=0.;
			for(jj=1;jj<=jx;jj++){
				nch1=ne[jj]-nb[jj]+1;
				kk = (nch1-nx[jj])/2;
				for(k=1;k<=kk;k++){
					for(i=ij;i<=ij+kk;i++) acorr[k+ij] += \
						res[i]*res[i+k-1];
					acorr[k+ij] /= (kk*chisq);
				}
				if(jj<jx) ij += nch1-nx[jj]+nx[jj+1];
			}
			out_file(fres,fres2,jx,cycle,iter,chisq,ma,a,alpha,yex,yem,ycal,res,acorr,ntot,mfit,lista,j);
			out_screen(jx,ma,a,alpha);
			//fflush(stdin);
			//plotfunction();

	        //plotfunction(j);
			getchar();
			
			dt=xi[1];
	
		}
		randomize_amps(jx,ma,a,pa);
		j++;
		//plotfunction();
		//getchar();
		if(j==cycle+1)
		{
			printf("/nDo You want Another 5 Analyses?[y/n]:");
//			fflush(stdin);
			ans=getchar();
			getchar();
			if(ans=='Y' || ans=='y') cycle += 5;
		    if(cycle>=5 && (strcmp(zim,"batch"))==0) nx[1]+=5;
		}
	  
     //plotfunction();
	  }while(j<=cycle);
  }while(jbat>0 && jbat<jjbat);
	fclose(fres);
	fclose(fres2);
	free_matrix(covar,1,ma,1,ma);
	free_matrix(alpha,1,ma,1,ma);
	free_matrix(il,1,jx,1,nchan);
	free_matrix(ie,1,jx,1,nchan);
	free_ivector(lista,1,ma);
	free_vector(yex,1,ntot);
	free_vector(yem,1,ntot);
	free_vector(ycal,1,ntot);
	free_vector(res,1,ntot);
	free_vector(acorr,1,ntot);
	free_vector(sigma,1,ntot);
	free_vector(temp,1,nchan);

	//plotfunction();

	printf("DONE\n");
	return(1);
}


/*void plotfunction(int c){

	int n = c;
	FILE *gnuplotPipe , *tempDataFile  ;
 //printf("Enter the value of n:\n");
  //scanf("%d",&n);
  
 // printf("Enter the file name : ");

  char FileName[7] = ;
 //scanf("%s",FileName);
  //FileName ;
  //tempDataFile = fopen(FileName,"r");

  
  gnuplotPipe = popen("gnuplot -persistent","w");
  
  fprintf(gnuplotPipe,"set multiplot \n");
  fprintf(gnuplotPipe, "set size 1.00,0.70 \n");
  fprintf(gnuplotPipe, "set origin 0.00,0.30\n");
  fprintf(gnuplotPipe, "set logscale y 10 \n");
  fprintf(gnuplotPipe, "set ylabel \"log counts \" \n");
  fprintf(gnuplotPipe, "set xlabel \"channel no.\" \n");
  //  fprintf(gnuplotPipe, "plot \"%s\" every:\"%s\" u 1:2 w lines title \"ex vs channel\",\"%s\" every:\"%s\"  u 1:3 w lines title \"em vs channel\",\"%s\" every:\"%s\" u 1:4 w lines title \"cal em vs channel\" \n ",FileName,n,FileName,n,FileName,n);
  if(n==1) fprintf(gnuplotPipe, "plot \"%s\" every:1 u 1:2 w lines title \"ex vs channel\",\"%s\" every:1  u 1:3 w lines title \"em vs channel\",\"%s\" every:1 u 1:4 w lines title \"cal em vs channel\" \n ",FileName,FileName,FileName);
  if(n==2) fprintf(gnuplotPipe, "plot \"%s\" every:2 u 1:2 w lines title \"ex vs channel\",\"%s\" every:2  u 1:3 w lines title \"em vs channel\",\"%s\" every:2 u 1:4 w lines title \"cal em vs channel\" \n ",FileName,FileName,FileName);
  if(n==3) fprintf(gnuplotPipe, "plot \"%s\" every:3 u 1:2 w lines title \"ex vs channel\",\"%s\" every:3  u 1:3 w lines title \"em vs channel\",\"%s\" every:3 u 1:4 w lines title \"cal em vs channel\" \n ",FileName,FileName,FileName);
  if(n==4) fprintf(gnuplotPipe, "plot \"%s\" every:4 u 1:2 w lines title \"ex vs channel\",\"%s\" every:4  u 1:3 w lines title \"em vs channel\",\"%s\" every:4 u 1:4 w lines title \"cal em vs channel\" \n ",FileName,FileName,FileName);
  if(n==5) fprintf(gnuplotPipe, "plot \"%s\" every:5 u 1:2 w lines title \"ex vs channel\",\"%s\" every:5  u 1:3 w lines title \"em vs channel\",\"%s\" every:5 u 1:4 w lines title \"cal em vs channel\" \n ",FileName,FileName,FileName);
  if(n==6) fprintf(gnuplotPipe, "plot \"%s\" every:6 u 1:2 w lines title \"ex vs channel\",\"%s\" every:6  u 1:3 w lines title \"em vs channel\",\"%s\" every:6 u 1:4 w lines title \"cal em vs channel\" \n ",FileName,FileName,FileName);
  if(n==7) fprintf(gnuplotPipe, "plot \"%s\" every:7 u 1:2 w lines title \"ex vs channel\",\"%s\" every:7  u 1:3 w lines title \"em vs channel\",\"%s\" every:7 u 1:4 w lines title \"cal em vs channel\" \n ",FileName,FileName,FileName);
  if(n==8) fprintf(gnuplotPipe, "plot \"%s\" every:8 u 1:2 w lines title \"ex vs channel\",\"%s\" every:8  u 1:3 w lines title \"em vs channel\",\"%s\" every:8 u 1:4 w lines title \"cal em vs channel\" \n ",FileName,FileName,FileName);
  if(n==9) fprintf(gnuplotPipe, "plot \"%s\" every:9 u 1:2 w lines title \"ex vs channel\",\"%s\" every:9  u 1:3 w lines title \"em vs channel\",\"%s\" every:9 u 1:4 w lines title \"cal em vs channel\" \n ",FileName,FileName,FileName);
  if(n==10) fprintf(gnuplotPipe, "plot \"%s\" every:10 u 1:2 w lines title \"ex vs channel\",\"%s\" every:10  u 1:3 w lines title \"em vs channel\",\"%s\" every:10 u 1:4 w lines title \"cal em vs channel\" \n ",FileName,FileName,FileName);
  //fprintf(gnuplotPipe, "plot \"%s\" u 1:2 w lines title \"ex vs channel\",\"%s\" u 1:3 w lines title \"em vs channel\",\"%s\" u 1:4 w lines title \"cal em vs channel\" \n ",FileName,FileName,FileName);
 // fflush(gnuplotPipe);
  fprintf(gnuplotPipe,"set size 1.00,0.30 \n");
  fprintf(gnuplotPipe,"set origin 0.00,0.00 \n");
  fprintf(gnuplotPipe,"unset logscale y \n");
  fprintf(gnuplotPipe,"set ylabel \"std deviation\" \n");
  fprintf(gnuplotPipe,"plot \"%s\" u 1:5 w lines title \"residual \" \n",FileName);
  //char();
  //fprintf(gnuplotPipe,"exit");
  //fclose(gnuplotPipe);
  fflush(gnuplotPipe);
  //fprintf(gnuplotPipe, "replot\n");
		//printf("\n----PLOT Option not yet ready-------\n"
}*/

void read_parameters(int *jx,int *npar,int *ma,int *nchan,int *iplt, int *cycle,float *xi\
	,float *pa,float *win,int *nb, int *ne, int *nx,int *ip, int *ib,float *bkgd,float *bkgl, char *fname)

 {
	int i,j,k;
	FILE *fd;
	int a;

	if((fd=fopen("datt","rt")) != NULL)
	{
		fscanf(fd,"%d%d%d%d%d%d",jx,npar,ma,nchan,iplt,cycle);
		printf("%d %d %d %d %d\n",*jx,*npar,*ma,*nchan,*iplt);
		for(i=1;i<=*jx;i++) fscanf(fd,"%f",&xi[i]);
		for(i=1;i<=*jx;i++) printf("%f\n",xi[i]);
		for(i=1;i<=*jx;i++) fscanf(fd,"%d%d%d",&nb[i],&ne[i],&nx[i]);
		for(i=1;i<=*jx;i++) printf("%d %d %d\n",nb[i],ne[i],nx[i]);
		for(i=1;i<=*ma;i++) fscanf(fd,"%f",&pa[i]);
		for(i=1;i<=*ma;i++) printf("%f ",pa[i]);
		i=1;
		while(i<=*ma+1){
			a = fgetc(fd);
			if ((a=='N')||(a=='n'))
				{
					win[i] = 0.0;
					i++;
				}
			else if ( ((a>='0')&&(a<='9')) || (a=='.'))
				{
					fseek(fd,-1L,SEEK_CUR);
					fscanf(fd,"%f",&win[i]);
					i++;
				}
		}
		fseek(fd,-1L,SEEK_CUR);
		for(i=1;i<=*ma;i++) printf("\n%f ",win[i]);
		fscanf(fd,"%d",ip);
		printf("\nip=%d \n",*ip);
		if (*ip>0)
		{
			for(i=1;i<=*ip;i++) fscanf(fd,"%d",&ib[i]);
			for(i=1;i<=*ip;i++) printf("%d  ",ib[i]);
		}
		for(i=1;i<=*jx;i++) fscanf(fd,"%f%f",&bkgl[i],&bkgd[i]);
		for(i=1;i<=*jx;i++) printf("\n%f %f",bkgl[i],bkgd[i]);
		fscanf(fd,"%s",fname);
		printf("\n%s\n Okay??\n",fname);

		getchar();
	}
	fclose(fd);
	return;
}

void read_exem_file(int jx,char *fname,int nchan, float **il, float **ie,\
    char *zim)
{
	int i,j,k,ii,jj;
	FILE *fd;
	int nlines;
	char dummy[9],dum[80];
	long dat[8];
	float *x;

	x = vector(1,nchan);
	nlines = nchan/8;
	printf("nchan=%d nlines=%d\n",nchan,nlines);
	if((fd=fopen(fname,"rt")) == NULL)
	{
		printf(" Data file '%s' does not exist. Quitting!\n",fname);
		exit(1);
	}
	else
	{
		for(jj=1;jj<=jx;jj++)
		{
			for(i=0;i<nlines;i++)
			{
				do
				{
					if((j=fgetc(fd))==EOF)
					{
						printf("\nCheck data file '%s'. Quitting.",fname);
						exit(1);
					}
					//printf("%d  ",j);
					//getchar();
				}while(j!='>');
				if((ii=fscanf(fd,"%ld%ld%ld%ld%ld%ld%ld%ld",&dat[0],\
				&dat[1],&dat[2],&dat[3],&dat[4],\
				&dat[5],&dat[6],&dat[7]))!=EOF)
				{
					for (j=1;j<=8;j++) il[jj][i*8+j] =  (float)dat[j-1];
					//printf("i=%d %ld %ld %f %f\n",i,dat[0],dat[7],il[jj][i*8+1],il[jj][i*8+8]);
					//getchar();
				}
				else
				{
					printf("\nNot enough data in '%s'; quitting!!\n",fname);
					exit(1);
				}
			}
				//getchar();
			for(i=0;i<nlines;i++)
			{
				do
				{
					if((j=fgetc(fd))==EOF)
					{
						printf("\nCheck data file '%s'. Quitting.",fname);
						exit(1);
					}
					//printf("%d  ",j);
					//getchar();
				}while(j!='>');
				if((ii=fscanf(fd,"%ld%ld%ld%ld%ld%ld%ld%ld",&dat[0],\
				&dat[1],&dat[2],&dat[3],&dat[4],\
				&dat[5],&dat[6],&dat[7]))!=EOF)
				{
					for (j=1;j<=8;j++) ie[jj][i*8+j] =  (float)dat[j-1];
					//printf("i=%d %ld %ld %f %f\n",i,dat[0],dat[7],ie[jj][i*8+1],il[jj][i*8+8]);
					//getchar();
				}
				else
				{
					printf("\nNot enough data in '%s'; quitting!!\n",fname);
					exit(1);
				}
			}
			if((strcmp(zim,"fs"))==0)
			{
			for(i=0;i<nlines;i++)
			{
				do
				{
					if((j=fgetc(fd))==EOF)
					{
					printf("\n file subtraction [fs] selected.\n");
						printf("\nCheck data file '%s'. Quitting.",fname);
						exit(1);
					}
					//printf("%d  ",j);
					//getchar();
				}while(j!='>');
				if((ii=fscanf(fd,"%ld%ld%ld%ld%ld%ld%ld%ld",&dat[0],\
				&dat[1],&dat[2],&dat[3],&dat[4],\
				&dat[5],&dat[6],&dat[7]))!=EOF)
				{
					for (j=1;j<=8;j++) x[i*8+j] =  (float)dat[j-1];
					//printf("i=%d %ld %ld %f %f\n",i,dat[0],dat[7],x[i*8+1],il[jj][i*8+8]);
					//getchar();
				}
				else
				{
					printf("\n file subtraction [fs] selected.\n");
					printf("\nNot enough data in '%s'; quitting!!\n",fname);
					exit(1);
				}
			}
			for(i=1;i<=nchan;i++) ie[jj][i] -= x[i];
			}
		}
	}
	fclose(fd);
}

void read_ex_em_file(int jx,char *fname1,char *fname2,int nchan,\
 float **il, float **ie,char *zim)
{
	int i,j,k,ii,jj;
	FILE *fd,*fdd;
	int nlines;
	char dummy[9],dum[80];
	long dat[8];
	float *x;


	x = vector(1,nchan);
	nlines = nchan/8;
	printf("nchan=%d nlines=%d\n",nchan,nlines);
	if((fd=fopen(fname1,"rt")) == NULL)
	{
	    printf(" Excitation file '%s' does not exist. Quitting!\n",fname1);
		exit(1);
	}
	else
	{
		for(jj=1;jj<=jx;jj++)
		{
			for(i=0;i<nlines;i++)
			{
				do
				{
					if((j=fgetc(fd))==EOF)
					{
			printf("\nCheck data file '%s'. Quitting.",fname1);
						exit(1);
					}
					//printf("%d  ",j);
					//getchar();
				}while(j!='>');
				if((ii=fscanf(fd,"%ld%ld%ld%ld%ld%ld%ld%ld",&dat[0],\
				&dat[1],&dat[2],&dat[3],&dat[4],\
				&dat[5],&dat[6],&dat[7]))!=EOF)
				{
			for (j=1;j<=8;j++) il[jj][i*8+j] = (float)dat[j-1];
		//printf("i=%d %ld %ld %f %f\n",i,dat[0],dat[7],il[jj][i*8+1],il[jj][i*8+8]);
					//getchar();
				}
				else
				{
		printf("\nNot enough data in '%s'; quitting!!\n",fname1);
					exit(1);
				}
			}

		}
	}
	fclose(fd);
	if((fd=fopen(fname2,"rt")) == NULL)
	{
	    printf(" Emission file '%s' does not exist. Quitting!\n",fname2);
		exit(1);
	}
	else
	{
		for(jj=1;jj<=jx;jj++)
		{
			for(i=0;i<nlines;i++)
			{
				do
				{
					if((j=fgetc(fd))==EOF)
					{
			printf("\nCheck data file '%s'. Quitting.",fname2);
						exit(1);
					}
					//printf("%d  ",j);
					//getchar();
				}while(j!='>');
			if((ii=fscanf(fd,"%ld%ld%ld%ld%ld%ld%ld%ld",&dat[0],\
				&dat[1],&dat[2],&dat[3],&dat[4],\
				&dat[5],&dat[6],&dat[7]))!=EOF)
				{
			for (j=1;j<=8;j++) ie[jj][i*8+j] = (float)dat[j-1];
		//printf("i=%d %ld %ld %f %f\n",i,dat[0],dat[7],ie[jj][i*8+1],ie[jj][i*8+8]);
					//getchar();
				}
				else
				{
		printf("\nNot enough data in '%s'; quitting!!\n",fname2);
					exit(1);
				}
			}

		}
	}
	fclose(fd);
}
void read_exem_melylab_file(int jx,char *fname,int nchan, float **il, float **ie,char *zim)
{
	int i,j,k,ii,jj;
	FILE *fd;
	int nlines;
	char dummy[9],dum[80];
	float dat[8];
	float *x;


	x = vector(1,nchan);
	nlines = nchan/8;
	//printf("nchan=%d nlines=%d\n",nchan,nlines);
	if((fd=fopen(fname,"rt")) == NULL)
	{
		printf(" Data file '%s' does not exist. Quitting!\n",fname);
		exit(1);
	}
	else
	{
		for(jj=1;jj<=jx;jj++)
		{
			do
				{
					k=j;
					if((j=fgetc(fd))==EOF)
					{
						printf("\nCheck data file '%s'. Quitting.",fname);
						exit(1);
					}
					//printf("%d  ",j);
					//getchar();
				}while(j!='X' || k!='U');
			for(i=0;i<nlines;i++)
			{
				if((ii=fscanf(fd,"%d%f%f%f%f%f%f%f%f",&k,&dat[0],\
				&dat[1],&dat[2],&dat[3],&dat[4],\
				&dat[5],&dat[6],&dat[7]))!=EOF)
				{
			for (j=1;j<=8;j++) il[jj][i*8+j] =  dat[j-1];
			//printf("k=%d,i=%d %f %f %f %f\n",k,i,dat[0],dat[7],il[jj][i*8+1],il[jj][i*8+8]);
			//getchar();
				}
				else
				{
					printf("\nNot enough data in '%s'; quitting!!\n",fname);
					exit(1);
				}
			}
				//getchar();
				do
				{
					k=j;
					if((j=fgetc(fd))==EOF)
					{
						printf("\nCheck data file '%s'. Quitting.",fname);
						exit(1);
					}
					//printf("%d  ",j);
					//getchar();
				}while(j!='X' || k!='U');
			for(i=0;i<nlines;i++)
			{
				if((ii=fscanf(fd,"%d%f%f%f%f%f%f%f%f",&k,&dat[0],\
				&dat[1],&dat[2],&dat[3],&dat[4],\
				&dat[5],&dat[6],&dat[7]))!=EOF)
				{
			for (j=1;j<=8;j++) ie[jj][i*8+j] =  dat[j-1];
			//printf("k=%d,i=%d %f %f %f %f\n",k,i,dat[0],dat[7],ie[jj][i*8+1],ie[jj][i*8+8]);
			//getchar();
				}
				else
				{
					printf("\nNot enough data in '%s'; quitting!!\n",fname);
					exit(1);
				}
			}
			if((strcmp(zim,"fs"))==0)
			{
			for(i=0;i<nlines;i++)
			{
				do
				{
					if((j=fgetc(fd))==EOF)
					{
					printf("\n file subtraction [fs] selected.\n");
						printf("\nCheck data file '%s'. Quitting.",fname);
						exit(1);
					}
					//printf("%d  ",j);
					//getchar();
				}while(j!='>');
				if((ii=fscanf(fd,"%f%f%f%f%f%f%f%f",&dat[0],\
				&dat[1],&dat[2],&dat[3],&dat[4],\
				&dat[5],&dat[6],&dat[7]))!=EOF)
				{
					for (j=1;j<=8;j++) x[i*8+j] =  (float)dat[j-1];
					//printf("i=%d %ld %ld %f %f\n",i,dat[0],dat[7],x[i*8+1],il[jj][i*8+8]);
					//getchar();
				}
				else
				{
					printf("\n file subtraction [fs] selected.\n");
					printf("\nNot enough data in '%s'; quitting!!\n",fname);
					exit(1);
				}
			}
			for(i=1;i<=nchan;i++) ie[jj][i] -= x[i];
			}
		}
	}
	fclose(fd);
	//getchar();
}
void covsrt(float **covar,int ma,int lista[],int mfit)
{
	int i=0,j=0;
	float swap=0;

	for (j=1;j<ma;j++) for (i=j+1;i<=ma;i++) covar[i][j]=0.0;
	for (i=1;i<mfit;i++)	for (j=i+1;j<=mfit;j++) 
	{
		if (lista[j] > lista[i])	covar[lista[j]][lista[i]]=covar[i][j];
		else covar[lista[i]][lista[j]]=covar[i][j];
	}
	swap=covar[1][1];
	for (j=1;j<=ma;j++)
	{
		covar[1][j]=covar[j][j];
		covar[j][j]=0.0;
	}
	covar[lista[1]][lista[1]]=swap;
	for (j=2;j<=mfit;j++) covar[lista[j]][lista[j]]=covar[1][j];
	for (j=2;j<=ma;j++)	for (i=1;i<=j-1;i++) covar[i][j]=covar[j][i];
}
void gaussj(float **a,int n,float **b,int m)
{
	int *indxc,*indxr,*ipiv;
	int i=0,icol=0,irow=0,j=0,k=0,l=0,ll=0;
	float big=0,dum=0,pivinv=0;
	indxc=ivector(1,n);
	indxr=ivector(1,n);
	ipiv=ivector(1,n);
	for (j=1;j<=n;j++) ipiv[j]=0;
	for (i=1;i<=n;i++) 
	{
		big=0.0;
		for (j=1;j<=n;j++) if (ipiv[j] != 1) for (k=1;k<=n;k++) 
		{
			if (ipiv[k] == 0) 
			{
				if (fabs(a[j][k]) >= big) 
				{
					big=fabs(a[j][k]);
					irow=j;
					icol=k;
				}
			} 
			else if (ipiv[k] > 1){ printf("GAUSSJ: Singular Matrix-1");exit(1);}
		}
		++(ipiv[icol]);
		if (irow != icol)
		{
				for (l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l]);
				for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l]);
		}
		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol][icol] == 0.0) a[icol][icol] += 1.e-20;
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=1;l<=n;l++) a[icol][l] *= pivinv;
		for (l=1;l<=m;l++) b[icol][l] *= pivinv;
		for (ll=1;ll<=n;ll++) if (ll != icol) 
		{
			dum=a[ll][icol];
			a[ll][icol]=0.0;
			for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
			for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
		}
	}
	for (l=n;l>=1;l--) if (indxr[l] != indxc[l]) for (k=1;k<=n;k++) SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	free_ivector(indxc,1,n);
	free_ivector(indxr,1,n);
	free_ivector(ipiv,1,n);
}

void core_job(float *yex,float *yval,float *ycal,int jx,int ndata,int npar,int ma, int mfit, int *nb, int *ne, int *nx,int *lista,float *a,float *xi, float **alpha,float *beta,float *chisq,float *sigma)
{
	int i=0,j=0,k=0;
	float ymod=0,wt=0,sig2i=0,dy=0;
	float **dyda;
	int nch=0;
	int jj=0,kount=0;
	dyda = matrix(1,ma,1,ndata);
	for(j=1;j<=ma;j++) for(i=1;i<=ndata;i++) dyda[j][i] = 0.;
	for(j=1;j<=ndata;j++) if(sigma[j]<=0.)
	{
		printf( "Error. j=%d, yval=%f,sigma=%f\n",j,yval[j],sigma[j]);
		sigma[j]=1.;
	}
	convolute(yex,jx,ndata,npar,ma,mfit,nb,ne,nx,lista,a,xi,ycal,dyda);
	for(jj=1;jj<=jx;jj++)
	{
		nch = ne[jj]-nb[jj]+1;
		for(i=kount+1;i<=kount+nch;i++)
		{
			if(i>kount+nx[jj])
			{
				sig2i=1.0/(sigma[i]*sigma[i]);
				dy=yval[i]-ycal[i];
				for (j=1;j<=mfit;j++)
				{
					wt=dyda[lista[j]][i]*sig2i;
					for (k=1;k<=j;k++)
					alpha[j][k] += wt*dyda[lista[k]][i];
					if(j==k && alpha[j][k]==0.) printf("(core-job) WARNING! Parameter %d to be checked.\n",lista[k]);
					beta[j] += dy*wt;
				}
				(*chisq) += dy*dy*sig2i;
			}
		}
		kount += nch;
	}
	free_matrix(dyda,1,ma,1,ndata);
	return;
}
void mrqcof(float *yex,float *yem,float *ycal,int jx,int ndata,int npar,int ma,	int mfit, int *nb, int *ne, int *nx,int *lista,float *a,float *xi, float **alpha,float *beta,float *chisq,float *sigma)
{
	int k=0,j=0,i=0;
	for (j=1;j<=mfit;j++)
	{
		for (k=1;k<=j;k++) alpha[j][k]=0.0;
		beta[j]=0.0;
	}
	*chisq=0.0;
	core_job(yex,yem,ycal,jx,ndata,npar, ma,mfit,nb,ne,nx,lista,a,xi,alpha,beta,chisq,sigma);
	for (j=2;j<=mfit;j++) for (k=1;k<=j-1;k++) alpha[k][j]=alpha[j][k];
	return;
}


void mrqmin(float *yex,float *yem,float *ycal,int jx,int ndata,int npar,int ma,	int mfit, int *nb, int *ne, int *nx,int *lista,float *a,float *xi, float **covar,float **alpha,float *chisq,float *alamda,float *sigma)
{
	int k=0,kk=0,j=0,ihit=0,i=0;
	static float *da,*atry,**oneda,*beta,ochisq;
	float x=0;
	if (*alamda < 0.0) 
	{
		oneda=matrix(1,mfit,1,1);
		atry=vector(1,ma);
		da=vector(1,ma);
		beta=vector(1,ma);
		kk=mfit+1;
		for (j=1;j<=ma;j++)
		{
			ihit=0;
			for (k=1;k<=mfit;k++) if (lista[k] == j) ihit++;
			if (ihit == 0) lista[kk++]=j;
			else if (ihit > 1){ printf("Bad LISTA permutation in MRQMIN-1");exit(1);}
		}
		if (kk != ma+1) {printf("Bad LISTA permutation in MRQMIN-2");exit(1);}
		*alamda=0.01;
		mrqcof(yex,yem,ycal,jx,ndata,npar,ma,mfit,nb,ne,nx,lista,a,xi,alpha,beta,chisq,sigma);
		ochisq=(*chisq);
	}
	for (j=1;j<=mfit;j++)
	{
		for (k=1;k<=mfit;k++) covar[j][k]=alpha[j][k];
		covar[j][j]=alpha[j][j]*(1.0+(*alamda));
		oneda[j][1]=beta[j];
	}
	gaussj(covar,mfit,oneda,1);
	for (j=1;j<=mfit;j++)
		da[j]=oneda[j][1];
	if (*alamda == 0.0)
	{
		for (j=1;j<=ma;j++) atry[j]=a[j];
		mrqcof(yex,yem,ycal,jx,ndata,npar,ma,mfit,nb,ne,nx,lista,atry,xi,alpha,beta,chisq,sigma);
		covsrt(covar,ma,lista,mfit);
		free_matrix(oneda,1,mfit,1,1);
		free_vector(atry,1,ma);
		free_vector(da,1,ma);
		free_vector(beta,1,ma);
		return;
	}
	for (j=1;j<=ma;j++) atry[j]=a[j];
	for (j=1;j<=mfit;j++)	atry[lista[j]] = a[lista[j]]+0.1*da[j];
	mrqcof(yex,yem,ycal,jx,ndata,npar,ma,mfit,nb,ne,nx,lista,atry,xi,alpha,beta,chisq,sigma);
	if (*chisq < ochisq)
	{
		if(*alamda>1e-10) *alamda *= 0.1;
		ochisq=(*chisq);
		for (j=1;j<=mfit;j++)	a[lista[j]]=atry[lista[j]];
	}
	else
	{
		if(*alamda<1000.) *alamda *= 10.0;
		*chisq=ochisq;
	}
	return;
}


void convolute(float *yex,int jx,int ndata,int npar,int ma,int mfit,\
	int *nb, int *ne, int *nx,int *lista,float *a,float *xi,float *ycal,\
	float **dyda)
{
	int i,j,k;
	int nch,jin,count=1;
	float **pd;
	float **b;
	//float *pd1,*pd2,*pd3,*pd4,*pd5;
	float *yls,*yll,delta=0.,dt,b2sq,delta_in;
	float *temp1,*temp2,*valu1,*valu2,amp1,amp2,tau1,tau2;
	float f[3],g[3],*v1,*v2;
	char shift_option='y';

	b = matrix(1,jx,1,npar);
	para_convert(jx,ma,b,a);
	//for(i=1;i<=jx;i++) printf(" b's aren %f %f %f\n", b[i][1],b[i][2],b[i][3]);
	//getchar();
	for(jin=1;jin<=jx;jin++)
	{
		nch = ne[jin]-nb[jin]+1;
		valu1= vector(1,nch);
		temp1= vector(1,nch);
		valu2= vector(1,nch);
		temp2= vector(1,nch);
		yls = vector(1,nch);
		v1 = vector (1,nch);
		v2 = vector (1,nch);
		//pd1 = vector(1,nch);
		//pd2 = vector(1,nch);
		//pd3 = vector (1,nch);
		//pd4 = vector(1,nch);
		//pd5 = vector (1,nch);
		//for(i=1;i<=nch;i++){
		//    pd1[i]=0.;pd2[i]=0.;pd3[i]=0.;pd4[i]=0.;pd5[i]=0.;
		//}
		pd = matrix(1,npar,1,nch);
		for(i=1;i<=npar;i++)
			for(j=1;j<=nch;j++)
				pd[i][j]=0.0;
		amp1= b[jin][1];
		tau1= b[jin][2];
		amp2= b[jin][3];
		tau2= b[jin][4];
		dt = xi[jin];
		f[1] = 0.5*dt* amp1;
		f[2] = (float) exp ( (double) (- dt/tau1));
		g[1] = 0.5*dt* amp2;
		g[2] = (float) exp ( (double) (- dt/tau2));
		//printf("mfit=%d,nch=%d amp1=%f tau1=%f delta=%f dt=%f\n",mfit,nch,amp1,tau1,delta,dt);
		//getchar();
		for(i=1;i<=nch;i++) yls[i] = yex[i+count-1];
		//for(k=1;k<=nch;k++){ printf("%d %f\n",k,yls[k]);getchar();}
		for(i=mfit+1;i<=ma;i++){
			if(lista[i]==jin*3+2)shift_option = 'n';
			 //printf("i=%d,lista[%d]=%d,  %c\n",i,i,lista[i],shift_option);
			 //getchar();
		}
		if(shift_option == 'y')
		{
			yll = vector(1,nch);
			for(i=1;i<=nch;i++) yll[i] = yex[i+count-1];
			delta = b[jin][5]/dt;
			delta_in = 0.1 + delta;
			time_shift( delta_in,yll, yls,nch);
			for(j=1; j<=nch; j++) v1[j] = f[1] * yls[j];
			for(j=1; j<=nch; j++) v2[j] = g[1] * yls[j];
			temp1[1] = v1[1];
			for(j=2;j<=nch;j++) temp1[j] = f[2]*(temp1[j-1]+v1[j-1])+v1[j];
			temp2[1] = v2[1];
			for(j=2;j<=nch;j++) temp2[j] = g[2]*(temp2[j-1]+v2[j-1])+v2[j];
			time_shift( delta,yll,yls,nch);
			for(j=1; j<=nch; j++) v1[j] = f[1] * yls[j];
			for(j=1; j<=nch; j++) v2[j] = g[1] * yls[j];
			valu1[1] = v1[1];
			for(j=2;j<=nch;j++) valu1[j] = f[2]*(valu1[j-1]+v1[j-1])+v1[j];
			valu2[1] = v2[1];
			for(j=2;j<=nch;j++) valu2[j] = g[2]*(valu2[j-1]+v2[j-1])+v2[j];
			free_vector(yll,1,nch);
		}
		else
		{
			for(j=1;j<=nch;j++){
				temp1[j]=0.;
				temp2[j]=0.;
			}
			for(j=1; j<=nch; j++) v1[j] = f[1] * yls[j];
			valu1[1] = v1[1];
			for(j=2;j<=nch;j++)
				valu1[j] = f[2] * (valu1[j-1]+v1[j-1])+v1[j];
			for(j=1; j<=nch; j++) v2[j] = g[1] * yls[j];
			valu2[1] = v2[1];
			for(j=2;j<=nch;j++)
				valu2[j] = g[2] * (valu2[j-1]+v2[j-1])+v2[j];
		}
		for(j=1;j<=nch;j++)
		{
			if(shift_option=='y') pd[5][j] = (temp1[j]+temp2[j]-valu1[j]-valu2[j])/(0.1*dt);
			else pd[5][j]=0.;
			pd[3][j] = valu1[j] / (amp1+1.e-20);
			pd[4][j] = valu2[j] / (amp2+1.e-20);
		}
		pd[1][1]=0.;
		pd[2][1]=0.;
		for(j=2;j<=nch;j++){
			pd[1][j] = f[2]*(pd[1][j-1]+dt*(valu1[j-1]+v1[j-1])/(tau1*tau1));
			pd[2][j] = g[2]*(pd[2][j-1]+dt*(valu2[j-1]+v2[j-1])/(tau2*tau2));
		}
		for(k=count,j=1;k<=(count+nch-1);k++,j++)
		{
			ycal[k] = valu1[j]+valu2[j];
			for(i=1;i<=mfit;i++)
			{
				if(lista[i]==1) dyda[lista[i]][k]=pd[1][j];
				if(lista[i]==2) dyda[lista[i]][k]=pd[2][j];
				if(lista[i]==3+(jin-1)*3) dyda[lista[i]][k]=pd[3][j];
				if(lista[i]==4+(jin-1)*3) dyda[lista[i]][k]=pd[4][j];
				if(lista[i]==5+(jin-1)*3) dyda[lista[i]][k] = pd[5][j];
			}
		}
		count += nch;
		//printf(" jin=%d, count=%d\n",jin,count);
		//getchar();
	free_vector(yls,1,nch);
	free_vector(valu1,1,nch);
	free_vector(v1,1,nch);
	free_vector(temp1,1,nch);
	free_vector(valu2,1,nch);
	free_vector(v2,1,nch);
	free_vector(temp2,1,nch);
	//free_vector(pd1,1,nch);
	//free_vector(pd2,1,nch);
	//free_vector(pd3,1,nch);
	//free_vector(pd4,1,nch);
	//free_vector(pd5,1,nch);
	free_matrix(pd,1,npar,1,nch);
	}
	free_matrix(b,1,jx,1,npar);
	return;
}
void time_shift(float d,float *yll,float *yls,int nch)
{
	int i=0,j=0,k=0;
	int nnb=0,ns=0;
	double x=0,y=0;
	float cshift=0, shift=0, *tem;
	for(i=1;i<=nch;i++) yls[i]=yll[i];
	if(d==0.0) return;
	tem = vector(1,nch);
	if(d>0.)
	{
		x = modf((double) d,&y);
		shift = (float) x;
		ns = (int) y;
		if(ns>=1)
		{
			for(i=1;i<=ns;i++) tem[i]=0.;
			for(i=1;i<=nch-ns;i++) tem[i+ns] = yll[i];
			for(i=1;i<=nch;i++) yls[i] = tem[i];
		}
		if(shift!= 0.)
		{
			cshift = 1. - shift;
			tem[ns+1] = yls[ns+1]*cshift;
			for(i=ns+2;i<=nch;i++) tem[i]=yls[i]*cshift+yls[i-1]*shift;
			for(i=1;i<=nch;i++) yls[i] = tem[i];
		}
	}
	else
	{
		d *= -1.;
		x = modf((double) d,&y);
		shift = (float) x;
		ns = (int) y;
		if(ns>=1)
		{
			for(i=nch;i>=nch-ns+1;i--) tem[i]=0.;
			for(i=1;i<=nch-ns;i++) tem[i] = yll[i+ns];
			for(i=1;i<=nch;i++) yls[i] = tem[i];
		}
		if(shift!=0.)
		{
			cshift = 1. - shift;
			tem[nch-ns] = yls[nch-ns]*shift;
			for(i=1;i<=nch-ns-1;i++) tem[i]=yls[i]*cshift+yls[i+1]*shift;
			for(i=1;i<=nch;i++) yls[i] = tem[i];
		}
	}
	free_vector(tem,1,nch);
	return;
}


void para_convert(int jx,int ma,float **b, float *a)
/*
	parameters a[] are as follows:
	a[1] lifetime 1;
	a[2] -> lifetime 2;
	a[3],a[4],a[5]  are amp1, amp2 and delta for jx=1;
	a[6],a[7], a[8] are amp1, amp2 and delata for jx=2; and so on.

	After conversion,

	b[][1] & b[][3] are amplitudes;
	b[][2] & b[][4] are the lifetimes in nanoseconds;
	b[][5] are the time shifts in nanoseconds.
*/

{
	int i;
	if(a[1]<0.) a[1] = (float) fabs((double) a[1]);
	if(a[2]<0.) a[2] = (float) fabs((double) a[2]);
	for(i=1;i<=ma;i++)
	{
		if(a[i]<par.low[i]) a[i]=par.low[i];
		if(a[i]>par.high[i]) a[i]=par.high[i];
	}
	for(i=1;i<=ma;i++) if(a[i]==0.) a[i]=1.e-10;
	for (i=1;i<=jx;i++)
	{
		   b[i][2] = a[1];
		   b[i][4] = a[2];
		   b[i][1] = a[(i-1)*3+3];
		   b[i][3] = a[(i-1)*3+4];
		   b[i][5] = a[(i-1)*3+5];
	}
	return;
}

void randomize_amps(int jx,int ma, float *a, float *pa)
{
	int i,j,k;
	float x;

		for(i=1;i<=ma;i++)a[i]=pa[i];
		for(i=1;i<=jx;i++)
		{
		k = rand();
		x = (float)k;
		while(x>1.) x /=2.;
		a[3+(i-1)*3] *=  x;
		k = rand();
		x = (float)k;
		while(x>1.) x /=2.;
		a[4+(i-1)*3] *=  x;
		return;
		}
}

void find_lifetime(int jx,int nchan,int *nb,int *ne,float *xi,\
    float **ie,float *pa)
{
	int i,j,k,nmax;
	float x,y;
	int jj=1;
	nmax=find_max_chan(jj,nchan,ie);
	x=ie[jj][nmax]/2.;
	i=nmax;
	while(x<ie[jj][i]){;i++;};
	y=(float)(i-nmax)*xi[1]/0.693;
	pa[1]=0.5*y;
	pa[2]=2.*y;
}

void scale_amps(int jx,int nchan,int *nb,int *ne,int *nx,float *xi,float **il,//
	float **ie,float *pa)
{
	int i,j,k;
	float sum,em_max;
	int nmax;
	float tau;

		for(i=1;i<=jx;i++)
		{
			nmax = find_max_chan(i,nchan,ie);
			nx[i] = nmax-nb[i]+nx[i];
			if(nx[i]<0){printf("i=%d,nmax=%d nb[i]=%d nx[i]=%d: choose 'nx' properly!\n",//
				i,nmax,nb[i],nx[i]);exit(0);}
			printf("i=%d nmax=%d nb[i]=%d nx[i]=%d\n",i,nmax,nb[i],nx[i]);
			em_max = ie[i][nmax];
			tau=pa[1]*pa[3+(i-1)*3]+pa[2]*pa[4+(i-1)*3];
			sum=0.;
			for(j=nb[i];j<=ne[i];j++)
				if(j<=nmax) sum +=  il[i][j]*(nmax-j)*xi[i]/tau;
			pa[3+(i-1)*3] *= em_max/(sum);
			pa[4+(i-1)*3] *= em_max/(sum);
		}
	//printf("pa's %f %f %f\n",pa[4],pa[5],pa[6]);
	//getchar();
}

void scale_amps_batch(int jx,int nchan,int *nb,int *ne,int *nx,float *xi,float **il,//
	float **ie,float *pa)
{
	int i,j,k;
	float sum,em_max;
	int nmax,nmx;

		for(i=1;i<=jx;i++)
		{
		    nmax=find_max_chan(i,nchan,il);
		    do
			{
			    nmax--;
			}while(il[i][nmax]>5 && nmax>1);
			if(nmax<10)nb[i]=nmax;
			if(nmax<20)nb[i]=nmax-10;
			if(nmax>20) nb[i]=nmax-15;
			nmax = find_max_chan(i,nchan,ie);
			nmx=nmax;
			do
			    {
				nmx--;
			    }while(ie[i][nmx]>5 && nmx>1);
			if(nmx<nb[i])nb[i]=nmx;
			nmx=nmax;
			do
			{
			    nmx++;
			    }while(ie[i][nmx]>10 && nmx<nchan-1);
			ne[i]=nmx;
			nx[i] = nmax-nb[i]+nx[i];
			if(nx[i]<0){printf("i=%d,nmax=%d nb[i]=%d nx[i]=%d: choose 'nx' properly!\n",//
				i,nmax,nb[i],nx[i]);exit(0);}
			printf("i=%d nmax=%d nb[i]=%d ne[i]=%d nx[i]=%d\n",i,nmax,nb[i],ne[i],nx[i]);
			em_max = ie[i][nmax];
			sum=0.;
			for(j=nb[i];j<=ne[i];j++)
				if(j<=nmax) sum +=  il[i][j];
			pa[3+(i-1)*3] *= em_max/(sum*0.5*xi[i]);
			pa[4+(i-1)*3] *= em_max/(sum*0.5*xi[i]);
		    nmx=find_max_chan(i,nchan,il);
		    k=nmx;
		    do
			{
			    k++;
			}while(il[i][k]>(0.5*il[i][nmx]));
		    if(k>nmax) pa[5+(i-1)*3]=(nmax-k+1)*xi[i];
		}
}

int find_max_chan(int jj, int nchan,float **ie)
{
	int i,j,k;
	float xmax=0.;

	for(i=1;i<=nchan;i++)
	{
		if(ie[jj][i]>xmax){ k=i;xmax=ie[jj][i];}
	}
	return k;
}

void out_screen(int jx,int ma, float *a, float **alpha)
{
	int i;
	float b[PARMS], sd[PARMS],c=0.,d;

	for(i=1;i<=ma;i++)sd[i];
	for(i=1;i<=ma;i++) if(alpha[i]>0)sd[i]=(float) sqrt((double)(1./alpha[i][i]));
	for(i=1;i<=ma;i++)b[i]=a[i];
	for(i=1;i<=jx;i++)
	{
		if(a[(i-1)*3+3]>0.)c+=a[(i-1)*3+3];
		if(a[(i-1)*3+4]>0.)c+=a[(i-1)*3+4];
		b[(i-1)*3+3] = a[(i-1)*3+3]/c;
		b[(i-1)*3+4] = a[(i-1)*3+4]/c;
	}


	 printf("\n\tfluor. lifetime1  = %2.3f ns; \n",b[1]);
	 for(i=0;i<=jx;i++)printf((i==0)?"\t\t amplitudes: ": "  %2.3f",b[3+(i-1)*3]);
	 printf("\n\tfluor. lifetime2  = %2.3f ns;\n",b[2]);
	 for(i=0;i<=jx;i++)printf((i==0)?"\t\t amplitudes: ": "  %2.3f",b[4+(i-1)*3]);
	 for(i=0;i<=jx;i++)
	 {
		d = b[1]*b[3+(i-1)*3] + b[2]*b[4+(i-1)*3];
	 printf((i==0)?"\n\n\tAverage Lifetimes: ":"  %2.3f",d);
	 }
	 for(i=0;i<=jx;i++)printf((i==0)?"\n\t\t shifts: ": "  %2.3f",b[5+(i-1)*3]);

	 for(i=3;i<=ma-jx;i++) if(b[i]<0.){
		printf("\n\nWARNING: Negative amp. returned.  \n\
			Model should be appropriate!\n");
		}
}

void out_file(FILE *fres,FILE *fres2,int jx,int cycle,int iter,float chisq,int ma,//
	float *a,float **alpha,float *yex,float *yem, float *ycal,float *res,float *acorr,//
	int ntot,int mfit,int *lista,int c_c)
{
	int i,j;
	float b[PARMS], sd[PARMS],c=0.,d;

	for(i=1;i<=ma;i++)sd[i]=0.;
	for(i=1;i<=mfit;i++) if(alpha[i]>0)sd[lista[i]]=(float) sqrt((double)(1./alpha[i][i]));
	for(i=1;i<=ma;i++)b[i]=a[i];
	for(i=1;i<=jx;i++)
	{
		if(a[(i-1)*3+3]>0.)c+=a[(i-1)*3+3];
		if(a[(i-1)*3+4]>0.)c+=a[(i-1)*3+4];
		b[(i-1)*3+3] = a[(i-1)*3+3]/c;
		b[(i-1)*3+4] = a[(i-1)*3+4]/c;
		sd[(i-1)*3+3] = sd[(i-1)*3+3]/c;
		sd[(i-1)*3+4] = sd[(i-1)*3+4]/c;
	}
	fprintf(fres,"\nCycle=%d, Iter=%d, chisq=%2.3f\n",cycle,iter,chisq);
	 fprintf(fres,"\nfluor. lifetime1  = %2.3f ns; \n",b[1]);
	 for(i=0;i<=jx;i++)fprintf(fres,(i==0)?"\t\t amplitudes: ": "  %2.3f",b[3+(i-1)*3]);
	 fprintf(fres,"\nfluor. lifetime2  = %2.3f ns;\n",b[2]);
	 for(i=0;i<=jx;i++)fprintf(fres,(i==0)?"\t\t amplitudes: ": "  %2.3f",b[4+(i-1)*3]);
	 for(i=0;i<=jx;i++)
	 {
		d = b[1]*b[3+(i-1)*3] + b[2]*b[4+(i-1)*3];
	 fprintf(fres,(i==0)?"\n\n\tAverage Lifetimes: ":"  %2.3f ns",d);
	 }
	 for(i=0;i<=jx;i++)fprintf(fres,(i==0)?"\n\t\t Shifts: ": "  %2.3f ns",b[5+(i-1)*3]);
	fprintf(fres,"\nStd.Devs:\n");
	 fprintf(fres,"\n lifetime1  = %2.3f ns; \n",sd[1]);
	 for(i=0;i<=jx;i++)fprintf(fres,(i==0)?"\t\t amplitudes: ": "  %2.3f",sd[3+(i-1)*3]);
	 fprintf(fres,"\n lifetime2  = %2.3f ns;\n",sd[2]);
	 for(i=0;i<=jx;i++)fprintf(fres,(i==0)?"\t\t amplitudes: ": "  %2.3f",sd[4+(i-1)*3]);
	 fprintf(fres,"\n------------------------------------------------\n");

	fprintf(fres2,"Cycle %d\nIter=%d, chisq=%2.3f\n",cycle,iter,chisq);
	fprintf(fres2,"\t\tfluor. lifetime1  = %2.3f ns; amp1=%2.3f\n",b[1],b[3]);
	fprintf(fres2,"\t\tfluor. lifetime2  = %2.3f ns; amp2=%2.3f\n",b[2],b[4]);
	fprintf(fres2,"\n\t  chan\t  ex\t    em\t     cal em\t    res\t   acorr\n");
	for(i=1;i<=ntot;i++) fprintf(fres2,"%10d %10.2f %10.2f %10.2f %10.2f %10.2f \n",\
	   i, yex[i],yem[i],ycal[i],res[i],acorr[i]);
		fflush(fres2);
	int n = c_c;
	//FILE *tempDataFile;
  char FileName[7] = "result2" ;
  //tempDataFile = fopen(FileName,"r");
  plot(FileName,n);
  
}

void simulate(float *yex,float *ycal,int jx,int ndata,int npar,int mfit,//
	int *nb,int *ne,int *nx,int *lista,float *a,float *xi)
{
		FILE *fsim;
		int i,j,k;
		char fname[12];
		int ma;
		float **dyda,x[513],y[513],p;

		ma = npar;
		dyda = matrix(1,ma,1,ndata);

		for(i=1;i<=512;i++)
		{
			if(i<=ndata){x[i]=yex[i];y[i]=yex[i];}
			else {x[i]=0.0;}
		}
		p=0.0;
		for(i=1;i<=512;i++) if(x[i]>p) { p=x[i]; k=i; }
		for(i=k+2;i<=510;i++)
		{
			p=0.0;
			for(j=-2;j<=2;j++) p += x[i+j];
			y[i] = p/5.;
		}
		for(i=511;i<=512;i++) y[i]=0.;
		for(i=1;i<=512;i++) x[i]=y[i];
		ndata=512;
		ne[1]=512;
		nb[1]=1;
		convolute(x,jx,ndata,npar,ma,mfit,nb,ne,nx,lista,a,xi,y,dyda);
		printf("\nGive filename for simulated data:");
		scanf("%s",fname);
		getchar();
		printf("%s\n",fname);
		if((fsim=fopen(fname,"w")) != NULL)
		{
			for(i=1;i<=512;i++)
			{
				//printf("%d\t%f\t%f\n",i,x[i],y[i]);
				fprintf(fsim,"%d\t%f\t%f\n",i,x[i],y[i]);
			}
		}
		else{ printf(" Error in opening file"); exit(1);}
		fclose(fsim);
		printf(" \nFile %s created. \n",fname);
		free_matrix(dyda,1,ma,1,ndata);
}
void invert_data(int jx,int nchan,float **il,float **ie)
{
	int i,j,k,ii,jj,iii,jjj;
	int nmax;
	float imax,imax_half,imax_tenth,*x;
	int ch;

	for(i=1;i<=jx;i++){
	nmax=find_max_chan(i,nchan,ie);
	imax_half=0.5*ie[i][nmax];
	imax_tenth=0.1*ie[i][nmax];
	ii=0;
	iii=0;
	for(j=1;j<=nchan;j++){
		if(ii==0 && ie[i][j]>imax_half) {ii=j-1;jj=j+1;}
		else if(ii!=0 && ie[i][j]>imax_half)jj=j+1;
		}
	for(j=1;j<=nchan;j++){
		if(iii==0 && ie[i][j]>imax_tenth) {iii=j-1;jjj=j+1;}
		else if(iii!=0 && ie[i][j]>imax_tenth)jjj=j+1;
		}
	//printf("\n nmax=%d,ii=%d,jj=%d iii=%d jjj=%d\n",nmax,ii,jj,iii,jjj);
	//getchar();
	//if((nmax-ii)==(jj-nmax)||(nmax-iii)==(jjj-nmax)){
	  //  printf("\n Decay has either large scatter or intense short lifetime.\n Is the data inverted [Y/N]?:");
		//ch=getchar();
	//if(ch=='Y'||ch=='y')ii--;
	//}
	if((nmax-ii)<(jj-nmax)&&(nmax-iii)<(jjj-nmax)){ /* normal data */
		}
	else if((nmax-ii)>(jj-nmax)&&(nmax-iii)>(jjj-nmax)){    /* invert data */
		x=vector(1,nchan);
		for(j=1;j<=nchan;j++)x[j]=il[i][j];
		for(j=1;j<=nchan;j++)il[i][j]=x[nchan-j+1];
		for(j=1;j<=nchan;j++)x[j]=ie[i][j];
		for(j=1;j<=nchan;j++)ie[i][j]=x[nchan-j+1];
		free_vector(x,1,nchan);
		}
	else{
		printf("\n Decay has either large scatter or intense short lifetime.\n Is the data inverted [Y/N]?:");
		ch=getchar();
		getchar();
		if(ch=='Y'||ch=='y'){
		   x=vector(1,nchan);
		   for(j=1;j<=nchan;j++)x[j]=il[i][j];
		   for(j=1;j<=nchan;j++)il[i][j]=x[nchan-j+1];
		   for(j=1;j<=nchan;j++)x[j]=ie[i][j];
		   for(j=1;j<=nchan;j++)ie[i][j]=x[nchan-j+1];
		   free_vector(x,1,nchan);
		   }
		}
	}
}