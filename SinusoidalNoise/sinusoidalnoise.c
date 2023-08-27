#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include <stdbool.h>

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

//===================================================================//

#define L 32           //Length of simulation box
#define Lh (L/2)        //Half Length of simulation box
#define LS (L*L)        //Area of simulation box
#define N (LS*2)         //number of swimmers
#define ROW0 (N/(1.0*L*L))  //Total density of the particles

#define R 1             //Radius of influence
#define V0 0.3        //Self-propulsion velocity: 0.003<V0<0.3:

#define ETA 0.2    //Noise Amplitude
#define dt 1            //time step
#define PI 3.14159265
#define SAT_TIME 100000  
#define TIME 100000
#define TTIME (SAT_TIME+TIME)
#define StepSize 10

static long iseed=-99999999;
float ran2(long *idum);

//====================================================================//

#define GNUPLOT "gnuplot -persist"

//====================================================================//


int main(void) 
{    
  FILE *gp,*fs;
//====================================================================//
  FILE *fp,*fp1,*fp2,*fp3;
  
  static char FILE_SNAP[30],FILENAME[30],FILENAME1[30],FILENAME2[30],FILENAME3[30];
  
  static long i,j,k,l,m,time,np,list[N],box,cbox,lbox,tbox,rbox,bbox,ltbox,rtbox,lbbox,rbbox,BOX[LS][N],NBOX[LS],frame;
  static int flag,n,X,Y,XCL,XCR,XLL,XLR,XRL,XRR,YCB,YCT,YBB,YBT,YTB,YTT;
  bool flag1=false;
  int count=0;
  double eta;
  static double x[N],y[N],theta[N][2],dx[N],dy[N];  
  static double PHI,thetaj,avgtheta,rantheta,Rtheta,phi,PHI2,PHI4,CHI;
  static double dr,drc,avgdx,avgdy,sx,sy,ln;
  static double vx,vy,va,va1,va2,va4,nsamp,BC,sigma;
  
  printf("For L=%d of N=%d ROW=%lf V0=%f\n",L,N,ROW0,V0);
  
  sprintf(FILENAME1,"PHI-eta_L%d.dat",L);//avearge velocity vs noise
  fp1=fopen(FILENAME1,"a");
  
  sprintf(FILENAME2,"BC-eta_L%d.dat",L);//Binder Cumulant
  fp2=fopen(FILENAME2,"a");
  
  sprintf(FILENAME3,"CHI-eta_L%d.dat",L);//Binder Cumulant
  fp3=fopen(FILENAME3,"a");
  
  sprintf(FILENAME,"TS_L%d_eta%5.3f.dat",L,ETA);//Time Series
  fp=fopen(FILENAME,"w");
      
  //============== INITIALIZATION ================
  for(i=0;i<LS;i++)
    {
      NBOX[i]=0; //No of the particles in the box i
      for(j=0;j<N;j++)BOX[i][j]=0; //index of the particles in the box i
    }
  
  for(i=0;i<N;i++)
    {
      x[i]=ran2(&iseed)*L;//x position
      y[i]=ran2(&iseed)*L;//y position
      theta[i][0]=ran2(&iseed)*2*PI-PI;//oriented angle
      theta[i][1]=0;//to avoid garbage value...
      list[i]=0;
    }
  //=====================================================
  for(i=0;i<N;i++)//loop over the particles
    {
      l=x[i];
      m=y[i];
      box=m*L+l;
      BOX[box][NBOX[box]]=i; //which particles are in the box
      NBOX[box]++;
    }
  
  //=====================================================
  frame=0;
  nsamp=0;
  PHI=0.0;
  PHI2=0.0;
  PHI4=0.0;
  for(time=0;time<=TTIME;time++)//TTIME 
    {
   
   if(count==StepSize){
	count=0;
	flag1=!flag1;
   }
   if(flag1){
	eta=ETA;
   }
   else{
	eta=-ETA;
   }
   count++;
   printf(" %1.1f\n",eta);

      for(i=0;i<N;i++) //Loop over swimmers for orientation update
	{
	  X=x[i];
	  Y=y[i];
	  cbox=Y*L+X;
	  //========================================//	
	  if(X==0)
	    lbox=cbox+L-1;
	  else
	    lbox=cbox-1;
	  
	  if(Y==(L-1))
	    tbox=cbox+L-L*L;
	  else
	    tbox=cbox+L;
	  
	  if(X==(L-1))
	    rbox=cbox+1-L;
	  else
	    rbox=cbox+1;
	  
	  if(Y==0)
	    bbox=cbox-L+L*L;
	  else
	    bbox=cbox-L;
	  //========================================//
	  if(X==0)
	    ltbox=tbox+L-1;
	  else
	    ltbox=tbox-1;;
	  
	  if(X==(L-1))
	    rtbox=tbox+1-L;
	  else
	    rtbox=tbox+1;
	  
	  if(X==0)
	    lbbox=bbox+L-1;
	  else
	    lbbox=bbox-1;;
	  
	  if(X==(L-1))
	    rbbox=bbox+1-L;
	  else
	    rbbox=bbox+1;
	  
	  
	  //---------------------------List the number of swimmers-------------------------------------
	  m=0;
	  for(j=0;j<NBOX[lbox];j++) //To find swimmers in nine boxes
	    {
	      list[m]=BOX[lbox][j];
	      m++;
	    }
	  for(j=0;j<NBOX[tbox];j++) //To find swimmers in nine boxes
	    {
	      list[m]=BOX[tbox][j];
	      m++;
	    }
	  for(j=0;j<NBOX[rbox];j++) //To find swimmers in nine boxes
	    {
	      list[m]=BOX[rbox][j];
	      m++;
	    }
	  for(j=0;j<NBOX[bbox];j++) //To find swimmers in nine boxes
	    {
	      list[m]=BOX[bbox][j];
	      m++;
	    }
	  
	  //==============================================================//
	  
	  for(j=0;j<NBOX[ltbox];j++) //To find swimmers in nine boxes
	    {
	      list[m]=BOX[ltbox][j];
	      m++;
	    }
	  for(j=0;j<NBOX[rtbox];j++) //To find swimmers in nine boxes
	    {
	      list[m]=BOX[rtbox][j];
	      m++;
	    }
	  for(j=0;j<NBOX[lbbox];j++) //To find swimmers in nine boxes
	    {
	      list[m]=BOX[lbbox][j];
	      m++;
	    }
	  for(j=0;j<NBOX[rbbox];j++) //To find swimmers in nine boxes
	    {
	      list[m]=BOX[rbbox][j];
	      m++;
	    }
	  
	  //==============================================================//
	  
	  avgdx=cos(theta[i][0]);//x component	
	  avgdy=sin(theta[i][0]);//y component
	  
	  np=1;//number of particles
	  
	  for(j=0;j<m;j++) //Loop over other swimmers in the list 
	    {
	      dr=sqrt(pow((x[i]-x[list[j]]),2)+pow((y[i]-y[list[j]]),2));//distance between the listed particles
	      if(dr>Lh)
		dr=L-dr; //corrected dr
	      if(dr<=R)
		{
		  avgdx+=cos(theta[list[j]][0]);
		  avgdy+=sin(theta[list[j]][0]);
		  np++; //# of swimmers inside the unit circle
		}
	    }
	  
	  rantheta=(ran2(&iseed)-0.5)*2*PI;
	  avgdx+=eta*np*cos(rantheta);
	  avgdy+=eta*np*sin(rantheta);
	  
	  
	  if(avgdy>0)
	    avgtheta=acos(avgdx/sqrt(pow(avgdx,2)+pow(avgdy,2)));
	  if(avgdy<0)
	    avgtheta=-acos(avgdx/sqrt(pow(avgdx,2)+pow(avgdy,2)));
	  if(avgdy==0)
	    {
	      if(avgdx>0)avgtheta=0;
	      if(avgdx<0)avgtheta=PI;
	    }
	  
	  theta[i][1]=avgtheta;	
	  //================================================================//
	}//End Loop Over Swimmers
      
      for(i=0;i<LS;i++)
	{
	  NBOX[i]=0;
	  for(j=0;j<N;j++)
	    BOX[i][j]=0;
	}
      
      for(i=0;i<N;i++) //Loop over swimmers for position update
	{
	  theta[i][0]=theta[i][1];
	  x[i]=x[i]+V0*cos(theta[i][1])*dt;
	  y[i]=y[i]+V0*sin(theta[i][1])*dt;
	  
	  if(x[i]>L)x[i]=x[i]-L; //PBC
	  if(x[i]<0)x[i]=L+x[i]; //PBC
	  if(y[i]>L)y[i]=y[i]-L; //PBC
	  if(y[i]<0)y[i]=L+y[i]; //PBC
	  
	  X=x[i];
	  Y=y[i];
	  box=Y*L+X;
	  
	  BOX[box][NBOX[box]]=i;
	  NBOX[box]++;
	}
      
      if(time%100==0)
	{
	  vx=0.0;
	  vy=0.0;
	  for(i=0;i<N;i++)
	    {
	      vx+=V0*cos(theta[i][0]);
	      vy+=V0*sin(theta[i][0]);
	    }
	  va=sqrt(pow(vx,2)+pow(vy,2))/(V0*N);
	  fprintf(fp,"%ld	%lf\n",time,va);
	  
	  if(time>SAT_TIME)
	    {
	      nsamp++;
	      PHI+=va;
	      PHI2+=pow(va,2);
	      PHI4+=pow(va,4);
	    }//If Saturation Time
	}
	if(time%10000==0){
		printf("time_%ld\n",time);
	}
	//If time%100
      
      //=============================================================//
      
 /*    if((time%1000)==0)
	{
	  frame++;
	  sprintf(FILE_SNAP,"Snap_eta%5.3f_time%ld.dat",ETA,time);
	  fs=fopen(FILE_SNAP,"w");
	  for(i=0;i<N;i++)
	    fprintf(fs,"%lf %lf %lf %lf\n",x[i],y[i],0.7*cos(theta[i][1]),0.7*sin(theta[i][1]));
	  fclose(fs);
	  
	 /*gp = popen(GNUPLOT,"w"); // 'gp' is the pipe descriptor 
	  fprintf(gp,"set term pngcairo enhanced color linewidth 1.5\n");
	  
	  if(frame<10)
	    fprintf(gp,"set output \"eta%5.3f_000%ld.png\" \n",ETA,frame);
	  if(frame>=10 && frame<100)
	    fprintf(gp,"set output \"eta%5.3f_00%ld.png\" \n",ETA,frame);
	  if(frame>=100 && frame<1000)
	    fprintf(gp,"set output \"eta%5.3f_0%ld.png\" \n",ETA,frame);
	  if(frame>=1000 && frame<10000)
	    fprintf(gp,"set output \"eta%5.3f_%ld.png\" \n",ETA,frame);
	  
	  fprintf(gp,"set size square\n");
	  fprintf(gp, "set xrange [%d:%d]\n",0,L);
	  fprintf(gp, "set yrange [%d:%d]\n",0,L);
	  fprintf(gp, "set title \"eta = %5.3f ; time=%ld\" \n",ETA,time);
	  fprintf(gp, "unset key\n");
	  fprintf(gp, "set tics scale 0\n");
	  fprintf(gp, "set ticslevel 0\n");
	  
	  fprintf(gp, "plot \"Snap_eta%5.3f.dat\" using 1:2:3:4 with vectors head filled lt 9 lw 0.25\n",ETA);
	  
	  fclose(gp);
	}*/
      
    }//loop over time
  
  PHI=PHI/nsamp;
  PHI2=PHI2/nsamp;
  PHI4=PHI4/nsamp;
  
  BC=1-PHI4/(3*pow(PHI2,2));//BINDER CUMULANT FORMULA SUBTITUTION
  CHI=PHI2-pow(PHI,2);
  
  fprintf(fp1,"%lf	%lf\n",ETA,PHI);
  fprintf(fp2,"%lf	%lf\n",ETA,BC);
  fprintf(fp3,"%lf	%lf\n",ETA,CHI);

 /*for(long i=0; i<TTIME; i=i+1000){ 
	gp = popen(GNUPLOT,"w"); // 'gp' is the pipe descriptor 
	  fprintf(gp,"set term pngcairo enhanced color linewidth 1.5\n");
	  fprintf(gp,"set size square\n");
	  fprintf(gp, "set xrange [%d:%d]\n",0,L);
	  fprintf(gp, "set yrange [%d:%d]\n",0,L);
	  fprintf(gp, "set title \"eta = %5.3f ; time=%ld\" \n",ETA,i);
	  fprintf(gp, "unset key\n");
	  fprintf(gp, "set tics scale 0\n");
	  fprintf(gp, "set ticslevel 0\n");
	  fprintf(gp,"set output \"eta%5.3f_%ld.png\"\n",ETA,i);
	  r
	  fprintf(gp, "plot \"Snap_eta%5.3f_time%ld.dat\" using 1:2:3:4 with vectors head filled lt 9 lw 0.25\n",ETA,i);
  }
  pclose(gp);*/
  
  return 0;      
}  //end main




//==== RANDOM NUMBER GENERATOR ====//

float ran2(long *idum)
{
  int j;
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  float temp;
  if (*idum <= 0) {
    if (-(*idum)<1) *idum = 1;
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) {
      k=(*idum)/IQ1;
      *idum = IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;
  *idum=IA1*(*idum-k*IQ1)-k*IR1;
  if (*idum < 0) *idum += IM1;
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV;
  iy=iv[j]-idum2;
  iv[j] = *idum;
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}