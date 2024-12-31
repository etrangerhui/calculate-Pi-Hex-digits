
#include<iostream>
#include <stdio.h>
 #define N3 12 
 static long long b3[][3]={{8,1,2},{6,1,1},{24,6,1},{1536,1280,1},{192,96,-1},{192,112,-1},{256,160,-1},{3072,2816,-1}};
 //static long long b3[][3]={{24,2,256/32},{24,3,192/32},{24,4,-256/32},{24,6,-96/32},{24,8,-96/32},{24*2,10*2,1},{24*8,12*8,-1},{24*32,15*32,-3},{24*16,16*16,-3},{24*16,18*16,-1},{24*32,20*32,-1}};
 #define N2 8 
 static long long b2[][3]={{16,1,4},{16,4,-2},{16,5,-1},{16,6,-1},{16*4,9*4,1},{16*8,12*8,-1},{16*16,13*16,-1},{16*16,14*16,-1}}; 
 #define N4 16
 static long long b4[][3]= {{32,1,4}, {32,4,-2}, {32,5,-1}, {32,6,-1}, {128,36,1}, {256,96,-1}, {512,208,-1}, {512,224,-1}, {2048,1088,1}, {4096,2560,-1}, {8192,5376,-1}, {8192,5632,-1}, {32768,25600,1}, {65536,57344,-1}, {131072,118784,-1}, {131072,122880,-1}};
 #define N5 30
 static long long b5[][3]= {{60.,3.,6},{120.,10.,5}, {60.,6.,3}, {240.,36.,3}, {480.,80.,5}, {3840.,960.,-1}, {3840.,1152.,-3}, {15360.,5376.,-3}, {122880.,51200.,-5}, {122880.,55296.,3}, {491520.,245760.,1}, {983040.,540672.,3}, {3932160.,2293760.,-5}, {7864320.,5111808.,-3}, {15728640.,11010048.,-3}, {125829120.,94371840.,-1}, {503316480.,419430400.,5}, {503316480.,427819008.,3}, {1006632960.,905969664.,3}, {4026531840.,3690987520.,5}, {4026531840.,3825205248.,3}};
 #define N1 4
 static long long b1[][3]={{8,1,4},{8,4,-2},{8,5,-1},{8,6,-1}};

void print(unsigned long long n,int width,int bits,long long t) 
{
 long long c=(long long(1) << bits);
 int m=width%bits;
 if(t%2==0)width-=m;
 for (int i = width-bits; i >= 0; i-=bits) 
 {
 	 long long x=(n>>i)%c;
 	 if(x>=0&&x<=9)putchar(x+'0');else putchar(x-10+'a');
 }
}

long long expm(long long p,long long a,long long base,int M,int N) //(base^p)%a //最大出现a*base
{
 long long r=1,c,d,q;
 if(a==1)r=0;
 else if(p==0)r=1;
 else if(p<=M)r=(long long(1)<<(p*N))%a;
 else 
  {
   long long bits[100]={},j=0; 
   for(j=0;p>0;j++,p>>=N)bits[j]=p%base; 		 	
   for(d=j-1;d>0;d--)
   {
    for(c=0;c<bits[d]*N;c++)r=(r+r)%a;
    for(int e=0;e<N;e++)
    {
     p=r;q=r;r=0;
     while(p>0){r=(r+(q*(p%base))%a)%a;q=(q*base)%a;p>>=N;}
    } 		 		
   }
   for(c=0;c<bits[0]*N;c++)r=(r+r)%a;
  }
 return r;
}
long long  addl(long long t,long long *b,int m,long long base,int M,int N,long long MAX)
{
          long long y=0;
#pragma omp parallel for reduction(+:y)  schedule(dynamic)
         	for(long long k=0;k<=t+M-1;k++)
         	{
         		 for(int i=0;i<m;i++)
         		 {
         		 	   long long p=t-k,a=(k)*b[3*i]+b[3*i+1];
 		 	           long long x,z=0;
 		 	           if(k<=t)x=expm(p,a,base,M,N);else x=1;
         		 	   if(k<=t)for(int c=0;c<M*N;c++){z+=x/a;x%=a;x+=x;z=(z+z)%MAX;} 
         		 	   //else for(int c=0;c<(M-(k-t))*N;c++){z+=x/a;x%=a;x+=x;z=(z+z)%MAX;}
         		 	   else z=(x<<(M-(k-t))*N)/a;
         		 	   if((N%4) && (k%2))z=-z;        		 	 
         		 	   y+=z*b[3*i+2];
         		   }  
         		}
         		y%=MAX;if(y<0)y+=MAX;
         		
         		return y;
  }

int main(int argc,char*argv[])
{
	if(argc<2)
		{
			printf("long long max: %lld\n",LLONG_MAX);
			printf("Usage:%s <digits>\n",argv[0]);
			return -1;
		}
	int N[]={N1,N2,N3,N4,N5};
	long long (*b[])[3]={b1,b2,b3,b4,b5};
	int m[]={sizeof(b1)/sizeof(b1[0]),sizeof(b2)/sizeof(b2[0]),sizeof(b3)/sizeof(b3[0]),sizeof(b4)/sizeof(b4[0]),sizeof(b5)/sizeof(b5[0])};
	long long j,n;
	sscanf(argv[1],"%lld",&n);j=n;
	for(int i=0;i<5;i++)
	{	
   const int M=int((sizeof(long long)*8-1)/N[i]); 
   const long long MAX=long long(1)<<(M*N[i]);
   const long base=long long(1)<<N[i];
	 printf("\n");
   //for( j=n;j<=n;j++)
    {
     long long t,c=N[i]/4;
	   if(N[i]%4)t=(j-1)*4.0/N[i]+1;else t=(j-1)/c+1;;
 	   long long y=addl(t-1,&b[i][0][0],m[i],base,M,N[i],MAX);
 	   long long x=y>>((M-1)*N[i]); 
 	   if(N[i]%4)print(x,N[i],4,t);
     else print(x,N[i],4,t);
    }
  }
}
 
 	
