#include "mpreal.h"
#include<iostream>
#include <stdio.h>
#include<iomanip>
using mpfr::mpreal;
using namespace std;

#include <gmp.h>  // GMP 的 C 接口头文件（可选，但通常包含以获取 GMP 定义的常量）
#include <gmpxx.h> // GMP 的 C++ 接口头文件
#if defined(THREE)
  #define N 12
  #define M 5
  //static int b[][3]={{8,1,2},{6,1,1},{24,6,1},{1536,1280,1},{192,96,-1},{192,112,-1},{256,160,-1},{3072,2816,-1}};
  static long long b[][3]=
{{24   ,     2 ,256/32 },
 {24   ,     3 ,192/32 },
 {24   ,     4 ,-256/32},
 {24   ,     6 ,-96/32 },
 {24   ,     8 ,-96/32 },
 {24*2 ,  10*2 ,1      },
 {24*8 ,  12*8 ,-1     },
 {24*32, 15*32 ,-3     },
 {24*16, 16*16 ,-3     },
 {24*16, 18*16 ,-1     },
 {24*32, 20*32 ,-1     }};
#elif defined(TWO)
  #define N 8
  #define M 7
  static long long b[][3]=
 {{16   ,     1 ,4 },
 {16   ,     4 ,-2 },
 {16   ,     5 ,-1},
 {16   ,     6 ,-1 },
 {16*4   , 9*4 ,1 },
 {16*8 ,  12*8 ,-1 },
 {16*16 ,13*16 ,-1 },
 {16*16, 14*16 ,-1 }}; 

#elif defined(FOUR)
  #define N 16
  #define M 3 

  static long long b[][3]=
         {{      32 ,          1    ,       4},
          {      32 ,          4    ,      -2},
          {      32 ,          5    ,      -1},
          {      32 ,          6    ,      -1},
          {     128 ,         36    ,       1},
          {     256 ,         96    ,      -1},
          {     512 ,        208    ,      -1},
          {     512 ,        224    ,      -1},
          {    2048 ,       1088    ,       1},
          {    4096 ,       2560    ,      -1},
          {    8192 ,       5376    ,      -1},
          {    8192 ,       5632    ,      -1},
          {   32768 ,      25600    ,       1},
          {   65536 ,      57344    ,      -1},
          {  131072 ,     118784    ,      -1},
          {  131072 ,     122880    ,      -1}};
#elif defined(FIVE)
  #define N 30
  #define M 2
  static long long b[][3]=
 {{ 0000000060. ,  0000000003.  , 6 },
  { 0000000120. ,  0000000010.  , 5 }, 
  { 0000000060. ,  0000000006.  , 3 }, 
  { 0000000240. ,  0000000036.  , 3 }, 
  { 0000000480. ,  0000000080.  , 5 }, 
  { 0000003840. ,  0000000960.  ,-1 }, 
  { 0000003840. ,  0000001152.  ,-3 }, 
  { 0000015360. ,  0000005376.  ,-3 }, 
  { 0000122880. ,  0000051200.  ,-5 }, 
  { 0000122880. ,  0000055296.  , 3 }, 
  { 0000491520. ,  0000245760.  , 1 }, 
  { 0000983040. ,  0000540672.  , 3 }, 
  { 0003932160. ,  0002293760.  ,-5 }, 
  { 0007864320. ,  0005111808.  ,-3 }, 
  { 0015728640. ,  0011010048.  ,-3 }, 
  { 0125829120. ,  0094371840.  ,-1 }, 
  { 0503316480. ,  0419430400.  , 5 }, 
  { 0503316480. ,  0427819008.  , 3 }, 
  { 1006632960. ,  0905969664.  , 3 }, 
  { 4026531840. ,  3690987520.  , 5 }, 
  { 4026531840. ,  3825205248.  , 3 }};
#else
  #define N 4
  #define M 13     
  static long long b[][3]={{8,1,4},{8,4,-2},{8,5,-1},{8,6,-1}};
#endif

          
const long long MAX=long long(1)<<(M*N);
const  long base=long long(1)<<N;
void printBinary(unsigned long long n, int width) 
{
    for (int i = width - 1; i >= 0; i--) {
        putchar((n & (long long(1) << i)) ? '1' : '0');
    }
    
}

long long expm(long long p,long long a,long long base)  //(base^p)%a
{
      long long r=1,c,d,q;
      if(a==1)r=0;
      else if(p==0)r=1;
      else if(p<=M)r=(long long(1)<<(p*N))%a;
      else 
      {
        long long bits[100]={0},j=0; 
        for(j=0;p>0;j++,p>>=N)bits[j]=p%base;         		 	
        for(d=j-1;d>0;d--)
        {
         for(c=0;c<bits[d]*N;c++)r=(r+r)%a;
         for(int e=0;e<N;e++)
         {
           p=r;q=r;r=0;
           while(p>0)
           {
             r=(r+(q*(p%base))%a)%a;
             q=(q*base)%a;
             p>>=N;
            }
          }         		 		
         }
         for(c=0;c<bits[0]*N;c++)r=(r+r)%a;
       }
    return r;
  }
long long  addl(long long t,long long *b,int m,long long base)
{
          long long y=0;
#pragma omp parallel for reduction(+:y)  schedule(dynamic)
         	for(long long k=0;k<=t+M-1;k++)
         	{
         		 for(int i=0;i<m;i++)
         		 {
         		 	   long long p=t-k,a=(k)*b[3*i]+b[3*i+1],z=0;
         		 	 //printf("%lld %lld %lld\n",t-k,a,r);
         		 	 if(k<=t)
         		 	 	{
         		 	   long long x=expm(p,a,base);
         		 	   for(int c=0;c<M*N;c++){z+=x/a;x%=a;x+=x;z=(z+z)%MAX;} 
         		 	 }
         		 	 else z=(1<<(M-(k-t))*N)/a;
         		 	   if((N%4) && (k%2))z=-z;        		 	 
         		 	   y+=z*b[3*i+2];
         		 	   //printf("%lld %lld %lld %lld %lld\n",base,k,p,x,y);
         		   }  
         		   //if(y<0)printf("%lld\n",y);
         		}
         		y%=MAX;if(y<0)y+=MAX;
         		
         		return y;
  }
 
int main(int argc,char*argv[])
{
	if(argc<3)
		{
			printf("Usage:%s <digits> <flag>\n",argv[0]);
			printf("flag:1-digits/3;2-digits;3-agm(log2;4-log10;5-chudnovsky;\n");
			printf("flag:6 hex*[1234]/bin*30[5] #if compiled with define ONE TWO THREE FOUR FIVE \n");
			return -1;
		}
	long long  i,k,t,n=atoi(argv[1]),flag=atoi(argv[2]);
  mpreal::set_default_prec(mpfr::digits2bits(n+10));
  	if(flag==1)
  {
  mpreal v(1),u=v;
  for(i=n;i>0;i--)u=1+u/(2+v/i);                    //slow 3 times
  u*=2;cout<<setprecision(n/3.3)<<u;
} 
if(flag==2)
	{
   mpreal u(8+5*n);                                //same
  for(i=n;i>0;i--)
  {
  	mpreal p=i;
  	mpreal d=(p+1)*(2*p+1)/3/(3*p+4)/(3*p+5);
  	u=8+5*(p-1)+d*u;
  }
  u=3+u/60;cout<<setprecision(n)<<u;
}
  if(flag==3)
  {
  mpreal a(1),b=1/sqrt(mpreal(2)),s=1/mpreal(4),p(1);  //agm
  
  mpreal m=log2(mpreal(n));
  for(k=1;k<m+1;k++)
  {
  	mpreal c=(a+b)/2;
  	b=sqrt(a*b);
  	s-=p*(c-a)*(c-a);
  	p=2*p;
  	a=c;
  }cout<<setprecision(n)<<a*a/s;
} 
if(flag==4)
	{
  mpreal a=1/mpreal(3),s=sqrt(mpreal(2))-1;  //Order 16
  for(k=1;k<=long long(log10(n))+1;k++)
  {
  	mpreal ss=pow(1-pow(s,4),0.25);
  	mpreal t=1+ss,m1=pow((1+s)/t,4),m2=pow(t,-4),u=pow(8*ss*(1+ss*ss),0.25);
  	s=pow(1-ss,4)/((t+u)*(t+u)*(t*t+u*u));
  	a=16*m1*a+pow(mpreal(4),2*k-1)/3*(1-12*m2-4*m1);  	
  }
  cout<<setprecision(n)<<1/a;
 } 
 if(flag==5)
 	{
  mpz_class one(1);
  for(i=0;i<n+10;i++)one=one*10;  
 	mpz_class k(1), a_k(one), a_sum(one), b_sum(0), C(640320);
 	mpz_class C3_OVER_24(C * C * C / 24), total, pi, sqrt_result, factor;
 
    while (true) {
        a_k *= -(6 * k - 5) * (2 * k - 1) * (6 * k - 1);
        a_k /= (k * k * k * C3_OVER_24);
        a_sum += a_k;
        b_sum += k * a_k;
        ++k;
 
        if (a_k == 0) {
            break;
        }
    }
 
     total = 13591409 * a_sum + 545140134 * b_sum;
     sqrt_result = sqrt(mpz_class(10005) * one * one); 
     pi = (426880 * sqrt_result * one) / total;
     cout<<"3."<<pi.get_str().substr(1,n-1);
   }
   if(flag==6)
   	{
   		  #ifdef FIVE
   			   printf("11.");
   			#else
   			   printf("3.");
   			#endif
   			 
     for(t=1;t<=n;t++)
      {
         	#if defined(THREE)
         	   long long y=addl(t-1,&b[0][0],11,base);         	   
         	   printf("%03x",(y>>(M-1)*N));
         	#elif defined(FOUR)
         	   long long y=addl(t-1,&b[0][0],16,base);
         	   printf("%04x",(y>>(M-1)*N));
         #elif defined(TWO)
         	   long long y=addl(t-1,&b[0][0],8,base);
         	   printf("%02x",(y>>(M-1)*N));   
         	#elif defined(FIVE)
         	   long long y=addl(t-1,&b[0][0],21,base);
         	   //printf("%07x",(y>>(M-1)*N+2)); 
         	   printBinary((y>>(M-1)*N),N);
         	#else         	
         	   long long y=addl(t-1,&b[0][0],4,base);         	   
         	   printf("%x", (y>>(M-1)*N));
         	  
         	#endif
      }
  }
 
 	
}