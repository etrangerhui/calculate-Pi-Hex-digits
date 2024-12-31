

import multiprocessing
from time import time
import math
from functools import reduce

N1 = 13;M1 = 1;b1=[[8,1,4],[8,4,-2],[8,5,-1],[8,6,-1]];
N2 = 7; M2 = 2;b2=[[16,1,4],[16,4,-2],[16,5,-1],[16,6,-1],[16*4,9*4,1],[16*8,12*8,-1],[16*16,13*16,-1],[16*16,14*16,-1]];
N3 = 5; M3 = 3;b3=[[24,2,256//32],[24,3,192//32],[24,4,-256//32],[24,6,-96//32],[24,8,-96//32],[24*2,10*2,1],[24*8,12*8,-1],[24*32,15*32,-3],[24*16,16*16,-3],[24*16,18*16,-1],[24*32,20*32,-1]];
#b3=[[8,1,2],[6,1,1],[24,6,1],[1536,1280,1],[192,96,-1],[192,112,-1],[256,160,-1],[3072,2816,-1]];
#为了补充huvent公式的舍入，从t项增加为t+1项求和，求余也从a变为a*base
N4 = 3; M4 = 4;b4=[[32,1,4],[32,4,-2],[32,5,-1],[32,6,-1],[128,36,1],[256,96,-1],[512,208,-1],[512,224,-1],[2048,1088,1],[4096,2560,-1],[8192,5376,-1],[8192,5632,-1],[32768,25600,1],[65536,57344,-1],[131072,118784,-1],[131072,122880,-1]]

numprocs = 6

 
def lcm_multiple(numbers):
    def lcm(a, b):
        return abs(a * b) // math.gcd(a, b)
    return reduce(lcm, numbers)

def process_chunk(args):
    k, t ,M,N,b= args
    MAX=1<<(N*M*4)  
    result=0
    for k in range(k, t+N, numprocs):
        p = t - k
        for i in range(len(b)):
            ak = (k * b[i][0]) + b[i][1]  
            if p>=0:
                s=expm(p,ak,1<<(M*4))                 
                result+=MAX  +  ((s*b[i][2])<<(N*M*4))//ak  
            else:
                result+= MAX  +((b[i][2])<<((N+p)*M*4))//ak   
            result%=MAX
    return result
 

def parallel_process(t,M,N,b):
    MAX=1<<(N*M*4)  
    with multiprocessing.Pool(processes=numprocs) as pool:
        tasks = [(k,t,M,N,b) for k in range(numprocs)]          
        results = pool.map(process_chunk, tasks)
        
    return sum(results)%MAX
def expm(p,a,base):
    r=1
    s=base
    #b=a*base
    while p>0:
        if(p%2):r=r*s%a;
        s=s*s%a;
        p>>=1;
    return r
def bbp(bits,M,N,b):
    r = ''
   
    for t in range(((bits-1)//M), ((bits-1)//M)+1): 
        totalSum=parallel_process(t,M,N,b)  
        result = hex(totalSum>>(N-1)*M*4)[2:]  # 去掉 '0x' 前缀
        while len(result)<M:result='0'+result            
        r+=result;
    
    return r
def find_unique_end_string_index(strings):
    
    end_char_indices = {}
    for index, string in enumerate(strings):
        end_char = string[-1]
        if end_char in end_char_indices:
            end_char_indices[end_char].append(index)
        else:            
            end_char_indices[end_char] = [index]
    for end_char, indices in end_char_indices.items():
        if len(indices) == 1:
            return indices[0]    
    return -1



if __name__ == '__main__':
    
    M=[M1,M2,M3,M4];N=[N1,N2,N3,N4];b=[b1,b2,b3,b4];
    name=["BBP","jason","huvent","Bellard"]
    # M=[M1,M2,M4];N=[N1,N2,N4];b=[b1,b2,b4];
    # name=["BBP","jason","Bellard"]
    n=50
    a=lcm_multiple(M)
    for bits in range(a,n*a+1,a):
        result=[]
        for i in range(len(b)):
            start =time()            
            result.append(bbp(bits,M[i],N[i],b[i]))
        index = find_unique_end_string_index(result)
        if index != -1:        
            print(f"The {name[index]} is: {result[index]}")
            print(bits,index,result)
        
            
        