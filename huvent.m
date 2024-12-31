
function v=huvent(bits)
    % for k=0:t
    %   a=4096.^(t-k).*(2/(8*k + 1) + 1/(6*k + 1) + 2^-1/(12*k + 3) + 2^-8/(6*k + 5) - 2^-5/(6*k + 3) - 2^-4/(12*k + 7) - 2^-5/(8*k + 5) - 2^-8/(12*k + 11));
    %   a0=a0+a;
    % end 
    ce4096=[8,1,2;6,1,1;12*2^1,3*2^1,1;6*2^8,5*2^8,1;6*2^5,3*2^5,-1;12*2^4,7*2^4,-1;8*2^5,5*2^5,-1;12*2^8,11*2^8,-1];
% 计算16^t * sum(1/16^k * (4/(8k+1) - 2/(8k+4) - 1/(8k+5) - 1/(8k+6)))
    ce16=[8,1,4;8,4,-2;8,5,-1;8,6,-1];
%浮点运算，可以减小舍入误差，特别是huventbad(4)舍入结果为5A2,实际是5A3
%与calpi基础BBP算法可以相互印证
    ce30=[60,3,         +6;...
60*2,5*2,               +5;...
60,6,                   +3;...
60*2^2,9*2^2,           +3;...
60*2^3,10*2^3,          +5;...
60*2^6,15*2^6,          -1;...
60*2^6,18*2^6,          -3;...
60*2^8,21*2^8,          -3;...
60*2^11,25*2^11,        -5;...
60*2^11,27*2^11,        +3;...
60*2^13,30*2^13,        +1;...
60*2^14,33*2^14,        +3;...
60*2^16,35*2^16,        -5;...
60*2^17,39*2^17,        -3;...
60*2^18,42*2^18,        -3;...
60*2^21,45*2^21,        -1;...
60*2^23,50*2^23,        +5;...
60*2^23,51*2^23,        +3;...
60*2^24,54*2^24,        +3;...
60*2^26,55*2^26,        +5;...
60*2^26,57*2^26,       +3];
b=[       32           1           4
          32           4          -2
          32           5          -1
          32           6          -1
         128          36           1
         256          96          -1
         512         208          -1
         512         224          -1
        2048        1088           1
        4096        2560          -1
        8192        5376          -1
        8192        5632          -1
       32768       25600           1
       65536       57344          -1
      131072      118784          -1
      131072      122880          -1];
    % p6=hpi(bits,b,4);p6(:)'
     % p30=hpi(bits,ce30,30);
     % p7=bin2hex(p30(:)')
    % p3=hpi(bits,ce4096,12);p3(:)'
    numprocs=1;
    calpi(bits,ce4096,12,numprocs)
   p0=calpi(bits,ce16,4,numprocs)
   p=hpi(bits,ce16,4)
   % P = agm_pi(100000);
   N1 = 13; M1 = 1; b1 = [8 1 4; 8 4 -2; 8 5 -1; 8 6 -1];
N2 = 7; M2 = 2; b2 = [16 1 4; 16 4 -2; 16 5 -1; 16 6 -1; 64 36 1; 128 96 -1; 4096 3328 -1; 4096 3584 -1]; % Adjusted b2 for MATLAB indexing
N3 = 5; M3 = 3; b3 = [24 2 8; 24 3 6; 24 4 -8; 24 6 -3; 24 8 -3; 48 20 1; 192 96 -1; 768 15*32 -3; 384 256 -3; 384 288 -1; 768 400 -1]; % Adjusted b3 for MATLAB indexing and division
% Note: b4 is already in a format suitable for MATLAB
N4 = 3; M4 = 4; b4 = [32 1 4; 32 4 -2; 32 5 -1; 32 6 -1; 128 36 1; 256 96 -1; 512 208 -1; 512 224 -1; 2048 1088 1; 4096 2560 -1; 8192 5376 -1; 8192 5632 -1; 32768 25600 1; 65536 57344 -1; 131072 118784 -1; 131072 122880 -1];
 
bits = 100; % Example bits value
r = bbp_matlab(bits, M1, N1, b1);
disp(r);
end
function P = agm_pi(d) 
 % AGM _ PI Arithmetic-geometric mean for pi. 
 % agm _ pi(d) produces d decimal digits. 
 digits(d) 
 a = vpa(1,d); 
 b = 1/sqrt(vpa(2,d));  
 s = 1/vpa(4,d); 
 p = 1; 
 n = ceil(log2(d)); 
 for k = 1:n 
     c = (a+b)/2; 
     b = sqrt(a*b); 
     s = s - p*(c-a)^2; 
     p = 2*p; a = c; 
 end 
 P = a^2/s; 
end
function h=bin2hex(b)
   k=fix(length(b)/4);
   for i=1:k
       j=(i-1)*4+1:i*4;
       h(i)=dec2hex(bin2dec(b(j)));
   end
end
 function v=hpi(bits,ce,bytes)
 buf=ones(2,1);buf(2)=bitshift(buf(1), bytes);
 for t=1:bits
   a0=0; c = 0; % 用于补偿舍入误差的变量
   for k=0:t
       for i=1:length(ce)
         ak = kahanSum(t,k, ce(i,1),ce(i,2),buf);
         if mod(bytes,4)==0
             y=ak*ce(i,3);
         else
             y=ak*ce(i,3)*(-1)^k;
         end
         a0 = mod(a0+y,buf(end)) ;      
         if a0<0,a0=a0+buf(end);end         
       end
      
   end 
 
     p=convert(a0,bytes);
     v(:,t)=p;
 end 
  
end
 
function p=convert(a0,bytes)
  if mod(bytes,4)==0
     result=dec2hex(fix(a0));  
     p=repmat('0', bytes/4, 1);
     for i=1:bytes/4+1
        if length(result) <i
           p(bytes/4-i+2:bytes/4)=result(1:i-1);
           break;
        end
     end
  elseif mod(bytes,2)==0
     result=dec2bin(fix(a0));  
     p=repmat('0', bytes, 1);
     for i=1:bytes+1
        if length(result) <i
           p(bytes-i+2:bytes)=result(1:i-1);
           break;
        end
     end
  end
 end
 function y = kahanSum(t, k,u,m,buf)
   
        a = (k * u) + m;
    % y = 4096^(t-k)/a ; 
    p=t-k;
    result = 1;
    base = buf(2);
    
    % 将 p 转换为二进制并遍历每一位
    while p > 0
        % 如果 p 的当前位为 1，则将 result 乘以 base 并取模 a*4096
        b=a * buf(2);
        if mod(p, 2) == 1
            result = mod(result * base, b);
        end
        base = mod(base * base, b);
        % 右移 p 的二进制表示
        p = bitshift(p,-1);
    end
    y=result/a;
    % end
 end   
 
function p=calpi(bits,ce,M,numprocs)
N=fix(63/M);%N<14会出现很多位数错误，但又不能大于14，否则溢出
newbuf=ones(N,1);for i=1:N+1,newbuf(i)=bitshift(1,M*(i-1));end
if mod(M,4)==0,p='3.';else p='11.';end
for t=1:bits
    pid=zeros(numprocs,1);
    for i=1:numprocs
       for j=length(ce)
         s = series(i-1, t, numprocs,ce(j,:),newbuf); 
         pid(i) = mod(pid(i)+s,newbuf(end));
       end
    end
    totalSum=0;
    for i=1:numprocs
      totalSum=totalSum+pid(i);
      totalSum = mod(totalSum,  newbuf(end));
    end
    result = convert(bitshift(totalSum,-(N-1)*M),M);
    p=[p,result'];
 end

end


function s = series( from, ic, inc,ce,newbuf)
    s = 0; m=log2(newbuf(end));
    for k = from:inc:ic+length(newbuf)-1
        ak = (k * ce(1)) + ce(2); % 计算分母
        if ic>=k
          p = ic - k; 
          t = expm(p, ak,newbuf); %计算(16^(ic-k) % ak)
          
          x = t;
          y = 0;
          for n=1:m  %放大2^n倍，如果是ffffffffffffff可能结果有错误，因此必然存在错误，不知道什么时候在哪里，除非N>48
            y=y+fix(x/ak); % 整数除法
            x=mod(x, ak); % 模运算
            x=x+x;
            y=mod(y+y,newbuf(end));
          end
            % y=fix(t*newbuf(end)/ak);
        else
            y=fix(newbuf(length(newbuf)-(k-ic))/ak);
        end
        s = mod(s + y*ce(3)+newbuf(end), newbuf(end));
    end
end

function r = expm(p, ak,newbuf)

    if ak == 1
        r = 0; 
        return;
    end
    if p < length(newbuf)
        r = mod(newbuf(p+1), ak);
        return;
    end
    n=log2(newbuf(2));
    p1 = p;
    i = 0;
    while p1 > 0
        buf(i+1) = mod(p1,newbuf(2)); 
        p1 = Floor(p1,newbuf(2)); 
        i = i + 1;
    end
    
    r=1;
    for j=i:-1:2
        for m = 1:buf(j)*n
            r = mod(r + r, ak);
        end
        for k = 1:n
            p1 = r;
            t = r;
            r = 0;
            while p1>0
               r = mod(r + mod(t * mod(p1,newbuf(2)), ak), ak); 
               t = mod(t * newbuf(2), ak);                     
               p1 = Floor(p1, newbuf(2));                     
            end
        end
    end
    for m = 1:buf(1)*n
        r = mod(r + r, ak);
    end
   
end

function count = Floor(xx, ak)
    count = 0;
    remaining = xx;
 
    while remaining >= ak
        remaining = remaining - ak;
        count = count + 1;
    end
end
 function r = bbp_matlab(bits, M, N, b)
    numprocs = 6;

    % Helper function to compute LCM of two numbers
    function lcm = lcm_matlab(a, b)
        lcm = abs(a * b) / gcd(a, b);
    end

    % Helper function to compute the modular exponentiation
    function r = expm_matlab(p, a, base)
        r = 1;
        s = base;
        while p > 0
            if mod(p, 2) == 1
                r = mod(r * s, a);
            end
            s = mod(s * s, a);
            p = bitshift(p, -1);
        end
    end

    % Helper function to process a chunk of data in parallel-like fashion
    function result = process_chunk_matlab(args)
        k = args{1};
        t = args{2};
        M = args{3};
        N = args{4};
        b = args{5};

        MAX = 2^(N * M * 4);
        result = 0;
        for k_idx = k:(t + N - 1)
            k_local = k_idx; % Adjust k_local for each iteration
            p = t - k_local;
            for i = 1:size(b, 1)
                ak = (k_local * b(i, 1)) + b(i, 2);
                if p >= 0
                    x = expm_matlab(p, ak, 2^(M * 4));y=0;
                    for n=1:M*N*4  %放大2^n倍，如果是ffffffffffffff可能结果有错误，因此必然存在错误，不知道什么时候在哪里，除非N>48
                        y=y+fix(x/ak); % 整数除法
                        x=mod(x, ak); % 模运算
                        x=x+x;
                        y=mod(y+y,MAX);
                    end
                    result = mod(result + MAX + y, MAX);
                else
                    result = mod(result + MAX + bitshift(b(i, 3), (N + p) * M * 4) / ak, MAX);
                end
            end
        end
    end

    % Main parallel processing function (simulated using a loop in MATLAB)
    function totalSum = parallel_process_matlab(t, M, N, b)
        MAX = 2^(N * M * 4);
        results = zeros(1, numprocs);
        for k = 1:numprocs
            args = {k, t, M, N, b};
            results(k) = process_chunk_matlab(args);
        end
        totalSum = mod(sum(results), MAX);
    end

    % Main BBP function
    r = '';
    for t = floor((bits - 1) / M):floor((bits - 1) / M)
        totalSum = parallel_process_matlab(t, M, N, b);
        result = dec2hex(bitshift(totalSum, -(N - 1) * M * 4), M);
        while length(result) < M
            result = ['0' result];
        end
        r = [r result];
    end
end