
function v=huvent(bits)
    % for k=0:t
    %   a=4096.^(t-k).*(2/(8*k + 1) + 1/(6*k + 1) + 2^-1/(12*k + 3) + 2^-8/(6*k + 5) - 2^-5/(6*k + 3) - 2^-4/(12*k + 7) - 2^-5/(8*k + 5) - 2^-8/(12*k + 11));
    %   a0=a0+a;
    % end 
M3= 12; b3=[8,1,2;6,1,1;12*2^1,3*2^1,1;6*2^8,5*2^8,1;6*2^5,3*2^5,-1;12*2^4,7*2^4,-1;8*2^5,5*2^5,-1;12*2^8,11*2^8,-1];
%b3 = [24 2 8; 24 3 6; 24 4 -8; 24 6 -3; 24 8 -3; 48 20 1; 192 96 -1; 768 480 -3; 384 256 -3; 384 288 -1; 768 640 -1];
M1 = 4; b1 = [8 1 4; 8 4 -2; 8 5 -1; 8 6 -1];
M5 = 30;b5=[60,3, +6;60*2,5*2, +5;60,6, +3;60*2^2,9*2^2, +3;60*2^3,10*2^3,+5;60*2^6,15*2^6,-1;60*2^6,18*2^6,-3;60*2^8,21*2^8,-3;60*2^11,25*2^11,-5;60*2^11,27*2^11,+3;60*2^13,30*2^13,+1;60*2^14,33*2^14,+3;60*2^16,35*2^16,-5;60*2^17,39*2^17,-3;60*2^18,42*2^18,-3;60*2^21,45*2^21,-1;60*2^23,50*2^23,+5;60*2^23,51*2^23,+3;60*2^24,54*2^24,+3;60*2^26,55*2^26,+5;60*2^26,57*2^26, +3]; 
 M2 = 8; b2 = [16 1 4; 16 4 -2; 16 5 -1; 16 6 -1; 64 36 1; 128 96 -1; 256 208 -1; 256 224 -1]; 
M4 = 16; b4 = [32 1 4; 32 4 -2; 32 5 -1; 32 6 -1; 128 36 1; 256 96 -1; 512 208 -1; 512 224 -1; 2048 1088 1; 4096 2560 -1; 8192 5376 -1; 8192 5632 -1; 32768 25600 1; 65536 57344 -1; 131072 118784 -1; 131072 122880 -1];
M6 = 10; b6=[8 2 -1;256 192 -1;10 1 4;10 3 -1;160 80 -1;160 112 -1;640 576 1] ;
[bbp(bits, M1, b1) ' ' bbp(bits, M2,  b2) ' ' bbp(bits, M3, b3) ' ' bbp(bits, M4, b4) ' ' bbp(bits, M5, b5) ' ' bbp(bits, M6, b6)]
end
function r = expm(p, ak,n)
    base=bitshift(1,n);
    if ak==1,r=0;return;end
    if p<2,r=mod(base^p, ak);return;end
    p1 = p;
    i = 0;
    while p1 > 0
        buf(i+1) = mod(p1,base); 
        p1 = bitshift(p1,-n); 
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
               r = mod(r + mod(t * mod(p1,base), ak), ak);  %最大出现ak*base
               t = mod(t * base, ak);                     
               p1 = bitshift(p1,-n);                     
            end
        end
    end
    for m = 1:buf(1)*n
        r = mod(r + r, ak);
    end
   
end
% function r = expm(p, ak, n)

%         if a==1,r=0;return;end
 %        base=bitshift(1,n);
%         if p<2,r=mod(base^p, ak);return;end
%         s = base;r = 1;
%         while p > 0
%             if mod(p, 2) == 1
%                 r = mod(r * s, ak);
%             end
%             s = mod(s * s, ak);  %最大出现ak*ak
%             p = bitshift(p, -1);
%         end
% end
 function result = process(t, M, N, b)

        MAX = 2^(N * M );
        result = 0;
        for k = 0:(t + N - 1)
            p = t - k;
            parfor i = 1:size(b, 1)
                ak =k * b(i, 1) + b(i, 2);
                y(i)=0;
                if p >= 0,x = expm(p, ak, M );m=M*N ;else, x=1;m=(N + p) * M  ;end
                for n=1:m  
                        y(i)=y(i)+fix(x/ak); 
                        x=mod(x,ak); 
                        x=bitshift(x,1);
                        y(i)=mod(bitshift(y(i),1),MAX);
                end

            end
            for i = 1:size(b, 1)
                if mod(M,4)>0 && mod(k,2)==1,s=-1;else s=1;end
                result = mod(result+MAX+s*y(i)*b(i, 3),MAX);
            end
        end
    end
 function h=bin2hex(b)
   k=fix(length(b)/4);
   for i=1:k
       j=(i-1)*4+1:i*4;
       h(i)=dec2hex(bin2dec(b(j)));
   end
end
 function r = bbp(bits, M, b)
    N=fix(63/M);
    m=fix((bits - 1)*4 / M);
    if m<0;beg=0;m=-m;else beg=m;end
    if beg==0&&mod(M,4)>0, r='';else r='';end
    for t =beg:m
        totalSum = process(t, M, N, b);
        if mod(t,2)==0&&mod(M,4)>0,k=mod(M,4);else k=0;end
        if beg==0&&mod(M,4)>0
            result=dec2bin(bitshift(totalSum, -(N - 1) * M));n=M;
        else
           n=fix(M/4);
           result = dec2hex(bitshift(totalSum, -(N - 1) * M - k ),n);
        end
        
        while length(result) < n
            result = ['0' result];
        end
         r = [r result];
        
    end
    if beg==0&&mod(M,4)>0,r=bin2hex(r);end
 end


%  function P = agm_pi(d) 
%  % AGM _ PI Arithmetic-geometric mean for pi. 
%  % agm _ pi(d) produces d decimal digits. 
%  digits(d) 
%  a = vpa(1,d); 
%  b = 1/sqrt(vpa(2,d));  
%  s = 1/vpa(4,d); 
%  p = 1; 
%  n = ceil(log2(d)); 
%  for k = 1:n 
%      c = (a+b)/2; 
%      b = sqrt(a*b); 
%      s = s - p*(c-a)^2; 
%      p = 2*p; a = c; 
%  end 
%  P = a^2/s; 
% end
% 
% %      p30=hpi(bits,ce30,30);
% %      p7=bin2hex(p30(:)')
%     % p3=hpi(bits,b3,12);p3(:)'
% %    p=hpi(bits,b1,4)
%  % p6=hpi(fix((bits - 1)*4 / M6)+1,b6,10);
% % p10=bin2hex(p6(:)')
% 
%  % p=hpi(bits,b1,4)
%    % P = agm_pi(100000);
%  function v=hpi(bits,ce,bytes)
%  buf=ones(2,1);buf(2)=bitshift(buf(1), bytes);
%  for t=1:bits
%    a0=0; 
%    for k=0:t
%        for i=1:length(ce)
%          ak = kahanSum(t,k, ce(i,1),ce(i,2),buf);
%          if mod(bytes,4)==0
%              y=ak*ce(i,3);
%          else
%              y=ak*ce(i,3)*(-1)^k;
%          end
%          a0 = mod(a0+y,buf(end)) ;      
%          if a0<0,a0=a0+buf(end);end         
%        end
%    end 
% 
%      p=convert(a0,bytes);
%      v(:,t)=p;
%  end 
% 
% end
% 
% function p=convert(a0,bytes)
%   if mod(bytes,4)==0
%      result=dec2hex(fix(a0));  
%      p=repmat('0', bytes/4, 1);
%      for i=1:bytes/4+1
%         if length(result) <i
%            p(bytes/4-i+2:bytes/4)=result(1:i-1);
%            break;
%         end
%      end
%   elseif mod(bytes,2)==0
%      result=dec2bin(fix(a0));  
%      p=repmat('0', bytes, 1);
%      for i=1:bytes+1
%         if length(result) <i
%            p(bytes-i+2:bytes)=result(1:i-1);
%            break;
%         end
%      end
%   end
%  end
%  function y = kahanSum(t, k,u,m,buf)
% 
%         a = (k * u) + m;
%     % y = 4096^(t-k)/a ; 
%     p=t-k;
%     result = 1;
%     base = buf(2);
% 
%     % 将 p 转换为二进制并遍历每一位
%     while p > 0
%         % 如果 p 的当前位为 1，则将 result 乘以 base 并取模 a*4096
%         b=a * buf(2);
%         if mod(p, 2) == 1
%             result = mod(result * base, b);
%         end
%         base = mod(base * base, b);
%         % 右移 p 的二进制表示
%         p = bitshift(p,-1);
%     end
%     y=result/a;
%     % end
%  end   