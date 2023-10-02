function [c3] = c3cum(K,L3,x)
%L3 =4;
%x= round(10* randn(1,100)) ; 
%K = 10;
N = length(x);

M = round(N/K) ;
A1 = [] ;
z = reshape(x,M,K)';
u =zeros(100,M);
r = zeros(2*L3 +1 , 2*L3 +1 , K);
c3 = zeros(2*L3 +1 , 2*L3 +1);

for i=1:K
    %u(i,:)=x((i-1)* M + 1 : i*M ) ;
    if (((i*M) - ((i-1) * M +1)) > 1 )
        
    m = mean(x(((i-1)*M + 1):(i*M) ));
    
    A1 = [A1 ; x(((i-1)*M + 1 ):(i*M)) - m];
    
    else
    
    A1 = [A1 ; x(((i-1)*M + 1 ):(i*M)) ];
    
    end
    
end


for i=1:K
  for t1 = -L3:1:L3
    for t2 = -L3:1:L3
        s1 = max(1,max(1-t1,1-t2));
        s2 = min(M,min(M-t1 , M-t2));
        %syms k ;
        %f = A(i,k)*A(i,k+t1)*A(i,k+t2) ;
        %r(t1 + L3 + 1 , t2 + L3 + 1 , i) = (1/M) * symsum(A(k,1),k,s1,s2);
        for l=s1:1:s2
            r(t1 + L3 + 1 , t2 + L3 + 1 , i)= r(t1 + L3 + 1 , t2 + L3 + 1 , i) + A1(i,l)*A1(i,l+t1)*A1(i,l+t2);
        end
        r(t1 + L3 + 1 , t2 + L3 + 1 , i) = (1/M) * r(t1 + L3 + 1 , t2 + L3 + 1 , i) ;
    end
  end
end

 c3 = (1/K)*sum(r,3);
 
end