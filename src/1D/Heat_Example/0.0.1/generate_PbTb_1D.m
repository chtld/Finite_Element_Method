function [Pb,Tb] = generate_PbTb_1D(left, right, h,basis_type)
%GENERATE_PBTB 此处显示有关此函数的摘要
%   此处显示详细说明
if basis_type==101             %线性元
    N = (right-left)/h;
    Pb = zeros(1,N+1);
    Tb = zeros(2,N);
    for i=1:1:N
        x(i) = left+(i-1)*h;
        Pb(1,i) = x(i);
        Tb(:,i) =[i, i+1]' ;
    end
    x(N+1) = left+N*h;
    Pb(1,N+1) = x(N+1);
 
elseif basis_type==102     %二次元
    N = (right-left)/h;
    Nb = 2*N+1;
    Pb = zeros(1,Nb);
    Tb = zeros(3,N);
    for i=1:1:Nb
        Pb(i)=left+(i-1)*h/2;
    end
    for i=1:1:N
        Tb(:,i)=[2*i-1 2*i+1 2*i]'; 
    end
end

end

