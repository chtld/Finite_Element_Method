function [Pb,Tb] = generate_PbTb_2D(P, T, basis_type)
%GENERATE_PBTB 此处显示有关此函数的摘要
%   此处显示详细说明
if basis_type==201
    Pb=P;
    Tb=T;
elseif basis_type==102
N = (right-left)/h;
Nb = 2*N+1;
for i=1:1:Nb
   Pb(i)=left+(i-1)*h/2;

end
for i=1:1:N
      Tb(:,i)=[2*i-1 2*i+1 2*i]'; 
end
end

end

