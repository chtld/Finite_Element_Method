function [ A, b ] = treat_Neumann_boundary( A, b, boundarynodes,coef_fun,Neumann_fun,Pb)
nbn=size(boundarynodes,2);
for k=1:nbn
   if boundarynodes(1,k) ==-2
      i = boundarynodes(2,k);
      b(i) = b(i)+feval(Neumann_fun,Pb(i))*feval(coef_fun,Pb(i));
   end
end

end

