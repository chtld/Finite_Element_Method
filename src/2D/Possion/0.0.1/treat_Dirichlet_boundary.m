function [A,b] = treat_Dirichlet_boundary(A, b, boundarynodes,Dirichlet_fun,Pb)
nbn=size(boundarynodes,2);
for k=1:nbn
   if boundarynodes(1,k) ==-1
      i = boundarynodes(2,k);
      A(i,:) = 0;
      A(i,i) = 1;
      b(i) = feval(Dirichlet_fun,Pb(1,i),Pb(2,i));
   end
end
end