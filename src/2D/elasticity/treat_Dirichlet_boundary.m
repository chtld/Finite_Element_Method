function [A,b] = treat_Dirichlet_boundary(A, b, boundarynodes, bdry, Pb)
nbn=size(boundarynodes,2);
row = 0.5*size(A,1);
for k=1:nbn
   if boundarynodes(1,k) ==-1
      i = boundarynodes(2,k);
      A(i,:) = 0; A(i + row,:) = 0;
      A(i,i) = 1; A(i + row,i + row) = 1;
      b(i) = feval(bdry.u1_dirichlet, Pb(1,i), Pb(2,i));
      b(i + row) = feval(bdry.u2_dirichlet, Pb(1,i), Pb(2,i));
   end
end
end