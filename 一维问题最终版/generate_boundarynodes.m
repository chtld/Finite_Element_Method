function boundarynodes = generate_boundarynodes(N,basis_type)
if basis_type == 101
    boundarynodes(1,1)=-1; 
    boundarynodes(2,1)=1; 
    boundarynodes(1,2)=-2; 
    boundarynodes(2,2)=N+1;
elseif basis_type == 102
    boundarynodes(1,1)=-1; 
    boundarynodes(2,1)=1; 
    boundarynodes(1,2)=-2; 
    boundarynodes(2,2)=2*N+1; 
else
    error('we have no this type basis!');
end
end