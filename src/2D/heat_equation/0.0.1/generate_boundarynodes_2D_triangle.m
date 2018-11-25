function boundarynodes = generate_boundarynodes_2D_triangle(Nx,Ny)
boundarynodes=zeros(2,2*Nx+2*Ny);
k=1;
for i=1:Nx
    boundarynodes(1,k)=-1;
    boundarynodes(2,k)=1+(Ny+1)*(i-1);
    k=k+1;
end
for j=1:Ny
    boundarynodes(1,k)=-1;
    boundarynodes(2,k)=1+(Ny+1)*Nx+(j-1);
    k=k+1;
end
for i=1:Nx
    boundarynodes(1,k)=-1;
    boundarynodes(2,k)=1+(Ny+1)*Nx+Ny-(i-1)*(Ny+1);
    k=k+1;
end
for j=1:1:Ny
    boundarynodes(1,k)=-1;
    boundarynodes(2,k)=1+(Ny+1)*Nx+Ny-Nx*(Ny+1)-(j-1);
    k=k+1;
end
end