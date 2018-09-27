function boundaryedges = generate_boundaryedges_2D_triangle(Nx,Ny)
boundaryedges=zeros(4,2*Nx+2*Ny);
k=1;
for i=1:Nx
    boundaryedges(1,k)=-1;
    boundaryedges(2,k)=1+2*Ny*(i-1);
    boundaryedges(3,k)=1+(Ny+1)*(i-1);
    boundaryedges(4,k)=1+(Ny+1)*i;
        k=k+1;
end
for j=1:Ny
    boundaryedges(1,k)=-1;
    boundaryedges(2,k)=1+2*Ny*(Nx-1)+1+(j-1)*2;
    boundaryedges(3,k)=1+(Ny+1)*Nx+(j-1);
    boundaryedges(4,k)=1+(Ny+1)*Nx+j;
    k=k+1;
end
for i=1:Nx
    boundaryedges(1,k)=-1;
    boundaryedges(2,k)=1+2*Ny*(Nx-1)+1+(Ny-1)*2-(i-1)*2*Ny;
    boundaryedges(3,k)=1+(Ny+1)*Nx+Ny-(i-1)*(Ny+1);
    boundaryedges(4,k)=1+(Ny+1)*Nx+Ny-i*(Ny+1);
    k=k+1;
end
for j=1:1:Ny
    boundaryedges(1,k)=-1;
    boundaryedges(2,k)=1+2*Ny*(Nx-1)+1+(Ny-1)*2-(Nx-1)*2*Ny-1-(j-1)*2;
    boundaryedges(3,k)=1+(Ny+1)*Nx+Ny-Nx*(Ny+1)-(j-1);
    boundaryedges(4,k)=1+(Ny+1)*Nx+Ny-Nx*(Ny+1)-j;
    k=k+1;
end
end