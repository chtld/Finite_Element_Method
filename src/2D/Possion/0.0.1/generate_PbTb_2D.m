function [Pb, Tb, Nx, Ny] = generate_PbTb_2D(P, T, interval_x, interval_y, h, basis_type)
%GENERATE_PBTB 此处显示有关此函数的摘要
%   此处显示详细说明
if basis_type==201
    Pb=P;
    Tb=T;
    Nx=(interval_x(2)-interval_x(1))/h(1);
    Ny=(interval_y(2)-interval_y(1))/h(2);
elseif basis_type==202
    h = h/2;
    Nx=(interval_x(2)-interval_x(1))/h(1);
    Ny=(interval_y(2)-interval_y(1))/h(2);
    Pb = zeros(2,(Nx+1)*(Ny+1));
    Tb = zeros(6,Nx*Ny/4);
    for i=1:1:Nx+1
        x(i) = interval_x(1)+(i-1)*h(1);
        for j=1:1:Ny+1
            y(j) = interval_y(1)+(j-1)*h(2);
            Pb(:,(i-1)*(Ny+1)+j) = [x(i) y(j)]';
        end
    end
    ele=1;
    for i=1:2:Nx
        for j=2:2:Ny+1
            k=(i-1)*(Ny+1)+ j + 1;
            Tb(:,ele)=[k-2 k+2*Ny-2+2 k k-2+Ny+1 k-1+Ny+1 k-1]';
            Tb(:,ele+1)=[k k+2*Ny-2+2 k+2*Ny+2 k+Ny-1+1 k+2*Ny-1+2 k+Ny+1]';
            ele=ele+2;
        end
    end
end
end

