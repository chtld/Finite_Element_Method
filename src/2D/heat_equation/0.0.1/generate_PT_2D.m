function [P,T] = generate_PT_2D(interval_x, interval_y, h ,basis_type)
Nx=(interval_x(2)-interval_x(1))/h(1);
Ny=(interval_y(2)-interval_y(1))/h(2);
if basis_type==201
    P = zeros(2,(Nx+1)*(Ny+1));
    T = zeros(3,2*Nx*Ny);
    for i=1:1:Nx+1
        x(i) = interval_x(1)+(i-1)*h(1);
        for j=1:1:Ny+1
            y(j) = interval_y(1)+(j-1)*h(2);
            P(:,(i-1)*(Ny+1)+j) = [x(i) y(j)]';
        end
    end
    ele=1;
    for i=1:1:Nx
        for j=2:1:Ny+1
            k=(i-1)*(Ny+1)+j;
            T(:,ele)=[k-1 k+Ny k]';
            T(:,ele+1)=[k k+Ny k+Ny+1]';
            ele=ele+2;
        end
    end

    
elseif basis_type==202
    P = zeros(2,(Nx+1)*(Ny+1));
    T = zeros(3,2*Nx*Ny);
    for i=1:1:Nx+1
        x(i) = interval_x(1)+(i-1)*h(1);
        for j=1:1:Ny+1
            y(j) = interval_y(1)+(j-1)*h(2);
            P(:,(i-1)*(Ny+1)+j) = [x(i) y(j)]';
        end
    end
    ele=1;
    for i=1:1:Nx
        for j=2:1:Ny+1
            k=(i-1)*(Ny+1)+j;
            T(:,ele)=[k-1 k+Ny k]';
            T(:,ele+1)=[k k+Ny k+Ny+1]';
            ele=ele+2;
        end
    end
else
    error('we have no this type basis');
end
end