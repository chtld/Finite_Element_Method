function [ result ] = FE_local_basis_2D( x,y,vertices,basis_type,basis_index,der_x,der_y )
%FE_LOCAL_BASIS_2D 此处显示有关此函数的摘要
%   此处显示详细说明
xn1=vertices(1,1);
xn2=vertices(1,2);
xn3=vertices(1,3);
yn1=vertices(2,1);
yn2=vertices(2,2);
yn3=vertices(2,3);
Jacobi=(xn2-xn1)*(yn3-yn1)-(xn3-xn1)*(yn2-yn1);
x_hat=((yn3-yn1)*(x-xn1)-(xn3-xn1)*(y-yn1))/Jacobi;
y_hat=(-(yn2-yn1)*(x-xn1)+(xn2-xn1)*(y-yn1))/Jacobi;
if der_x==0&&der_y==0
    result =  FE_reference_basis_2D( x_hat,y_hat,basis_type,basis_index,der_x,der_y );
elseif der_x==1&&der_y==0
    result = (yn3-yn1)/Jacobi*FE_reference_basis_2D( x_hat,y_hat,basis_type,basis_index,der_x,der_y )...
            +(yn1-yn2)/Jacobi*FE_reference_basis_2D(x_hat,y_hat,basis_type,basis_index,der_y,der_x);
elseif der_y==1&&der_x==0
    result = (xn1-xn3)/Jacobi*FE_reference_basis_2D(x_hat,y_hat,basis_type,basis_index,der_y,der_x)...
            +(xn2-xn1)/Jacobi*FE_reference_basis_2D(x_hat,y_hat,basis_type,basis_index,der_x,der_y);
elseif der_x+der_y>=2
    result = 0;
else
    error('you are wrong');
end
end

