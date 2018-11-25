function [Gauss_coefficient_local_triangle,Gauss_point_local_triangle]=generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,vertices_triangle)
%Xiaoming He, 04/14/2009.
%Generate the Gauss coefficients and Gauss points on the local triangle by using affine tranformation
%Some details are in my 2 dimensional bilinear IFE notes(part 8)
%Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle: the Gauss coefficients and Gauss points on the reference triangle
%vertices_triangle: the coordinates of the vertices of the local triangle 
%Gauss_coefficient_local_triangle,Gauss_point_local_triangle:the Gauss coefficients and Gauss points on the local triangle.

x1=vertices_triangle(1,1);
y1=vertices_triangle(2,1);
x2=vertices_triangle(1,2);
y2=vertices_triangle(2,2);
x3=vertices_triangle(1,3);
y3=vertices_triangle(2,3);
Jacobi=abs((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1));
Gauss_coefficient_local_triangle=Gauss_coefficient_reference_triangle*Jacobi;
Gauss_point_local_triangle(:,1)=x1+(x2-x1)*Gauss_point_reference_triangle(:,1)+(x3-x1)*Gauss_point_reference_triangle(:,2);
Gauss_point_local_triangle(:,2)=y1+(y2-y1)*Gauss_point_reference_triangle(:,1)+(y3-y1)*Gauss_point_reference_triangle(:,2);