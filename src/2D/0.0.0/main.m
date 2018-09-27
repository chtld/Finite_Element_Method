clc
clear
close all
format long
%%
%给定区域
left = -1;                  %左边界
right = 1;                  %右边界
bottom = -1;                %下边界
top = 1;                    %上边界
Gauss_type=4;               %Gauss积分的类型
basis_type_trial = 201;     %试探函数的基函数类型
basis_type_test = 201;      %测试函数的基函数类型
%%
fprintf('infinitenorm  L2norm   H1seminorm\n');
for i=2:1:6
    hx=1/2^i;
    hy=1/2^i;
    [solution,error] = Poisson_solver_2D(left,right,bottom,top,hx,hy,basis_type_trial,basis_type_test,Gauss_type);
end
 