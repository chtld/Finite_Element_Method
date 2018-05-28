% clc
% clear
% close all
format long
%初始化各项参数
left = 0;                          %区间左端点
right = 1;                         %区间右端点
Gauss_type=4;                      %高斯积分的类型，选择几点的高斯积分
basis_type_trial = 101;            %试探函数基函数类型
basis_type_test = 101;             %测试函数基函数类型
der_trial = 1;                     %试探函数导数阶
der_test = 1;                      %测试函数的导数阶

fprintf('infinitenorm   L2norm   H1seminorm\n');
for i=2:1:7
    h=1/2^i;
    [solution,error] = Poisson_solver_1D(left,right,h,basis_type_trial,basis_type_test,der_trial,der_test,Gauss_type);
end
 