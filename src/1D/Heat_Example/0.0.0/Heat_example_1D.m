clc
clear
close all
%其实我想在外面定义边界条件，在外面定义区域以及PT矩阵等，因为如果用其他软件
%生成PT矩阵必须从这里就开始导入，边界条件也应该是在这里就定义好然后传给求解器
%这一切都属于用户的需求，都需要放在外面，但PbTb矩阵的形成却涉及到了有限元空间
%的概念，因此需要在内部完成
%初始化各项参数
left = 0;                          %区间左端点
right = 1;                         %区间右端点

initial = 0;                       %初始时刻，终止时刻，时间步长
final = 0.2;
dt = 0.1;
                                   %在时间方向使用theta格式
theta = 1.0;

Gauss_type = 4;                    %高斯积分的类型，选择几点的高斯积分
basis_type_trial = 101;            %试探函数基函数类型
basis_type_test = 101;             %测试函数基函数类型
der_trial = 1;                     %试探函数导数阶
der_test = 1;                      %测试函数的导数阶
fprintf('infinitenorm   L2norm   H1seminorm\n');
for i=2:1:7
    h=1/2^i;
    [solution,error] = Heat_Solver_1D(left, right, h, initial, final, dt, theta,basis_type_trial,basis_type_test,der_trial,der_test,Gauss_type);
end
 