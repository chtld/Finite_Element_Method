clc
clear
close all
% 区域描述
left = 0;                                %区间左端点
right = 1;                               %区间右端点
% 计算时间描述
initial = 0;                             %初始时刻，终止时刻，时间步长
final = 0.2;
dt = 0.1;                                  
theta = 1;                               %在时间方向使用theta格式
% 方程信息
pde.exact_solution = @(x, t) 1 + x^2 + t;
pde.gradient = @(x, t) 2*x;
pde.initial_fun = @(x, t) 1 + x^2;
pde.right_hand_side = @(x, t) -1;
pde.coefficient1 = @(x, t) 1;
pde.coefficient2 = @(x, t) 1;
% 有限元空间信息
FE.Gauss_type = 4;                       %高斯积分的类型，选择几点的高斯积分
FE.basis_type_trial = 101;               %试探函数基函数类型
FE.basis_type_test = 101;                %测试函数基函数类型

fprintf('infinitenorm   L2norm   H1seminorm\n');
for i=2:1:7
    h=1/2^i;
    FE  = generate_FE_space( FE, left, right, h );
    [pde, FE] = Heat_Solver_1D(pde, FE, initial, final, dt, theta);
    error = compute_error( pde, FE, final );
end
 