clc
clear
close all
format long
%%
%网格信息
mesh.type = 'triangular';
mesh.left = 0;
mesh.right = 1;
mesh.bottom = 0;
mesh.top = 1;
mesh.hx = 0.25;
mesh.hy = 0.25;

%%
%方程信息
pde.exact_sol.u1 = @(x, y) exp(x + y - 1);
pde.exact_sol.u2 = @(x, y) exp(x - y + 1);
pde.exact_sol.u1_x = @(x, y) exp(x + y - 1);
pde.exact_sol.u1_y = @(x, y) exp(x + y - 1);
pde.exact_sol.u2_x = @(x, y) exp(x - y + 1);
pde.exact_sol.u2_y = @(x, y) -exp(x - y + 1);
pde.rhs1 = @(x, y) (1.75 .* exp(1 + x - y) - 3.75 .* exp(-1 + x + y));
pde.rhs2 = @(x, y) (-3.75 .* exp(1 + x - y) - 1.75 .* exp(-1 + x +y));
pde.coef.lambda = @(x, y) 0.75;
pde.coef.mu = @(x, y) 1;
pde.bdry.num = 1;
s1_1 = @(x, y) exp(y - 1);
s2_1 = @(x, y) exp(y);
s3_1 = @(x, y) exp(x - 1);
s4_1 = @(x, y) exp(x);
b1_1 = @(x, y) ((x == mesh.left && y ~= mesh.bottom && y ~= mesh.top) .* s1_1(x, y));
b2_1 = @(x, y) ((x == mesh.right && y ~= mesh.bottom && y ~= mesh.top) .* s2_1(x,y));
b3_1 = @(x, y) ((y == mesh.bottom) .* s3_1(x, y));
b4_1 = @(x, y) ((y == mesh.top) .* s4_1(x, y));
pde.bdry.u1_dirichlet = @(x, y) (b1_1(x, y) + b2_1(x, y) + b3_1(x, y) + b4_1(x, y));
s1_2 = @(x, y) exp(1 - y);
s2_2 = @(x, y) exp(2 - y);
s3_2 = @(x, y) exp(x + 1);
s4_2 = @(x, y) exp(x);
b1_2 = @(x, y) ((x == mesh.left && y ~= mesh.bottom && y ~= mesh.top) .* s1_2(x, y));
b2_2 = @(x, y) ((x == mesh.right && y ~= mesh.bottom && y ~= mesh.top) .* s2_2(x, y));
b3_2 = @(x, y) ((y == mesh.bottom) .* s3_2(x, y));
b4_2 = @(x, y) ((y == mesh.top) .* s4_2(x, y));
pde.bdry.u2_dirichlet = @(x, y) (b1_2(x, y) + b2_2(x, y) + b3_2(x, y) + b4_2(x, y));
%%
FE.basis_type_test = 201;
FE.basis_type_trial = 201;
FE.Gauss_type = 4;

%%
ns = [4, 8, 16, 32, 64, 128];
fprintf('h        L_inf_err         L2_err        H1_err          ');
fprintf('h        L_inf_err         L2_err        H1_err\n');
for i = 1: length(ns)
    mesh.hx = 1.0 / ns(i);
    mesh.hy = 1.0 / ns(i);
    [result1, result2, solution, exact] = elasticity_solver_2D(pde, mesh, FE);
    fprintf('1/%d  \t%e   \t%e   \t%e', ns(i), result1.error.L_inf, result1.errorl.L2, result1.error.H1);
    fprintf('\t1/%d  \t%e   \t%e   \t%e\n', ns(i), result2.error.L_inf, result2.errorl.L2, result2.error.H1);
end
%%
num = size(solution, 1)/2;
x = 0: mesh.hx: 1;
y = 0: mesh.hy: 1;
[X, Y] = meshgrid(x, y);
z11 = reshape(exact(1:num,:),size(X));
z1 = reshape(solution(1:num,:),size(X));
figure(1);
subplot(121);
surf(x, y, z1);
subplot(122);
surf(x, y, z11);
figure(2);
z2 = reshape(solution(num+1: end,1),size(X));
z22 = reshape(exact(num+1: end,1),size(X));
subplot(121);
surf(x, y, z2);
subplot(122);
surf(x, y, z22);
