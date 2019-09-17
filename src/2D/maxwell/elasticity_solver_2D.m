function [result1, result2, solution, exact] = elasticity_solver_2D(pde, mesh, FE)
% 二维possion方程求解器
interval_x=[mesh.left mesh.right];
interval_y=[mesh.bottom mesh.top];
h=[mesh.hx,mesh.hy];
basis_type = FE.basis_type_trial;
Nx=(interval_x(2)-interval_x(1))/h(1);
Ny=(interval_y(2)-interval_y(1))/h(2);
[mesh.P, mesh.T] = generate_PT_2D(interval_x, interval_y, h, basis_type); %生成网格节点及编号
[FE.Pb, FE.Tb] = generate_PbTb_2D(mesh.P, mesh.T, basis_type);                               %生成有限元节点及编号
Tb = FE.Tb;
Pb = FE.Pb;
boundaryedges = generate_boundaryedges_2D_triangle(Nx,Ny);                 %生成边界边
boundarynodes = generate_boundarynodes_2D_triangle(Nx,Ny);                 %生成边界点
%-----------------------组装刚度矩阵----------------------------------------%
%组装刚度矩阵
matrix_size = [size(FE.Pb,2), size(FE.Pb,2)];                                    %初始化刚度矩阵
A1 = assemble_matrix_volume_in_2D(matrix_size, pde.coef.lambda, mesh, FE, 1, 0, 1, 0);
A2 = assemble_matrix_volume_in_2D(matrix_size, pde.coef.mu, mesh, FE, 1, 0, 1, 0);
A3 = assemble_matrix_volume_in_2D(matrix_size, pde.coef.mu, mesh, FE, 0, 1, 0, 1);
A4 = assemble_matrix_volume_in_2D(matrix_size, pde.coef.lambda, mesh, FE, 0, 1, 1, 0);
A5 = assemble_matrix_volume_in_2D(matrix_size, pde.coef.mu, mesh, FE, 1, 0, 0, 1);
A6 = assemble_matrix_volume_in_2D(matrix_size, pde.coef.lambda, mesh, FE, 1, 0, 0, 1);
A7 = assemble_matrix_volume_in_2D(matrix_size, pde.coef.mu, mesh, FE, 0, 1, 1, 0);
A8 = assemble_matrix_volume_in_2D(matrix_size, pde.coef.lambda, mesh, FE, 0, 1, 0, 1);
A9 = assemble_matrix_volume_in_2D(matrix_size, pde.coef.mu, mesh, FE, 0, 1, 0, 1);
A10 = assemble_matrix_volume_in_2D(matrix_size, pde.coef.mu, mesh, FE, 1, 0, 1, 0);
A = [A1+2*A2+A3 A4+A5;
     A6+A7 A8+2*A9+A10];
%-----------------------组装右端项量----------------------------------------%
%组装右端向量
vector_size = size(FE.Pb,2);                                                  %初始化右端向量
b1 = assemble_vector_2D_volume_in_2D(vector_size, pde.rhs1, mesh, FE, 0,0);
b2 = assemble_vector_2D_volume_in_2D(vector_size, pde.rhs2, mesh, FE, 0,0);
b = [b1; b2];
%----------------------处理边界条件-----------------------------------------%
%处理边界条件
%[ A, b ] = treat_Neumann_boundary( A, b, boundarynodes,'function_c','exact1',Pb);    
[A,b] = treat_Dirichlet_boundary(A, b, boundarynodes,pde.bdry,FE.Pb);


%------------------------求解线性方程组-------------------------------------%
%求解线性方程组
solution = A\b;
exact1 = compute_exact_solution(pde.exact_sol.u1, Pb,Tb);
exact2 = compute_exact_solution(pde.exact_sol.u2, Pb,Tb);
exact = [exact1; exact2];
%------------------------误差计算-------------------------------------------%
%error = compute_error_on_nodes('exact',Pb,solution);
result1 = compute_error(pde.exact_sol.u1, pde.exact_sol.u1_x, pde.exact_sol.u1_y,mesh, FE, solution(1:vector_size,:));
result2 = compute_error(pde.exact_sol.u2, pde.exact_sol.u2_x, pde.exact_sol.u2_y, mesh, FE, solution(vector_size+1:end,:));
end