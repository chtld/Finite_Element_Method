function [pde, FE] = Heat_Solver_1D( pde, FE, initial, final, dt, theta)
% 一维热方程求解器
% 主要来求解如下形式的方程
% $u_t-\nabla \cdot c(x) \nabla u(x) = f(x),a\le x\le b$
%$u(a)=g_a$,
%$u'(b)+q_bu(b)=p_b$
Mt = (final - initial)/dt;
Pb = FE.Pb;
boundarynodes = FE.boundarynodes;                                      %基函数类型
%------------------------------组装器--------------------------------------%
%矩阵组装器
matrix_size = [size(Pb,2), size(Pb,2)];                                    %刚度矩阵的大小
x_now = generate_initial_vector(pde, Pb, initial);
vector_size = size(Pb, 2);                                                 %右端向量的大小
for i = 0: 1: Mt-1
    t = initial+(i+1)*dt;
    M = assemble_matrix_volume_1D(pde.coefficient1 , FE, matrix_size, 0, 0, t-dt);
    %-----------------------------矩阵组装器--------------------------------%
    A_now = assemble_matrix_volume_1D(pde.coefficient2, FE, matrix_size, 1, 1, t-dt);
    A_next = assemble_matrix_volume_1D(pde.coefficient2, FE, matrix_size, 1, 1, t);
    %向量组装器
    b_now = assemble_vector_1D(pde.right_hand_side, FE, vector_size, 0, t-dt);
    b_next = assemble_vector_1D(pde.right_hand_side, FE, vector_size, 0, t);
    %----------------------------解线性方程组-------------------------------%
    A_tilde = M/dt + theta*A_next;
    b_tilde = theta*b_next+(1-theta)*b_now+M/dt*x_now-(1-theta)*A_now*x_now;
    %----------------------------处理边界条件-------------------------------%
    [A_tilde, b_tilde] = treat_Dirichlet_boundary(A_tilde, b_tilde, boundarynodes, pde.exact_solution, Pb, t);           %处理Dirichlet边界
    x_next = A_tilde\b_tilde;
    x_now = x_next;
end
solution  = x_next;
pde.solution = solution;
end