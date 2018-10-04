function  FE  = generate_FE_space( FE, left, right, h )
%GENERATE_FE_SPACE 此处显示有关此函数的摘要
%   此处显示详细说明
N = (right - left)/h;                                                      %区间段数
basis_type = FE.basis_type_trial;                                          %基函数类型
[FE.P, FE.T] = generate_PT_1D(left, right, h,basis_type);                  %P矩阵存储节点编号及坐标，T矩阵存储每个单元的节点号，按一定的规则排序
[FE.Pb, FE.Tb] = generate_PbTb_1D(left, right, h,basis_type);              %Pb矩阵存储有限元节点及坐标，Tb矩阵存储每个单元的有限元节点号，按一定的规则排序
FE.boundarynodes = generate_boundarynodes(N, basis_type);                  %生成边界点，存储对应的边界点类型以及边界点的有限元编号

end

