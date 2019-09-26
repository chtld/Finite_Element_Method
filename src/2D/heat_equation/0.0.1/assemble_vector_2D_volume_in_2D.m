function b = assemble_vector_2D_volume_in_2D(vector_size, pde, mesh, FE, der_test_x,der_test_y, time)
%ASSEMBLE_VECTOR_1D_VOLUME_IN_1D 右端向量组装器
%   此处显示详细说明
    b = zeros(vector_size,1);
    P = mesh.P;
    T = mesh.T;
    Tb = FE.Tb;
    number_of_elements=size(T,2);%得到单元个数
    number_of_local_basis=get_number_of_local_basis(FE.basis_type_test);%得到测试函数的局部基函数个数
    [Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle]=generate_Gauss_reference_triangle(FE.Gauss_type);%得到参考单元的高斯权系数和高斯点
    for n=1: number_of_elements%循环每个单元
        vertices=P(:,T(:,n));%得到每个单元上点的坐标
        [Gauss_weight,Gauss_point]=generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,vertices);%得到当前单元的高斯权系数和高斯点
        for beta = 1:number_of_local_basis%循环测试函数的局部基函数
            int_value = Gauss_int_vector_test(pde.rhs, Gauss_weight, Gauss_point, vertices, FE.basis_type_test, beta, der_test_x, der_test_y, time);%计算高斯积分
             b(Tb(beta,n))=b(Tb(beta,n))+int_value;%加到全局的右端向量上
        end
    end


end

