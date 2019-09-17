function A = assemble_matrix_volume_in_2D(matrix_size, coef, mesh, FE, der_trial_x,der_trial_y,der_test_x,der_test_y)
%矩阵组装器,需要输入系数函数,矩阵大小,P,T,Tb矩阵,试探函数的基函数类型及对x,y的导数阶,测试函数的基函数类型及导数阶
    A = sparse(matrix_size(1),matrix_size(2));
    P = mesh.P;
    T = mesh.T;
    Tb = FE.Tb;
    number_of_elements=size(T,2);%得到单元个数
    number_of_local_basis=get_number_of_local_basis(FE.basis_type_test);%得到局部的基函数的个数
    [Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle]=generate_Gauss_reference_triangle(FE.Gauss_type);%得到参考单元的高斯权系数和高斯点
    for n=1: number_of_elements%循环每个单元
        vertices=P(:,T(:,n));  %得到每个单元上的点的坐标
        [Gauss_weight,Gauss_point]=generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,vertices);%得到当前三角形单元的高斯权系数和高斯点
        for alpha = 1:number_of_local_basis%循环试探函数
            for beta = 1:number_of_local_basis%循环测试函数
                int_value = Gauss_vol_int_trial_test_2D(coef,Gauss_weight,Gauss_point,vertices,FE.basis_type_trial,alpha,der_trial_x,der_trial_y,...
                                                                                               FE.basis_type_test,beta,der_test_x,der_test_y);%计算高斯积分
                A(Tb(beta,n),Tb(alpha,n))=A(Tb(beta,n),Tb(alpha,n))+int_value;
            end
        end
    end
end