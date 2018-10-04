function A = assemble_matrix_volume_1D(coef_fun, FE, matrix_size, der_trial, der_test, t)
    A = zeros(matrix_size(1),matrix_size(2));%初始化刚度矩阵
    T = FE.T;
    P = FE.P;
    Tb = FE.Tb;
    Gauss_type = FE.Gauss_type;
    basis_type_test = FE.basis_type_test;
    basis_type_trial = FE.basis_type_trial;
    number_of_elements=size(T,2);            %T矩阵的列数即单元个数
    number_of_local_basis=get_number_of_local_basis(basis_type_test);%局部基函数的个数
    for n=1: number_of_elements              %循环每个单元
        vertices=P(:,T(:,n));                %得到每个单元对应的节点坐标，即区间的左右端点坐标
        [Gauss_weight,Gauss_point]=generate_Gauss_formula(vertices, Gauss_type);%得到该区间的高斯点和高斯权系数
        for alpha = 1:number_of_local_basis                           %循环测试函数（test function）
            for beta = 1:number_of_local_basis                        %循环试探函数(trial function)
                %计算高斯积分，需要输入高斯点，高斯权，试探函数类型及导数阶，测试函数类型及导数阶以及它们是第几个局部基函数
                int_value = Gauss_vol_int_trial_test(coef_fun, Gauss_weight,Gauss_point,vertices,basis_type_trial,alpha,der_trial,...
                                                                                               basis_type_test,beta,der_test, t);
                %将每一次的结果加入总刚度矩阵，由Tb矩阵找到它们在总刚度矩阵的位置
                A(Tb(beta,n),Tb(alpha,n))=A(Tb(beta,n),Tb(alpha,n))+int_value;
            end
        end
    end
end