function [P,T] = generate_PT_1D(left, right, h,basis_type)
if basis_type==101              %一维线性元
    N = (right-left)/h;         %单元个数
    P = zeros(1,N+1);           %列数代表N+1个节点
    T = zeros(2,N);             %列数代表N个单元
    for i=1:1:N                 %循环每个单元
        x(i) = left+(i-1)*h;    %计算节点坐标
        P(1,i) = x(i);          %将节点坐标放入P矩阵的每一列
        T(:,i) =[i, i+1]' ;     %将每个单元的节点号放入T矩阵的每一列
    end
    x(N+1) = left+N*h;          %最后一个节点
    P(1,N+1) = x(N+1);

    
elseif basis_type==102          %一维二次元
    N = (right-left)/h;         %单元个数
    P = zeros(1,N+1);           %列数代表N+1个节点
    T = zeros(2,N);             %列数代表N个单元
    for i=1:1:N                 %循环每个单元
        x(i) = left+(i-1)*h;    %计算节点坐标
        P(1,i) = x(i);          %将节点坐标放入P矩阵的每一列
        T(:,i) =[i, i+1]' ;     %将每个单元的节点号放入T矩阵的每一列
    end
    x(N+1) = left+N*h;          %最后一个节点
    P(1,N+1) = x(N+1);
else
    error('we have no this type basis');
end



end