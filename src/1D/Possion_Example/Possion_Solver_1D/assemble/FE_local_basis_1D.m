function result = FE_local_basis_1D(x,vertices,basis_type,basis_index,der)    %nodel
h = vertices(2)-vertices(1);
if basis_type == 101 %101:1D linear nodal basis线性元
    if basis_index == 1    %第1个基函数
        if der ==0         %0阶导
           result = (vertices(2)-x)/h; 
        elseif der==1      %一阶导
            result = -1/h;
        elseif der>=2      %二阶导
            result = 0;
        else
            error('derivative order is wrong');
        end
    elseif basis_index ==2 %第二个基函数
        if der ==0         %0阶导
           result = (x-vertices(1))/h; 
        elseif der==1      %一阶导
            result = 1/h;
        elseif der>=2      %二阶导
            result = 0;
        else
            error('derivative order is wrong');
        end
    else
        error('wrong input for basis_index');
    end
elseif basis_type == 102    %二次元
    if basis_index == 1
       if der==0
           result = 2*((x-vertices(1))/h)^2-3*(x-vertices(1))/h+1;
       elseif der == 1
           result = 4*(x-vertices(1))/h^2-3/h;
       elseif der == 2
           result = 4/h^2;
       elseif der >= 3
           result = 0;
       else
           error('basis_index wrong');
       end
    elseif basis_index == 2
        
        if der == 0
            result = 2*((x-vertices(1))/h)^2-(x-vertices(1))/h;
        elseif der == 1
            result = 4*(x-vertices(1))/h^2-1/h;
        elseif der == 2
            result = 4/h^2;
        elseif der >= 3
            result = 0;
        else
            error('basis_index wrong');            
        end
        
    elseif  basis_index == 3
        if der == 0
            result = -4*((x-vertices(1))/h)^2+4*(x-vertices(1))/h;
        elseif der == 1
            result = -8*(x-vertices(1))/h^2+4/h;
        elseif der == 2
            result = -8/h^2;
        elseif der >= 3
            result = 0;
        else
            error('basis_index wrong');            
        end
    else
        error('no such basis_index');
    end
else
    error('basis_type wrong');
    
end
end