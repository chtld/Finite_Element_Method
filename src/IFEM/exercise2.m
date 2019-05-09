node = [1,0; 1,1; 0,1; -1,1; -1,0; -1,-1; 0,-1; 0,0];  % coordinates
elem = [1,2,8; 3,8,2; 8,3,5; 4,5,3; 7,8,6; 5,6,8];     % connectivity
for i = 1:3
    [node,elem] = uniformrefine(node,elem);
end
showmesh(node,elem);