function Q = my_laplacian(n)
    %construct the nxn discrete laplacian matrix
    Q = circshift(eye(n),-1,2) + -2*eye(n) + circshift(eye(n),1,2);
    Q(1,end) = Q(1,end)-1; %delete unwanted 1 in top right corner
    Q(end,1) = Q(end,1)-1; %delete unwanted 1 in bottom right corner
end