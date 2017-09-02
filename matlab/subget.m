function [row,col] = subget(m,j,n,k, dim)
% figures out indexing given m,j,n,k where m,n dictate columns and j,k
% dictate columns


col = (m-1)*dim^2 + n;
row = (j-1)*dim^2 + k;
end