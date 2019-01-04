function [ td ] = hashX( x, nbrs )
%hashX() - calulate #x the number of disagreeing neighbors
%   INPUT
%   x - binary array
%   nbrs - cell array giving nbrs
%   OUTPUT
%   td - td is #x the number of disagreeing nbrs

M=size(x,1); N=size(x,2); 
td=0; for (k=1:M*N), td=td+sum(x(nbrs{k})~=x(k)); end; td=td/2;

end

