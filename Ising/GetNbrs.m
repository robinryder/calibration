function nbrs=GetNbrs(M,N,bc)
%nbrs=GetNbrs(M,N,bc) - build a (M*N) by 1 cell array giving nbrs on lattice
%INPUT
%   M,N - side lengths
%   bc - 'free', anything else is periodic BC's (refered to here and there
%         as cylindrical, but not in fact cylindrical) 
%OUTPUT
%   nbrs - cell array giving nbrs
if nargin==1, N=M; end

nbrs=cell(1,M*N);

FREE=1;
if (nargin==3 && ~strcmp(bc,'free')), FREE=0; end

if FREE
    %free
    for k=1:M*N
        [i,j]=ind2sub([M,N],k);
        subnbr=repmat([i;j],1,4)+[0 0 -1 1;1 -1 0 0];
        II=find(subnbr(1,:)>M | subnbr(1,:)<1);
        JJ=find(subnbr(2,:)>N | subnbr(2,:)<1);
        subnbr(:,union(II,JJ))=[];
        nbrs{k}=sub2ind([M,N],subnbr(1,:),subnbr(2,:));
    end
else
    %cylindrical A
    for k=1:M*N
        [i,j]=ind2sub([M,N],k);
        NS=[1*(i==M)+(i+1)*(i<M),M*(i==1)+(i-1)*(i>1),i,i];
        EW=[j,j,1*(j==N)+(j+1)*(j<N),N*(j==1)+(j-1)*(j>1)];
        nbrs{k}=sub2ind([M,N],NS,EW);
    end
    
%    %cylindrical B
%     for k=1:M*N
%         nbrs{k}=[M*N*(k==1)+(k>1)*(k-1),1*(k==M*N)+(k+1)*(k<M*N),...
%                  (k-M)*(k>M)+(M*N-(M-k))*(k<=M),(M-(M*N-k))*(k>M*(N-1))+(k+M)*(k<=M*(N-1))];
%     end 
end