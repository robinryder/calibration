function [D,x,nbrs]=ising(theta,L,SS,M,nbrs,x)
%[D,x,nbrs]=ising(theta,L,SS,M,nbrs,x) - sample Ising model using simple MCMC 
%INPUT
%   theta - scalar smoothing parameter exp(-theta*#x) where #x is num
%     disagreeing nbrs in image x
%   L - run length in single pixel updates
%   SS - subsample interval
%   M - lattice side 
%   nbrs - cell array giving nbrhood structure (periodic if omitted)
%   x - optional starting config (random if omitted)
%OUTPUT
%   D - a vector of L/SS #x values sampled along the run
%   x - final state
%   nbrs - the nbrhood structure used - useful output if omitted in call

%a trivial edit would allow rectangular

N = M; 
if nargin==4, nbrs=GetNbrs(M,N,'cylindrical'); end

%initialise MCMC X_0=x to agree with g where data is present
if nargin<6, x = round(rand(M,N)); end
%C=zeros(M,N,floor(L/SS));
D=zeros(1,floor(L/SS));
%figure(3);hf=imagesc(x);colormap(gray);axis square;drawnow;
%td=sum(sum(abs(diff(x'))))+sum(sum(abs(diff(x))));
td=hashX(x,nbrs);
for s=1:L
   pixel = ceil(rand*M*N); 
   hashfn=(x(nbrs{pixel})~=x(pixel));
   disagree = sum(hashfn);
   agree = sum(~hashfn);
   if log(rand) < theta*(disagree-agree)
      x(pixel) = 1-x(pixel); 
      td=td+(agree-disagree);
   end
   if rem(s,SS)==0,
      %set(hf,'CData',255*x); drawnow;
      %C(:,:,floor(s/SS))=x;
      D(1,floor(s/SS))=td;
   end
end


