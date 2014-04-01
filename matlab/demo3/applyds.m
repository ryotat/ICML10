function out=applyds(Xte, C)

out = pack_matrices(C.W)'*pack_matrices(Xte)+C.bias;


% $$$ function [Xv,sz]=pack_matrices(X)
% $$$ 
% $$$ for ii=1:length(X)
% $$$   [R,C,m]=size(X{ii});
% $$$   mm(ii)=m;
% $$$   nn(ii)=R*C;
% $$$   sz{ii}=[R,C];
% $$$ end
% $$$ 
% $$$ if length(unique(mm))>1
% $$$   error('Sample size mismatch');
% $$$ end
% $$$ 
% $$$ mm=mm(1);
% $$$ 
% $$$ Xv=zeros(sum(nn),mm);
% $$$ ix0=0;
% $$$ for ii=1:length(X)
% $$$   [R,C,m]=size(X{ii});
% $$$   I=ix0+(1:R*C);
% $$$   ix0=I(end);
% $$$   Xv(I,:) = reshape(X{ii},[R*C,m]);
% $$$ end
% $$$ 
