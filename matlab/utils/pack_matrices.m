function [Xv,blks,mm]=pack_matrices(X)

blks=zeros(length(X),2);

for ii=1:length(X)
  [R,C,m]=size(X{ii});
  mm(ii)=m;
  nn(ii)=R*C;
  blks(ii,:)=[R,C];
end

if length(unique(mm))>1
  error('Sample size mismatch');
end

mm=mm(1);

Xv=zeros(sum(nn),mm);
ix0=0;
for ii=1:length(X)
  [R,C,m]=size(X{ii});
  I=ix0+(1:R*C);
  ix0=I(end);
  Xv(I,:) = reshape(X{ii},[R*C,m]);
end

