function X=unpack_matrices(Xv, blks)

mm=size(Xv,2);

nblks=size(blks,1);

X=cell(1,nblks);

ix0=0;
for ii=1:nblks
  sz=blks(ii,:);
  I=ix0+(1:prod(sz));
  ix0=I(end);
  X{ii}=reshape(Xv(I,:),[sz,mm]);
end
