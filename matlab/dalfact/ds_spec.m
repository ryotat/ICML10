function nm=ds_spec(ww,blks)

nn=min(size(ww.U,1),size(ww.V,1));
nm=zeros(nn,1);

ix0=0;
for kk=1:size(blks,1)
  R=size(ww.U,1);
  C=size(ww.V,1);
  rr=size(ww.U,2);
  J=ix0+(1:min(R,C));
  ix0=J(end);
  nm(J(1:rr))=ww.ss;
end
