function plotsvd(ww, blks)

mm=size(blks,1);

ix0=0;
for ii=1:mm
  I=ix0+(1:prod(blks(ii,:)));
  ix0=I(end);
  subplot(mm,1,ii);
  plot(svd(reshape(ww(I),blks(ii,:))),'-x');
end
