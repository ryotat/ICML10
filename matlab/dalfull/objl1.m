function varargout=objl1(loss, xx, A, yy, lambda)

fnc = str2func(['loss_' loss 'p']);

[fval,gg]=fnc(A*xx, yy);

switch(loss)
case 'sq'
 hh=ones(size(yy));
case 'lr'
 pp=-gg.*yy;
 hh=pp.*(1-pp);
end

fval=fval+lambda*sum(abs(xx));
gg  =A'*gg+lambda*sign(xx);

I1=find(xx==0 & gg>0);
I2=find(xx==0 & gg<0);

gg(I1)=gg(I1)+max(-lambda, -gg(I1));
gg(I2)=gg(I2)+min(lambda, -gg(I2));


varargout{1}=fval;
if nargout>=2
  varargout{2}=gg;
  if nargout>=3
    H=A'*diag(hh)*A;
    varargout{3}=H;
  end
end
