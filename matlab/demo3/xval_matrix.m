function memo = xval_matrix(Xtr, Ytr, Xte, Yte, lambda, solver)


if ~iscell(Xtr)
  Xtr={Xtr};
  Xte={Xte};
end
% $$$ 
% $$$ 
% $$$ 
% $$$ C=cell(size(Xtr)+1);
% $$$ for jj=1:length(Xtr)
% $$$   [R,C,n]=size(Xtr{jj});
% $$$   C{jj}=zeros(R,C);
% $$$ end
% $$$ C{end}=0;

W0=[];
b0=0;
for ii=1:length(lambda)
  lmd = lambda(ii);
  C   = trainds(Xtr, Ytr, lmd, 'solver', solver, 'W0', W0, 'b0',b0);
  W0=C.W;
  b0=C.bias;
  out = applyds(Xte, C)';
  acc =  100*(1-mean(loss_0_1(Yte, out)));
  
  nm = normc(C,'ds');

  time = C.status.time(end);
  niter = C.status.niter;
    
  fprintf('lambda=%g (nm=%g) fval=%g dval=%g acc=%g\n', lambda(ii), ...
          nm, C.status.fval(end), C.status.dval(end), acc);

  memo(ii)=archive('lmd','solver','C','out','acc','nm','time','niter');
end



function loss=loss_0_1(yy, out)

loss = double(yy.*out<0);
loss(yy.*out==0)=0.5;
