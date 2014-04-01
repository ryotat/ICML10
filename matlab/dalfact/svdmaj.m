function [U,S,V]=svdmaj(A, lambda, varargin)
opt=propertylist2struct(varargin{:});
opt=set_defaults(opt, 'kinit', 10, 'kstep', 2);

opt.kinit=max(opt.kinit,10);


if isnumeric(A)
  MM=min(size(A));
else
  MM=min(A{3:4});
end


if iscell(A) || max(size(A))>=100  && opt.kinit<MM
  mm = inf;

  fprintf('[svdmaj]\n');

  kk=opt.kinit/opt.kstep;
  while mm>lambda && kk<MM
    kk=min(kk*opt.kstep, MM);
    if iscell(A)
      [U,S,V]=lansvd(A{:},kk,'L',struct('s_target',lambda));
    else
      [U,S,V]=lansvd(A, kk,'L',struct('s_target',lambda));
    end
    mm=min(diag(S));
    fprintf('kk=%d lambda=%g mm=%g\n',kk,lambda,min(diag(S)));
  end
else
  [U,S,V]=svd(A);
end

