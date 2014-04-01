function C=trainds(Xtr, Ytr, lambda, varargin)

opt=propertylist2struct(varargin{:});
opt=set_defaults(opt, 'W0', [],...
                      'b0', 0,...
                      'eta', 1,...
                      'display',1,...
                      'tol', 1e-3,...
                      'solver', 'dal',...
                      'vv', []);

nn=0;
for ii=1:length(Xtr)
  nn=nn+size(Xtr{ii},1)*size(Xtr{ii},2);
end


if isstruct(opt.W0)
  opt.b0=opt.W0.bias;
  opt.W0=opt.W0.W;
end

if iscell(opt.W0)
  fprintf('trainds: lambda=%g Initial norm=%g\n',lambda, normc(opt.W0,'ds'));
  W0 = pack_matrices(opt.W0);
else  
  fprintf('trainds: lambda=%g\n',lambda);
  W0 = zeros(nn,1);
end


switch(opt.solver)
 case 'dal'
  [X,blks]=pack_matrices(Xtr);
  X=X';
  [ww,bias,status]=dallrds(W0,opt.b0,X,Ytr,lambda,'blks',blks,'eta',opt.eta,'tol',opt.tol,'display',opt.display); 

  status.dval=status.fval.*(1-status.res);

  C=struct('W', {unpack_matrices(ww,blks)},...
           'bias', bias,...
           'status', status);
 case 'ag'
  [X,blks]=pack_matrices(Xtr);
  X=X';
  if iscell(opt.W0) 
    opt_init={'W_init', {opt.W0}, 'b_init', opt.b0};
  else
   opt_init={'b_init',opt.b0};
  end
  
  [Wp,b,fval_vec,time_vec,itr_counter]=accel_grad_mmc(X,Ytr',lambda,struct('blks',blks,'epsilon',0,'max_itr',10000,opt_init{:}));
  status=struct('time',time_vec(1:itr_counter),...
                'fval',fval_vec,...
                'dval',evaldualobj(X,Ytr,pack_matrices(Wp),b,lambda,blks),...
                'niter',itr_counter);

  C=struct('W', {Wp},...
           'bias',b,...
           'status',status);
  case 'projgrad'
  if length(lambda)==1
    lambda = lambda*ones(1,length(Xtr));
  end
  
  problem = packvars(Xtr, Ytr', lambda, 'constraint', 'lr', 'ds');
  nblks=length(problem.ns);
  blks = [problem.ns; problem.nc]';
  
  xx0 = [W0; opt.b0];
 
  [xx, status] = feval(problem.solver, xx0, problem, 'epsf', opt.tol, ...
                       'display', 0,'maxiter',100000);

  bias = xx(end);
  xx = xx(1:end-1);
  
  C=struct('W', {unpack_matrices(xx,blks)},...
           'bias', bias,...
           'status', status);
end


% cls.v = vv;

fprintf('Elapsed time=%g\n', status.time(end));






