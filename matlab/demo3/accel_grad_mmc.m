% September 25, 2009
% written by Shuiwang Ji and Jieping Ye

% This function implements the accelerated gradient algorithm for matrix
% classification described in Ji and Ye (ICML 2009). The original
% formulation is described in Tomioka and Aihara (ICML 2007) in which
% each input data sample consists of a matrix and the outputs are binary
% classes. Note that the current version of code only works for binary-class
% problems.

% References:
%Ji, S. and Ye, J. 2009. An accelerated gradient method for trace norm minimization. 
%In Proceedings of the 26th Annual international Conference on Machine Learning 
%(Montreal, Quebec, Canada, June 14 - 18, 2009). ICML '09, vol. 382. ACM, New York, 
%NY, 457-464.

%Tomioka, R. and Aihara, K. 2007. Classifying matrices with a spectral regularization. 
%In Proceedings of the 24th international Conference on Machine Learning 
%(Corvalis, Oregon, June 20 - 24, 2007). Z. Ghahramani, Ed. ICML '07, vol. 
%227. ACM, New York, NY, 895-902.


%[Wp,b,fval_vec,itr_counter] =
%accel_grad_mc(Xtrain,Ytrain,lambda,opt)

% required inputs:
% Xtrain: C x C x N array where each data point is a CxC matrix and N is the sample size
% Ytrain: N x 1 vector of class labels
% lambda: regularization parameter

% optional inputs:
% opt.L0: Initial guess for the Lipschitz constant
% opt.gamma: the multiplicative factor for Lipschitz constant
% opt.W_init: initial weight matrix
% opt.b_init: initial value for bias
% opt.epsilon: precision for termination
% opt.max_itr: maximum number of iterations
% opt.loss_prim: value of the primal objective function for termination

% outputs:
% Wp: the computed weight matrix
% b: the computed bias
% fval_vec: a vector for the sequence of function values
% itr_counter: number of iterations executed

function [Wp,b,fval_vec,time_vec,itr_counter,res] = accel_grad_mmc(Xtrain,Ytrain,lambda,opt)

if nargin<4
    opt = [];
end

if isfield(opt, 'L0')
    L0 = opt.L0;
else
    L0 = 100;
end

if isfield(opt, 'gamma')
    gamma = opt.gamma;
else
    gamma = 1.1;
end

if isfield(opt, 'W_init')
    W_init = opt.W_init;
else
  if isfield(opt,'blks')
    if ~isequal(size(Xtrain,2)*size(Xtrain,3),sum(prod(opt.blks, ...
                                                       2)))
      error('size of Xtrain does not match opt.blks');
    end
    W_init = cell(1,size(opt.blks,1));
    for jj=1:length(W_init)
      W_init{jj}=zeros(opt.blks(jj,:));
    end
  else
    W_init = {zeros(size(Xtrain,2),size(Xtrain,3))};
  end
end

if isfield(opt, 'b_init')
    b_init = opt.b_init;
else
    b_init = 0;
end

if isfield(opt, 'epsilon')
    epsilon = opt.epsilon;
else
    epsilon = 10^-5;
end

if isfield(opt, 'max_itr')
    max_itr = opt.max_itr;
else
    max_itr = 100;
end

if isfield(opt, 'loss_prim')
    loss_prim = opt.loss_prim;
else
    loss_prim = -1;
end

if ndims(Xtrain)>2
  Xtrain=Xtrain(:,:);
end

fval_vec = [];
time_vec = zeros(1,max_itr);
time0 = cputime;

W_old = W_init;
L = L0;
fval_old = rand(1,1);
c_init = b_init;
fval = loss_prim+100;
itr_counter = 0;
Z_old = W_old;
alpha = 1;
loss_prim = loss_prim+10^-2;

while (fval>loss_prim)&&(abs((fval_old-fval)/fval_old)>epsilon)&&(itr_counter<max_itr)
    itr_counter = itr_counter+1;
    fval_old = fval;
    [Wp,b1,P,sval,delta_W,delta_loss] = ComputeQP(Xtrain,Ytrain,Z_old,c_init,L,lambda);
    f = ComputeFun(Xtrain,Ytrain,Wp,b1);
    fval = f+lambda*sval;
    Q = P+lambda*sval;
    
    while fval>Q
        fprintf('Searching step size (fval = %f, Q = %f)...\n',fval,Q);
 %       fprintf('norm W = %f\n',norm(Wp,'fro'));
        L = L*gamma;
        [Wp,b1,P,sval,delta_W,delta_loss] = ComputeQP(Xtrain,Ytrain,Z_old,c_init,L,lambda);
        f = ComputeFun(Xtrain,Ytrain,Wp,b1);
        fval = f+lambda*sval;
        Q = P+lambda*sval;
    end

    fval_vec = [fval_vec,fval];
    time_vec(itr_counter) = cputime-time0;

    res(itr_counter) = 1 - ComputeDual(delta_loss, Xtrain, Ytrain, lambda, opt.blks)/fval;
    
    alpha_old = alpha;
    alpha = (1+sqrt(1+4*alpha_old^2))/2;
    for jj=1:length(Wp)
      Z_old{jj} = Wp{jj}+((alpha_old-1)/alpha)*(Wp{jj}-W_old{jj});
    end
    c_init = b1+((alpha_old-1)/alpha)*(b1-b_init);
    b_init = b1;
    W_old = Wp;
    
    if mod(itr_counter,50)==0
        fprintf('Iteration = %8d,  objective = %f gap = %g\n',itr_counter, ...
                fval, res(itr_counter));
    end
    
    if res(itr_counter)<1e-3
      break;
    end
end
b = b1;
return;

function [Wp,b1,P,sval,delta_W,delta_loss] = ComputeQP(X,Y,W,b,L,lambda)

[W1,b1,delta_W,delta_b,f,delta_loss] = ComputeGradStep(X,Y,W,L,b);

sval=0;
Wp=cell(size(W));
for jj=1:length(W)
  [U,D,V] = svd(W1{jj},0);

  %[U,D] = eig(W1);

  D = diag(D);
  D = D-(lambda/L);
  idx = find(D>0);
  sval =sval+sum(D(idx));

  Wp{jj} = U(:,idx)*diag(D(idx))*V(:,idx)';
end

%Wp = U(:,idx)*diag(D(idx))*U(:,idx)';

%P = f+trace(delta_W'*(Wp-W))+delta_b*(b1-b)+0.5*L*(norm(Wp-W,'fro')^2+(b1-b)^2);

%disp('Compute P');

P = f+delta_b*(b1-b)+(b1-b)^2;

for jj=1:length(W)
  %  ComputeProdTrace(delta_W,(Wp-W))+0.5*L*(norm(Wp-W,'fro')^2
  Wpdiff = Wp{jj}-W{jj};
  P = P + delta_W{jj}(:)'*Wpdiff(:)+0.5*L*sum(Wpdiff(:).^2);
end


return;

function [W1,b1,delta_W,delta_b,f,delta_loss] = ComputeGradStep(X,Y,W,L,b)

[delta_W,delta_b,f,delta_loss] = ComputeDerivative(X,Y,W,b);
W1=cell(size(W));
for jj=1:length(W)
  W1{jj} = W{jj}-(1/L)*delta_W{jj};
end
b1 = b-(1/L)*delta_b;
return;

function Wvec = Vectorize(W,n)

Wvec=zeros(n,1);
ix0=0;
for jj=1:length(W)
  I=ix0+(1:prod(size(W{jj})));
  ix0=I(end);
  Wvec(I)=W{jj}(:);
end


function [dev,delta_b,f,delta_loss] = ComputeDerivative(X,Y,W,b)

Wvec=Vectorize(W,size(X,2));
tmp = exp(-Y.*(X*Wvec+b)');

delta_loss = -Y.*tmp./(1+tmp);

dev = cell(size(W));

f      = sum(log(1+tmp));
delta_b= sum(delta_loss);

ix0=0;
for jj=1:length(dev)
  I=ix0+(1:prod(size(W{jj})));
  ix0=I(end);
  dev{jj} = zeros(size(W{jj}));
  dev{jj}(:) = delta_loss*X(:,I);
end


return;

function f = ComputeFun(X,Y,W,b)

Wvec=Vectorize(W,size(X,2));
f=sum(log(1+exp(-Y'.*(X*Wvec+b))));

return;

function dval = ComputeDual(gg, X, Y, lambda, blks)
aa=-gg';
aa=aa-mean(aa);

dnm = ds_dnorm(X'*aa,blks);

aa  = min(1, lambda/dnm)*aa;

dval = -loss_lrd(aa, Y');


