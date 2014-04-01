function memo=test_matrix_classification(m, n, k, lambda, algos, nrep)

ix1=find(strcmp(algos,'lrds'));
ix2=find(strcmp(algos,'ag'));
ix3=find(strcmp(algos,'projgrad'));

if ~isempty(ix2) && (isempty(ix1) || ix2<ix1)
  error('AG method must be called after lrds');
end

if ~isempty(ix3) && (isempty(ix1) || ix3<ix1)
  error('projgrad method must be called after lrds');
end


for ii=1:nrep
  A=zeros(m,n,n);

  for kk=1:m
    X=randn(n,n);
    A(kk,:,:)=X'*X;
  end

  W0=randn(n,n);
  W0=(W0+W0')/2;

  [U,D]=eig(W0); I=[1:k/2, n:-1:n-k/2+1];

  W0=U(:,I)*D(I,I)*U(:,I)';

  yy = sign(A(:,:)*W0(:));

  %% Apply algorithms
  for jj=1:length(algos)
    algo=algos{jj};
    
    fprintf('Running algorithm [%s]\n',algo);
    switch(algo)
     case 'dal'
      [ww,bb,stat]=dallrds(zeros(n,n),0,A(:,:),yy,lambda);
      spec =svd(reshape(ww,[n,n]));
      res  =stat.res;
      time =stat.time;
      fval =stat.fval;
     case 'lrds'
      [W,bias,z,stat]=lrds_dual(permute(A,[2,3,1]), yy', lambda,'relgap',1,'tol',1e-3);
      spec=svd(W);
      res =stat.gap;
      time=stat.time;
      fval=stat.fval;
     case 'ag'
      opt=struct('max_itr',1000);
      [Wp,b,fval_vec,time_vec,itr_counter]=accel_grad_mc(permute(A,[2,3, 1]),yy',lambda,opt);
      %[Wp,b,fval_vec,time_vec,itr_counter]=accel_grad_mc(A,yy',lambda,opt);
      spec=svd(Wp);
      res = 1+stat.obj(end)./fval_vec;
      if length(time_vec)>length(fval_vec)
        time_vec(length(fval_vec)+1:end)=[];
      end
      time=time_vec;
      fval=fval_vec;
     case 'projgrad'
      problem = packvars(permute(A,[2,3,1]), yy', 1/sum(svd(W)),...
                         'constraint', 'lr', 'ds');

      xx0=zeros(n*n+1,1);
      
      [xx, stat] = projgrad(xx0, problem, 'epsf',1e-3,'display',2);

      spec =svd(reshape(xx(1:end-1),[n,n]));
      res  =stat.res;
      time =stat.time;
      fval =stat.fval;
    end
    fprintf('Algorithm [%s] converged\n',algo);

    memo(ii,jj)=archive('m','n','k','lambda','algo','spec','time','res','fval');
  end
end
