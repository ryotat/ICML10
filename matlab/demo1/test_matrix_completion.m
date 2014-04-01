function memo= ...
    test_matrix_completion(n,m,k,lambda,eta,nrep)



R=n;
C=n;

for ii=1:nrep
  W0=struct('U',randn(R,k),...
            'ss',(k:-1:1)',...
            'V',randn(C,k),...
            'D',[]);

  ind=randperm(R*C)';
  ind=ind(1:m);

  [I,J]=ind2sub([R,C],ind);

  fA=@(X)fobsX(X,I,J);
  fAT=@(aa)fobsXadj(aa,I,J,R,C);

  yy=fA(W0);

  [U,S,V]=lansvd(@(x)multlr(x,W0),@(x)multlrt(x,W0),R,C,k);

  ww=zeros(R,C);

  for jj=1:length(lambda)
    lmd=lambda(jj);
    profile on;
    [ww,stat]=dalsqds(ww,{fA,fAT,n,R*C},yy,lmd,'solver','qn','display',3,'eta',eta);

    pf=profile('info');
    names=getfieldarray(pf.FunctionTable, 'FunctionName');
    
    res=stat.res(end);
    time=stat.time(end);
    
    ix=find(strcmp(names,'lbfgs'));
    niter_out=pf.FunctionTable(ix).NumCalls;

    ix=find(strcmp(names,'lbfgs>linesearch_backtracking'));
    niter_in=pf.FunctionTable(ix).NumCalls;
    
    ix=find(strcmp(names,'lansvd'));
    nsvd=pf.FunctionTable(ix).NumCalls;
    
    rank=size(ww.U,2);
    dU=abs(U'*ww.U)-eye(k,rank);
    dV=abs(V'*ww.V)-eye(k,rank);
    error = sqrt(mean(dU(:).^2)+mean(dV(:).^2));  

    memo(ii,jj)=archive('lmd','ww','res','time','niter_out','niter_in','nsvd','rank','error');
  end
end
