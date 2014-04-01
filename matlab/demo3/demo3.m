addpath ../dalfull/
addpath ../PROPACK/
addpath ../opt08/
addpath ../utils/

load 'bci.mat'

lambda=exp(linspace(log(10),log(0.001),20));

res1=xval_matrix(Xtr, Ytr, Xte, Yte, lambda, 'dal');
res2=xval_matrix(Xtr, Ytr, Xte, Yte, lambda, 'ag');

lambda_pg = 1./cell2mat(getfieldarray(res1,'nm'));
lambda_pg(1)=[];

res3=xval_matrix(Xtr, Ytr, Xte, Yte, lambda_pg, 'projgrad');



ix=1+(1:4:19);

C0=struct('status', struct('time',0,'niter',0,'dval',0,'fval',0));
r0=struct('lmd',nan,'solver','projgrad','C',C0,'out',[],'acc',nan,'nm',[],'time',0,'niter',0);


memo = [res2; [r0 res3]; res1];


lambda=cell2mat(getfieldarray(memo,'lmd')); lambda=lambda(1,:);
time=cell2mat(getfieldarray(memo,'time'));
niter=cell2mat(getfieldarray(memo,'niter'));
acc=cell2mat(getfieldarray(memo,'acc'));

subplot(1,3,1);
loglog(lambda', cumsum(niter([3,2,1],:)'),'-o','linewidth',2)
set(gca,'xdir','reverse','fontsize',16)
set(gca,'xtick',[0.001 0.01 0.1 1 10])
axis tight;
xlabel('Regularization constant')        
ylabel('# iterations')
grid on;

subplot(1,3,2);
loglog(lambda', cumsum(time([3,2,1],:)'),'-o','linewidth',2)
set(gca,'xdir','reverse','fontsize',16)
set(gca,'xtick',[0.001 0.01 0.1 1 10])
axis tight;
xlabel('Regularization constant')        
ylabel('Time (s)')
grid on;

subplot(1,3,3);
semilogx(lambda', acc([3,2,1],:)','-o','linewidth',2)
set(gca,'xdir','reverse','fontsize',16)
set(gca,'xtick',[0.001 0.01 0.1 1 10])
axis tight;
xlabel('Regularization constant')        
ylabel('Accuracy')
grid on;

%%%%%%

time=zeros(1,3); niter=zeros(1,3);

for ii=1:size(memo,2)
  show = ~isempty(find(ix==ii));

  lmd   = memo(3,ii).lmd;
  
  for jj=1:size(memo,1)
    stat  = memo(jj,ii).C.status;
    time(jj)  =time(jj)+stat.time(end);
    niter(jj) =niter(jj)+stat.niter;
    
    res   = (1-stat.dval(end)/stat.fval(end));
    acc   =memo(jj,ii).acc;
    
    if show
      if jj==1
        fprintf('%.3g \t', lmd);
      end
      
      fprintf('%g \t %s \t %s \t %g', niter(jj),...
              format_num_fixed(time(jj),3), format_num_fixed(res,3), acc);
    
      if jj<size(memo,1)
        fprintf('\t');
      end
    end
  end
  fprintf('\n');
end

