addpath ../dalfull/
addpath ../PROPACK/
addpath ../opt08/
addpath ../utils/

m = 1000;      % Number of samples
n = 64;        % Dimension
k = 16;        % True rank
lambda = 800;  % Regularization constant
algos  = {'dal','lrds','ag','projgrad',};
               % Algorithms
nrep   = 10;   % Number of repetitions

memo = test_matrix_classification(m, n, k, lambda, algos, nrep);



memo=memo(:,[1,4,2,3]); 


timec=getfieldarray(memo,'time');
resc=getfieldarray(memo,'res');

ns=100;
for ii=1:size(memo,2)
  rr=mean(cell2mat(foreach(@rangeof,resc(:,ii))));
  res{ii}=exp(linspace(log(rr(2)),log(rr(1)),ns));

  time_tmp=nan*ones(size(memo,1),ns);
  for jj=1:ns
    ix=foreach(@(x)min(find(x<res{ii}(jj))),resc(:,ii));
    for kk=1:length(ix)
      if ~isempty(ix{kk})
        time_tmp(kk,jj)=timec{kk,ii}(ix{kk});
      end
    end
  end
  mtime{ii}=nanmean(time_tmp,1);
  stime{ii}=nanstd(time_tmp,[],1);
end


figure
col=get(gca,'colororder');
for ii=1:length(mtime)
  hold on;

  colsat=([1 1 1]+col(ii,:))/2;
  
  upper=mtime{ii}+stime{ii};
  lower=mtime{ii}-stime{ii}; 
  patch([lower, fliplr(upper)],[log(res{ii}), log(fliplr(res{ii}))],  colsat,'EdgeColor','none')

  plot(mtime{ii}, log(res{ii}), 'color', col(ii,:), 'linewidth',2);
end

set(gca,'fontsize',16);
grid on;
ytic=[1e-5, 1e-4, 0.001, 0.01, 0.1, 1];
set(gca,'ytick',log(ytic), 'yticklabel',num2cell(ytic))
h=get(gca,'children')
legend(h([1,3,5,7]),{'AG','IP','PG','DAL'});
xlabel('CPU time (s)');
ylabel('Relative duality gap');


