addpath ../dalfact/
addpath ../PROPACK/
addpath ../utils/


n=10000;    % size of the unknown matrix
m=1200000;  % number of observations
k=10;       % true rank of the unknown matrix
lambda = [1000 700 500 300 200 150 100];
            % regularization constant
eta=10;     % step-size parameter
nrep=10;    % number of repetition

memo=test_matrix_completion(n,m,k,lambda,eta,nrep);

% m=2400000;
% k=20;
% lambda=[2000 1500 1200 1000 700 500 400 300 200 150 120 100];
% 
% memo=test_matrix_completion(n,m,k,lambda,eta,nrep);



% Print out table

[nrep,nlmd]=size(memo);

flds = {'lmd','time','niter_out','niter_in','nsvd','rank','error'};

for jj=1:length(flds)
  assignin('caller',flds{jj},cell2mat(getfieldarray(memo,flds{jj})));
end

lambda = lmd(1,:);

time      = cumsum(time')';
niter_out = cumsum(niter_out')';
niter_in  = cumsum(niter_in')';
nsvd      = cumsum(nsvd')';


fprintf('lambda\t time(s)\t #iter_out\t #iter_in\t #svd\t\t rank\t\t error\n');
fprintf('----------------------------------------------------------------------------------------------------------------\n');
for ii=1:nlmd
  fprintf('%g\t %s\t %s\t %s\t %s\t %s\t %s\n',...
          lambda(ii),...
          format_num_std(time(:,ii)),...
          format_num_std(niter_out(:,ii)),...
          format_num_std(niter_in(:,ii)),...
          format_num_std(nsvd(:,ii)),...
          format_num_std(rank(:,ii)),...
          format_num_std(error(:,ii)));
  
end
