% ds_softth - soft threshold function for DS regularization
%
% Copyright(c) 2009 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt

function [vv,ss,info]=ds_softth(vv,lambda,info)

R=size(vv.U,1);
C=size(vv.V,1);
ss=min(R,C);

ixs=0;

for kk=1:size(info.blks,1)
  blk=[R,C];
  J=ixs+(1:min(blk));
  ixs=J(end);
  vcell={@(x)multlr(x,vv),@(x)multlrt(x,vv),R,C};
  if strcmp(info.solver,'cg')
    [U,S,V]=lansvd(vcell{:},min(blk));
  else
    [U,S,V]=svdmaj(vcell(:), lambda, 'kinit', round(info.nsv(kk)*1.05));
  end
  dd=diag(S);
  K=find(dd>lambda);
  ssk=dd(K)-lambda;
  ss(J(1:length(ssk)))=ssk;
  %  vvk=U(:,K)*diag(ssk)*V(:,K)';
  if length(K)>0
    vv(kk,1)=struct('U',U(:,K),...
                    'ss',ssk,...
                    'V',V(:,K),...
                    'D',[]);
  else
    vv(kk,1)=struct('U',zeros(R,1),'ss',0,'V',zeros(C,1),'D',[]);
  end
  info.nsv(kk)=length(ssk);
  
  info.U{kk}=U;
  info.S{kk}=S;
  info.V{kk}=V;

end
