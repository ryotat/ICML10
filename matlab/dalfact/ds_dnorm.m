% ds_dnorm - conjugate of the dual spectral norm regularizer
%
% Copyright(c) 2009 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt
function [nm,ishard]=ds_dnorm(ww,blks)

nm=0;
for kk=1:size(blks,1)
  nm=max(nm,lansvd(ww,1));
end
ishard=1;