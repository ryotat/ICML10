% This is a dummy function for "flops" that only
% existed until MATLAB 5.
function out = flops(in)

global flops_counter;

if exist('in','var') && in==0
  flops_counter=0;
  out = 0;
else
  out = cputime-flops_counter;
end
