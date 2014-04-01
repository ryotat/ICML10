% objdall1 - objective function of DAL with L1 regularization
%
% Copyright(c) 2009 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt

function varargout=objdall1fminunc(aa, prob, ww, uu, A, B, lambda, eta, tol)

m = length(aa);
n = length(ww);
vv = A'*aa+ww/eta;

Ip = find(vv>lambda);
In = find(vv<-lambda);

if nargout<=2
  [floss, gloss, hmin]=feval(prob.floss.d,aa, prob.floss.args{:});
else
  [floss, gloss, hloss, hmin]=prob.floss.d(aa, prob.floss.args{:});
end


vsth = l1_softth(vv,lambda);


fval = floss+0.5*eta*sum(vsth.^2);
if ~isempty(uu)
  u1   = uu/eta+B'*aa;
  fval = fval + 0.5*eta*sum(u1.^2);
end

varargout{1}=fval;

if nargout>=2
  gg  = gloss+eta*(A*vsth);
  soc = sum((vsth-ww/eta).^2);
  if ~isempty(uu)
    gg  = gg+eta*(B*u1);
    soc = soc+sum((B'*aa).^2);
  end

  if soc>0
    ginfo = norm(gg)/(sqrt(eta*hmin*soc));
  else
    ginfo = inf;
  end
  
  if ginfo>tol
    varargout{2} = gg;
  else
    varargout{2} = zeros(size(gg));
  end

  if nargout==3
    I = sort([Ip; In]);
    AF = A(:,I);
    if length(I)>0
      varargout{3} = hloss+eta*AF*AF';
    else
      varargout{3} = hloss;
    end
    if ~isempty(uu)
      varargout{3} = varargout{3}+eta*B*B';
    end
  end
end

