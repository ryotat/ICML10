function Z=funcXAYadj(aa,X,Y,I,J)

R=size(X,1);
C=size(Y,1);

A=sparse(I,J,aa,R,C);

Z=X*A*Y';


% $$$ R=size(X,2);
% $$$ C=size(Y,1);
% $$$ m=length(aa);
% $$$ 
% $$$ 
% $$$ Z=zeros(R,C);
% $$$ for ii=1:m
% $$$   Z=Z+X(I(ii),:)'*(aa(ii)*Y(:,J(ii))');
% $$$ end

function zz = funcXAYadjL(xx,aa,X,Y,I,J)
% zz=zeros(size(X,2),1);
% $$$ for ii=1:length(I)
% $$$   zz=zz+X(I(ii),:)'*(aa(ii)*(Y(:,J(ii))'*xx));
% $$$ end
zz = X(I,:)'*(spdiag(aa)*Y(:,J)'*xx);

function zz = funcXAYadjR(yy,aa,X,Y,I,J)
% $$$ zz=zeros(size(Y,1),1);
% $$$ for ii=1:length(I)
% $$$   zz=zz+Y(J(ii),:)'*(aa(ii)*(X(:,I(ii))'*yy));
% $$$ end
zz = Y(:,J)*(spdiag(aa)*X(I,:)*yy);
