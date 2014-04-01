function zz=funcXAY(A,X,Y,ind)

Xu=A.U'*X;
Yv=A.V'*Y;

Z =Xu'*diag(A.ss)*Yv;

zz=Z(ind);

% zz=zeros(n,1);
% for ii=1:n
%  zz(ii)=sum(Xu(:,I(ii)).*A.ss.*Yv(:,J(ii)));
  
% $$$   if ~isempty(A.D)
% $$$     zz(ii)=zz(ii)+X(I(ii),:)*(A.D*Y(:,J(ii)));
% $$$   end
%end

