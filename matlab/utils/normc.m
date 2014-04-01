function nm=normc(W,mode,varargin)

opt=propertylist2struct(varargin{:});
opt=set_defaults(opt, 'v',[],'div', 0, 'dosum', 1);

if isstruct(W)
  W=W.W;
end

v = opt.v;
if isempty(v)
  v=ones(1,length(W));
end
orig_v = v;


if opt.div==1
  v=v/v(1);
end

nm=[];
for ii=1:length(W);
  nm = [nm,svd(W{ii})'*sqrt(v(ii)/orig_v(ii))];
end

if opt. dosum
  switch(mode)
   case 'ds'
    nm = sum(nm);
   case 'fro'
    nm = sqrt(sum(nm.^2));
   otherwise 
    error('mode: %s not supported.',mode);
  end
end
