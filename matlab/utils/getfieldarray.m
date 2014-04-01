function C=getfieldarray(A, field)

C=cell(size(A));

for i=1:prod(size(A))
  C{i}=getfield(A(i), field);
end
