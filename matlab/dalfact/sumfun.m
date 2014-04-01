% fun = sumfun(a, fun1, b, fun2)
function fun = sumfun(a, fun1, b, fun2)

fun=cell(size(fun1));
for kk=1:size(fun1,1)
  fun{kk,1}=@(x)a*fun1{kk,1}(x)+b*fun2{kk,1}(x);
  fun{kk,2}=@(y)a*fun1{kk,2}(y)+b*fun2{kk,2}(y);
  fun{kk,3}=fun1{kk,3};
  fun{kk,4}=fun1{kk,4};
  fun{kk,5}=fun1{kk,5}+fun2{kk,5};
end
