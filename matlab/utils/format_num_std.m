function str = format_num_std(val)

mval=mean(val);
sval=std(val);

str= [format_num_fixed(mval,3), '(', char(177), ' ',...
       format_num_fixed(sval,3) ') '];




