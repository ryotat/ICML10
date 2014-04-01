function str=format_num_fixed(val,width)

p=floor(log10(val));

if val<1
  fmt=sprintf('%%.%dg',width);
  str=sprintf(fmt,val);
else
  if val==round(val)
    str=sprintf('%d',val);
  else
    fmt=sprintf('%%.%df',max(0,width-p-1));
    str=sprintf(fmt,val);
  end
end


