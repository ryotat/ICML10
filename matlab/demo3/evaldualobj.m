function dval = evaldualobj(X, Y, ww, bb, lambda, blks)

[floss,gg]=loss_lrp(X*ww+bb, Y);

aa = -gg;

aa=aa-mean(aa);

dnm = ds_dnorm(X'*aa,blks);

aa  = min(1, lambda/dnm)*aa;

dval = -loss_lrd(aa, Y);


