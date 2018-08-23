function r = drchrnd(n)

a=ones(1,n);

p = length(a);
r = gamrnd(repmat(a,n,1),1,n,p);
r = r ./ repmat(sum(r,2),1,p);

r=r(1,:);
r=r';