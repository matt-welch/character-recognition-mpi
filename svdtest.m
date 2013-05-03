m=8;n=m;
A = ones(m,1);
B = ceil(rand(m,n)*10); 
B_mn = mean(B,2);
C = B*A ./ m
B_mn - C