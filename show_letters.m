% calculate mean vector of reference set */
% subtract the mean from the test image and the reference set (mean center) */
% perform svd on the normalized A to get basis matrices, U & V */
% pdgesvd_() */
% Multiply X=U'A  and  x=U'T (corrected from x=UU'T) */
% Find mimimum:
% Subtract x from each column of X to make D */
% caculate sqrt(D'D) */
% find the minimum diagonal element, this is the matching character */
% print number of matching character */
clear all
close all
m=72;
n=128*128; % 16384

fidr = fopen('referenceset.bin');
A = fread(fidr, m*n,'double');
A = reshape(A, n, m);
for(i=1:m)
    As(:,:,i) = reshape(A(:,i),128,128);
%    imshow(As(:,:,i));
%     pause(0.25);
end

fidr = fopen('testcharacter.bin');
T = fread(fidr, n,'double');
Ts = reshape(T,128,128);
imshow(Ts(:,:));

% now we do the SVD

A_mean = mean(A,2);
T_mc = T - A_mean;
A_mean = repmat(A_mean,1,72);
A_mc = A - A_mean;

imshow(reshape(T_mc,128,128));

% need a waitbar here - this step takes forever
[U,S,V]=svd(A_mc);

X = U' * A_mc;
x = U' * T_mc;

x_mat = repmat(x,1,72);    
D = X - x_mat;

final = D'*D;
values = zeros(1,m);
for(i=1:m)
    values(i)=final(i,i);
end

minval = find(values == min(values));
figure; imshow(Ts); title('Testchar');
figure; imshow(As(:,:,minval)); title('refchar');
