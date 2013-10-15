%% hamming filter first:
% generate 1D filter with CSI matrix resolution:
w1 = hamming(11);
[x,y,z] = meshgrid(-5:1:5);
r = sqrt(x.^2 + y.^2 + z.^2);
w = zeros(size(r));
w(r<=5) = interp1(linspace(-5,5,11),w1,r(r<=5));
%imagesc(w(:,:,6));
% fill the matrix to 12x12x12
w = padarray(w,[1 1 1],'post');
sgmnts{6,1} = sgmnts{5,1} .* w; % multiply the hamming filter with fourier transform
% zero filling to 16x16x16
sgmnts{7,1} = abs(ifftn(padarray(sgmnts{6,1},[2 2 2])));







sgmnts{6,1} = abs(ifftn(padarray(sgmnts{5,1},[2 2 2])));
%sgmnts{7,2} = abs(ifftn(padarray(sgmnts{6,2},[2 2 2])));
%% hamming filter:
% generate 1D filter with CSI matrix resolution:
w1 = hamming(15);
[x,y,z] = meshgrid(-7:1:7);
r = sqrt(x.^2 + y.^2 + z.^2);
w = zeros(size(r));
w(r<=7) = interp1(linspace(-7,7,15),w1,r(r<=7));
%imagesc(w(:,:,6));
% fill the matrix to 12x12x12
w = padarray(w,[1 1 1],'post');
sgmnts{7,1} = sgmnts{6,1} .* w; % multiply the hamming filter with fourier transform

