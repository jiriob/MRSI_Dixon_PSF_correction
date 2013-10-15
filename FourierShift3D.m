function y = FourierShift2D(x, delta)
%
% y = FourierShift(x, [delta_x delta_y])
%
% Shifts x by delta cyclically. Uses the fourier shift theorem.
%
% Real inputs should give real outputs.
%
% By Tim Hutt, 26/03/2009
% Small fix thanks to Brian Krause, 11/02/2010

% The size of the matrix.
[N, M, O] = size(x);

%% FFT of our possibly padded input signal.
X = fftn(x);

%% The mathsy bit. The floors take care of odd-length signals.
x_shift = exp(-1i * 2 * pi * delta(1) * [0:floor(N/2)-1 floor(-N/2):-1]' / N);
y_shift = exp(-1i * 2 * pi * delta(2) * [0:floor(M/2)-1 floor(-M/2):-1] / M);
z_shift = exp(-1i * 2 * pi * delta(3) * [0:floor(O/2)-1 floor(-O/2):-1] / O);
z_shift = permute(z_shift,[1,3,2]);
%a=[0:floor(N/2)-1 floor(-N/2):-1];

%% Force conjugate symmetry. Otherwise this frequency component has no
% corresponding negative frequency to cancel out its imaginary part.
if mod(N, 2) == 0
	x_shift(N/2+1) = real(x_shift(N/2+1));
end
if mod(M, 2) == 0
	y_shift(M/2+1) = real(y_shift(M/2+1));
end
if mod(O, 2) == 0
	z_shift(O/2+1) = real(z_shift(O/2+1));
end
%% repmat everything
x_shift3d = repmat(x_shift,[1,size(y_shift,2),size(z_shift,3)]);
y_shift3d = repmat(y_shift,[size(x_shift,1),1,size(z_shift,3)]);
z_shift3d = repmat(z_shift,[size(x_shift,1),size(y_shift,2),1]);

%%
Y = X .* (x_shift3d .* y_shift3d .* z_shift3d);

%% Invert the FFT.
y = ifftn(Y);

% There should be no imaginary component (for real input
% signals) but due to numerical effects some remnants remain.
% if isreal(x)
%     y = real(y);
% end

end