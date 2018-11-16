function [Xjw_single_side, f_vec] = cal_FFT(X, fs)

Ts = 1/fs;
L  = length(X);

n = 2^nextpow2(L);
Xjw = fft(X,n);

Xjw_single_side = Xjw( 1 : n/2+1 );

f_vec = fs*( 0 : (n/2) )'/n;