function out = qpsk_modem(in, modem)
% QPSK modulation and demodulation
% Inputs	in -- input vector;
%		modem = 1 -- modulation, = others -- demodulation.
% Output	out -- output vector.

N = length(in);
if modem == 1          % modulation
    if mod(N,2)~=0
        error('Error: Input length should be a multiple of 2')
    end
    out = (in(1:2:N-1)*2-1) + 1j*(in(2:2:end)*2-1);
    out = out*sqrt(0.5);
elseif modem==0        % demodulation
    out = zeros(1,N*2);
    out(1:2:2*N-1) = (real(in)>0);
    out(2:2:2*N) = (imag(in)>0);
else
    error('modem should be 1: mod or 0: demod')
end