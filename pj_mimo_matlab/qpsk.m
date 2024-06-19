function data = qpsk(bits, flag)
% modulate bits into QPSK symbols or
% demodulate qpsk symbols into bits
% symbol is normalize to be with engery 1;
% Author Yi Jiang (SAL,UF) Jan. 17, 2004
N = length(bits);
if flag == 1  %modulate
    data = zeros(N/2,1);
   for m = 1:N/2
      tmp = bits(m*2-1:m*2);      if tmp == [0; 0]
         data(m) = 1+j;
      elseif tmp == [0; 1]
         data(m) = 1-j;
      elseif tmp ==[1; 1]
         data(m) = -1-j;
      elseif tmp ==[1; 0]
         data(m) = -1+j;
      else
         disp('error, input bits should be 0 or 1');
      end
   end
%   data = data/sqrt(2);
else %demodulate
   data = zeros(2*N,1);
   symbs = bits;
   data(1:2:2*length(symbs)-1) = sign(real(symbs));
   data(2:2:2*length(symbs)) = sign(imag(symbs));
   data = (1-data)/2;
end