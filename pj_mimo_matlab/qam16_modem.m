function out = qam16_modem(in, modem)
%Qam16		16QAM modulation and demodulation (hard output) for IEEE802.11.
% Inputs	in -- column vector;
%		modem = 1 -- modulation, = others -- demodulation.
% Output	out -- column vector.


N = length(in);

if modem == 1          % do modulation
   if ~(rem(N,4)==0)
      error('Error: Input length should be a multiple of 4')
   end
   
   NN = N/4;
   out = zeros(NN,1);
   
   for i = 1 : NN
      bits = in(4*i-3:4*i);
      if        bits(1:2) == [0  0], Inphase = -3;
         elseif bits(1:2) == [0 1], Inphase = -1;
         elseif bits(1:2) == [1 1], Inphase =  1;
         elseif bits(1:2) == [1 0], Inphase =  3;
      else                        
          error('Wrong data bits');
      end

      if        bits(3:4) == [0 0], Quad = -3;
         elseif bits(3:4) == [0 1], Quad = -1;
         elseif bits(3:4) == [1 1], Quad =  1;
         elseif bits(3:4) == [1 0], Quad =  3;
      else                        
          error('Wrong data bits');
      end
      
      out(i) = Inphase + 1j*Quad;
   end
   
   return;
else                   % do demodulation
   out = zeros(N*4,1);
   
   for k = 1 : N
      sym = in(k);
      Rsym = real(sym);
      Isym = imag(sym);
      
      if        Rsym <  -2,            out(k*4-3:k*4-2) = [0; 0];
         elseif Rsym >= -2 && Rsym < 0, out(k*4-3:k*4-2) = [0; 1];
         elseif Rsym >=  0 && Rsym < 2, out(k*4-3:k*4-2) = [1; 1];
         else                          out(k*4-3:k*4-2) = [1; 0];
      end
      
      if        Isym <  -2,            out(k*4-1:k*4) = [0; 0];
         elseif Isym >= -2 && Isym < 0, out(k*4-1:k*4) = [0; 1];
         elseif Isym >=  0 && Isym < 2, out(k*4-1:k*4) = [1; 1];
         else                          out(k*4-1:k*4) = [1; 0];
      end
   end

   return;
end
