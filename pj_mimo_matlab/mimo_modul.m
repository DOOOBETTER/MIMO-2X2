function [data] = mimo_modul(bits, modem, numOfStream)

% MIMO Modulation for the non-ofdm communications.
% Data modulation will be BPSK, if modem == 1,
%                         QPSK, if modem == 2,
%                        QAM16, if modem == 4,
%                        QAM64, if modem == 6.
% The output data in numOfStream rows

if modem == 2
   temp = qpsk_modem(bits, 1);
   data = temp/sqrt(2);
elseif modem == 4
   temp = qam16(bits, 1);
   data = temp/sqrt(10);
elseif modem == 6
   temp = qam64(bits, 1);
   data = temp/sqrt(42);
elseif modem == 8
    temp = qam256(bits, 1);
    data = temp/sqrt(170);
else
   error('Error: Unknown modulation mode.')
end

data = reshape(data,numOfStream,[]);
return;