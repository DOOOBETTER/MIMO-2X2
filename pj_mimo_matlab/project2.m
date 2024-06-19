% Compare performance of UCD-VBLAST, UCD-DP, sphere decoder,
clc,clear
close all
numMonteCarlo = 100;
nR = 2;     nT = 2;    % channel dimensionality
K = min(nR,nT);
% H = randn(nR,nT)+1j*randn(nR,nT);
% H = [1 -1j;exp(1j*pi/4) exp(-1j*pi/3)];
numBits = 960; % number of bits for each channel realizations.
modem = 2;
% Data modulation will be BPSK, if modem == 1,
%                         QPSK, if modem == 2,
%                        QAM16, if modem == 4,
%                        QAM64, if modem == 6.

M = 2^(modem/2); % PAM constellation size
if modem == 2   % constellation point spacing
    d = 2/sqrt(2);
elseif modem == 4
    d = 2/sqrt(10);
elseif modem == 6
    d = 2/sqrt(42);
end
snrDbSet = 0:3:30;
berMmse = zeros(1,length(snrDbSet)); % MMSE equalizer
berZf = zeros(1,length(snrDbSet)); % ZF equalizer
berSd = zeros(1,length(snrDbSet)); % sphere decoder
berUcdVb = zeros(1,length(snrDbSet)); % UCD-VBLAST
berAlamouti = zeros(1,length(snrDbSet)); % Alamouti code & ML receiver

for pp = 1:numMonteCarlo
    H = sqrt(0.5)*(rand(nR,nT)+1j*randn(nR,nT));%[1 -1j;exp(1j*pi/4) exp(1j*pi/3)];
    bits = round(rand(1,numBits));
    temp = qpsk_modem(bits, 1); %modulate 1/0 bits into QPSK symbols
    qpskData = reshape(temp/sqrt(2),K,[]);
    temp = qam16_modem(bits, 1); %modulate 1/0 bits into QPSK symbols
    % qam16Data = qammod(bits.',16,'InputType','bit').'; % modulate bits into QAM16 symbols
    qam16Data = reshape(temp/sqrt(2),1,[]);
    [U, S, V] = svd(H);
    ns = randn(K,size(qpskData,2))+1j*randn(K,size(qpskData,2));
    for nn = 1:length(snrDbSet)
        rho = db2pow(snrDbSet(nn)); % SNR
        sigma2 = K/rho;
        ns1 = ns*sqrt(.5*K/rho);
        decData = zeros(K,numBits/K/modem);
        %----------------------- Sphere Decoding uncoded  QPSKx2
        y = H*qpskData+ns1;
        decoded = zeros(size(y));
        for ii = 1:size(y,2)
            hatx = sphere_decode(y(:,ii), H, modem);
            decoded(:,ii) = hatx(1:nT)+1j*hatx(1+nT:end);
        end
        for k = K:-1:1
            decBits(k,:) = qpsk_modem(decoded(k,:), 0);
        end
        decBit = zeros(size(bits));
        for n = 1:size(decBits,2)/modem
            decBit((n-1)*modem*K+(1:modem*K)) = reshape(decBits(:,(n-1)*modem+(1:modem)).',1,[]);
        end
        berSd(nn) = berSd(nn)+sum(bits ~= decBit)/numBits/numMonteCarlo;

        %------------------------ UCD-VBLAST uncoded QPSKx2
        decBits = zeros(K,numBits/K);
        [P,W] = ucd(U,S,V,K,sigma2,1);
        Htil = H*P;
        y = Htil*qpskData+ns1;
        z = zeros(K,size(y,2));
        for k = K:-1:1
            if k < K
                y = y-Htil(:,k+1)*decData(k+1,:);
            end
            z(k,:) = W(:,k)'*y;
            tmp = qpsk_modem(z(k,:), 0);
            decBits(k,:) = tmp.';
            decData(k,:) = qpsk_modem(tmp,1);
        end
        decBit = zeros(size(bits));
        for n = 1:size(decBits,2)/modem
            decBit((n-1)*modem*K+(1:modem*K)) = reshape(decBits(:,(n-1)*modem+(1:modem)).',1,[]);
        end
        berUcdVb(nn) = berUcdVb(nn)+sum(bits ~= decBit)/numBits/numMonteCarlo;

        %----------------------- ZF 均衡 QPSKx2 uncoded
        %  30pt
        y = H*qpskData+ns1;
        W_ZF = inv(conj(H)*H)*conj(H);
        % H_pinv = pinv(H); % 计算伪逆
        x_hat = W_ZF * y; % ZF 均衡后的信号
        for k = K:-1:1
            decBits(k,:) = qpsk_modem(x_hat(k,:), 0);
        end
        decBit = zeros(size(bits));
        for n = 1:size(decBits,2)/modem
            decBit((n-1)*modem*K+(1:modem*K)) = reshape(decBits(:,(n-1)*modem+(1:modem)).',1,[]);
        end
        berZf(nn) = berZf(nn) + sum(bits ~= decBit) / numBits / numMonteCarlo;
        %----------------------- MMSE 均衡 QPSKx2 uncoded
        %  20pt
        y = H*qpskData+ns1;
        W_MMSE = inv(conj(H)*H+sigma2*eye(nR,nT))*conj(H);
        x_hat = W_MMSE * y;% MMSE 均衡后的信号
        for k = K:-1:1
            decBits(k,:) = qpsk_modem(x_hat(k,:), 0);
        end
        decBit = zeros(size(bits));
        for n = 1:size(decBits,2)/modem
            decBit((n-1)*modem*K+(1:modem*K)) = reshape(decBits(:,(n-1)*modem+(1:modem)).',1,[]);
        end
        berMmse(nn) = berMmse(nn) + sum(bits ~= decBit) / numBits / numMonteCarlo;
        %----------------------- Alamouti 16QAMx1 coded rate 1/2
        %  30pt
        alamoutiData = alamouti_code(qam16Data); % alamouti 编码 15pt
        tmp = zeros(size(qam16Data));
        h11 = H(1, 1);
        h12 = H(1, 2);
        h21 = H(2, 1);
        h22 = H(2, 2);
        normH = sum(sum(abs(H).^2));
        for i = 1:2:size(alamoutiData,2)
            n11 = ns1( 1 , i );
            n12 = ns1( 1 , i+1 );  
            n21 = ns1( 2 , i );
            n22 = ns1( 2 , i+1 );
            s1 = alamoutiData( 1 , i );
            s2 = alamoutiData( 2 , i );
            y11 = s1 * h11 + s2 * h12 + n11;
            y12 = -conj(s2) * h11 + conj(s1) * h12 + n12;
            y21 = s1 * h21 + s2 * h22 + n21;
            y22 = -conj(s2) * h21 + conj(s1) * h22 + n22;
            % alamouti 解码  15pt
            s_hat1 = (conj(h11) * y11 + h12 * conj(y12) + conj(h21) * y21 + h22*conj(y22))/normH;
            s_hat2 = (conj(h12) * y11 - h11 * conj(y12) + conj(h22) * y21 - h21*conj(y22))/normH;           
            tmp(i) = s_hat1;
            tmp(i+1)= s_hat2;
        end
        demulatedBits = qam16_modem(tmp, 0);
        decBit = demulatedBits';
        berAlamouti(nn) = berAlamouti(nn) + sum(bits ~= decBit) / numBits / numMonteCarlo;
    end

    if ~mod(pp,10)
        fprintf('%d \n',pp)
        %         berSd*numMonteC/pp
        %         berUcdVb*numMonteC/pp
    end
end
%%
semilogy(snrDbSet,berUcdVb,'r-^',snrDbSet,berZf,'g-^',snrDbSet,berMmse,'b-s',snrDbSet,berAlamouti,'m-*',snrDbSet,berSd,'k');
xlabel('SNR (dB)')
ylabel('BER')
legend('UCD-VB', 'ZF','MMSE','Alamouti','Sphere Decoder')
return
