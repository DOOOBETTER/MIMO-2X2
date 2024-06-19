function [P,W,rhobar] = ucd(U,S,V,L,alpha,waterfill)
% apply Even Channel Decomposition algorithm
% Input: U,S,V: the svd of the Channel matrix H
%        L: the # of subchannels
%        alpha = sigma2_n/sigma2_x = sigma2_n = L/rho (L/SNR)
% Output: P, precoder unitary matrix
%         W, receiver weight
%         rhobar, SNR of each subchannel (the subchannels are the same)
% Yi Jiang, July 28, 2004

if ( nargin < 6 )
    waterfill = 1 ;
end
if size(V,2) > L
    V = V(:,1:L);
end
K = size(V,2);
sigma2_n = alpha;
% Step 1: waterfilling...
if waterfill
    sv2 = diag(S).^2;
    for k = K:-1:1
        mu = L/k+sum(sigma2_n./sv2(1:k))/k;
        if mu > sigma2_n/sv2(k)
            break
        end      
    end
    ldpow = zeros(L,1);
    
    ldpow(1:k) = mu-sigma2_n./sv2(1:k);
else
    ldpow = ones(L,1);
end
Sigm = zeros(L,1);
Sigm(1:K) = sqrt(ldpow(1:K)).*diag(S);
SigmAlpha = sqrt(Sigm.^2+alpha);
[W,R,P] = gmd([U(:,1:K)*diag(Sigm(1:K)./SigmAlpha(1:K)),zeros(size(U,1),L-K)], diag(SigmAlpha),[V*diag(sqrt(ldpow(1:K))),zeros(size(V,1),L-K)]);
W = W/diag(diag(R));
rhobar = R(1,1)^2/alpha-1;