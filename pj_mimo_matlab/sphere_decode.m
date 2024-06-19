function hatx = sphere_decode(y, H, modem, x)
% Sphere_decoding	Sphere decoding
% Inputs	y -- received column vector;
%		H -- channel matrix; 
% 		cons -- the constellation of the PAM code
% 		which is in the following format:
%		L, L+2, ..., -3, -1, 1, 3, ..., U-2, U;
% 		modem -- the mode for the modulation.
%  Outputs	x -- grided (decoded) output;

% To follow the sphere packing convention, we 
% need to convert the relationship from
% y = H*x + n to r = M*u + v, where
% r in R^{2Rx X 1} is the received real vector
% (r = [Re(y^T) Im(y^T)]^T),
% u in Z^{2Tx X 1} is the transmitted real codes
% (u = [Re(x^T) Im(x^T)]^T),
% M in R^{2Tx X 2Rx} is the generator matrix
% (M = [Re(H) -Im(H); Im(H) Re(H)]),
% and v is the noise.
% (Note that we only consider the case of Tx <= Rx.)
%
% The output hatx is still in the usuall form, i.e.,
% hatx in C^{Tx X 1}, with contents being the elements
% in the constellation.

% Author: Yi Jiang, Univ. of Colorado Boulder, 3/8/2006

% Convert the data first
if modem == 2
   y = y*sqrt(2);
   constel = [-1 1];
elseif modem == 4
   y = y*sqrt(10);
   constel = [-3 -1 1 3];
elseif modem == 6
   y = y*sqrt(42);
   constel = [-7 -5 -3 -1 1 3 5 7];
else
   error('Wrong demodulation mode.');
end

% Change y in order to change the constelation.
% A line of QAM 64 is changed to -4 -3 -2 -1 0 1 2 3 .
q = max(constel);
y = (y+H*ones(size(H,2),1)*q*(1+1j))/2; % convert the data such that the constel is [0,1,...q]
% Form the needed vector and matrix first
y = [real(y); imag(y)];
M = [real(H) -imag(H); imag(H) real(H)];
[Q,R] = qr(M);
r = diag(R);
R = diag(sign(r))*R;
y = diag(sign(r))*Q'*y;
m = length(y);
%------------ The major story begins ...
dc = m*10;
T = 0;  Epsi = 0;
A = zeros(m,1); B = zeros(m,1);
A(m) = ceil((y(m)-Epsi-sqrt(dc-T))/R(m,m));
A(m) = max(0,A(m));
B(m) = floor((y(m)-Epsi+sqrt(dc-T))/R(m,m));
B(m) = min(q,B(m));
[X,d] = recur_sd(y,zeros(m,1),Epsi,dc,T,R,q,m,A,B);
[~,indx] = min(d);
hatx = X(:,indx);
hatx = hatx*2-q;

if modem == 2
   hatx = hatx/sqrt(2);
elseif modem == 4
   hatx = hatx/sqrt(10);
elseif modem == 6
   hatx = hatx/sqrt(42);
end
