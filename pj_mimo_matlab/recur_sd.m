function [X,D,A,B] = recur_sd(y,x,Epsi,dc,T,R,q,i,A,B);
% Apply Sphere Decoding in a recursive way.
% Algorithm I in Damen et. al. IT 2003 October paper

% Author: Yi Jiang, Univ. of Colorado Boulder, 3/9/2006
m = length(y);
X = [];
D = [];
if dc < T % step 4
    i = min(i+1,m);
    return;
end
A(i) = ceil((y(i)-Epsi-sqrt(dc-T))/R(i,i));
A(i) = max(0,A(i));
B(i) = floor((y(i)-Epsi+sqrt(dc-T))/R(i,i));
B(i) = min(q,B(i));
for t = A(i):B(i) % step 3
    x(i) = t;
    if i>1 % step 5
        Epsi1 = dot(R(i-1,i:end),x(i:end));
        
        T1 = T + (y(i)-Epsi-R(i,i)*x(i))^2; 
        [xtmp,d,A,B] = recur_sd(y,x,Epsi1,dc,T1,R,q,i-1,A,B);
        X = [X,xtmp];
        D = [D,d];
    else % step 6
        dhat = T+(y(1)-Epsi-R(1,1)*x(1))^2;
        if dhat < dc  % a valid point is found
            dc = dhat;
            D = [D,dhat];
            X = [X,x];
            for j = 1:m
                B(j) = floor((y(j)-Epsi+sqrt(dc-T))/R(j,j));
                B(j) = min(q,B(j));
            end
        end
    end
end
i = min(i+1,m); % step 4