function  c =  W3SCcw_ADMM( y, D, S, W1, W2, Par )
% This routine solves the following weighted nuclear norm optimization problem with column weights,
%
% min |W1(y-DSc)W2|_2,2 + |z|_1 s.t.  c=z
% inputs:
%        y -- d*1 data matrix, d is the data dimension, and M is the number
%             of image patches.
%        W1 -- d*d matrix of row weights
%        W2 -- 1 weight

% tol = 1e-10;
% Initializing optimization variables
c = zeros(size(D, 2), 1);
z = zeros(size(c));
u = zeros(size(c));
% Start main loop
iter = 0;
% stopCZ = zeros(Par.maxIter, 1);
% stopC = zeros(Par.maxIter, 1);
% stopZ = zeros(Par.maxIter, 1);
while iter < Par.maxIter
    iter = iter + 1;
%     cpre = c;
%     zpre = z;
    %% update C, fix Z and U
    % min_{c} W2^2 * ||W1 * (y - DSc)||_2^2 + 0.5 * rho * ||c - z + 1/rho * u||_2^2
    % The solution is equal to solve A * X + X * B = E
    Temp1 = S' * D' * diag(W1.^2);
    A = Temp1 * D * S;
    W2inv = (1/W2)^2;
    B = 0.5 * Par.rho * W2inv * eye(size(S, 2));
    E = Temp1 * y + 0.5 * (Par.rho * z - u) * W2inv;
    c = (A + B) \ E;
    
    %% update Z, fix X and D
    % min_{z} 0.5 * rho * ||z - (c + 1/rho * u)||_2^2 + ||z||_1
    Temp2 = c + u/Par.rho;
    z = sign(Temp2) .* max( abs(Temp2) - 1/Par.rho, 0 );
    
%     %% check the convergence conditions
%     stopCZ(iter) = max(max(abs(c - z)));
%     stopC(iter) = max(max(abs(c - cpre)));
%     stopZ(iter) = max(max(abs(z - zpre)));
%     if Par.display %&& (iter==1 || mod(iter,10)==0 || stopC<tol)
%         disp(['iter ' num2str(iter) ', mu=' num2str(Par.mu,'%2.1e') ...
%             ', max(||c-z||)=' num2str(stopCZ(iter),'%2.3e') ...
%             ', max(||c-cpre||)=' num2str(stopC(iter),'%2.3e') ...
%             ', max(||z-zpre||)=' num2str(stopZ(iter),'%2.3e')]);
%     end
%     if stopCZ(iter) < tol && stopC(iter) < tol && stopZ(iter) < tol
%         break;
%     else
        %% update the augmented multiplier D, fix Z and X
        u = u + Par.rho * (c - z);
        Par.rho = min(Par.maxrho, Par.mu * Par.rho);
%     end
end
return;
