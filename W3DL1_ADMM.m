function  [C] =  W3DL1_ADMM( Y, D, S, W1, W2, Par )
% This routine solves the following weighted nuclear norm optimization problem with column weights,
%
% min |W1(Y-DSC)W2|_1 + |Z|_1 s.t.  C=Z
% inputs:
%        Y -- d*M data matrix, d is the data dimension, and M is the number
%             of image patches.
%        W1 -- d*d matrix of column weights
%        W2 -- M*M matrix of column weights

tol = 1e-8;
Par.maxrho = 100;
Par.rho = 0.5;
Par.mu = 1.1;
Par.display = 1;
Par.maxIter = 100;
% Initializing optimization variables
C = zeros(size(S, 1), size(Y, 2));
A = zeros(size(C));
E = zeros(size(Y));
Delta = zeros(size(Y));
delta = zeros(size(C));
% Start main loop
iter = 0;
stopCA = zeros(Par.maxIter, 1);
stopC = zeros(Par.maxIter, 1);
stopA = zeros(Par.maxIter, 1);
stopEW = zeros(Par.maxIter, 1);
stopE = zeros(Par.maxIter, 1);
while iter < Par.maxIter
    iter = iter + 1;
    Cpre = C;
    Apre = A;
    Epre = E;
    %% update A, fix C and E
    % min_{C} ||W1 * (Y - DSC) * W2||_F^2 + 0.5 * rho * ||C - Z + 1/rho * U||_F^2
    % The solution is equal to solve A * X + X * B = E
    AA = S' * D' * diag(diag(W1).^2) * D * S;
    BB = diag(1./(W2.^2));
    CC = (S' * D' * W1 * (W1 * Y * diag(W2) - E - 1/Par.rho * Delta) * diag(W2) + C - 1/Par.rho * delta) * BB;
    A = sylvester(AA, BB, CC);
    
    %     % faster solution
    %     [Ua, Sa, ~] = svd(A);
    %     I1 = eye(size(A, 2));
    %     I2 = eye(size(B, 1));
    %     K = kron(I1, A) + kron(B', I2);
    %     invK = 1./diag(K);
    %
    %     UTE = Ua'*E;
    %     vecUTE = UTE(:);
    %
    %     vecUTC = invK .* vecUTE;
    %     MatvecUTC = reshape(vecUTC, [size(UTE, 1) size(UTE, 2)]);
    %     C = Ua*MatvecUTC;
    
    %% update C, fix A and E
    % min_{Z} 0.5 * rho * ||Z - (C + 1/rho * U)||_F^2 + ||Z||_1
    Temp1 = A + delta/Par.rho;
    C = sign(Temp1) .* max( abs(Temp1) - 1/Par.rho, 0 );
    
    %% update E, fix A and C
    % min_{Z} 0.5 * rho * ||Z - (C + 1/rho * U)||_F^2 + ||Z||_1
    Temp2 = W1 * (Y - D * S * A) * diag(W2);
    Temp3 = Temp2 - Delta/Par.rho;
    E = sign(Temp3) .* max( abs(Temp3) - 1/Par.rho, 0 );
    
    %% check the convergence conditions
    stopCA(iter) = max(max(abs(C - A)));
    stopA(iter) = max(max(abs(A - Apre)));
    stopC(iter) = max(max(abs(C - Cpre)));
    stopEW(iter) = max(max(abs(E - Temp2)));
    stopE(iter) = max(max(abs(E - Epre)));
    if Par.display %&& (iter==1 || mod(iter,10)==0 || stopC<tol)
        disp(['iter ' num2str(iter) ', mu=' num2str(Par.mu,'%2.1e') ...
            ', max(||A-C||)=' num2str(stopCA(iter),'%2.3e') ...
            ', max(||A-Apre||)=' num2str(stopA(iter),'%2.3e') ...
            ', max(||C-Cpre||)=' num2str(stopC(iter),'%2.3e') ...
            ', max(||E-W1(Y - DSA)W2||)=' num2str(stopEW(iter),'%2.3e') ...
            ', max(||E-Epre||)=' num2str(stopE(iter),'%2.3e')]);
    end
    if stopCA(iter) < tol && stopC(iter) < tol && stopA(iter) < tol
        break;
    else
        %% update the augmented multipliers Delta, delta
        Delta = Delta + Par.rho * (E - Temp2);
        delta = delta + Par.rho * (A - C);
        Par.rho = min(Par.maxrho, Par.mu * Par.rho);
    end
end
return;
