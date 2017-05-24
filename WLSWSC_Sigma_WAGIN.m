function  [im_out, par] = WLSWSC_Sigma_WAGIN(par)
im_out    =   par.pim;
% parameters for noisy image
[h,  w, ch]      =  size(im_out);
par.h = h;
par.w = w;
par.ch = ch;
par = SearchNeighborIndex( par );
% original noisy image to patches
NY = Image2PatchNew( par.nim, par );
PY = Image2PatchNew( par.pim, par );
for ite  =  1 : par.outerIter
    % iterative regularization
    im_out = im_out + par.delta * (par.pim - im_out);
    % image to patches
    Y = Image2PatchNew( im_out, par );
    % estimate local noise variance, par.lambdals is put here since the MAP
    % and Bayesian rules
    Sigma = par.lambda2 * sqrt(abs(repmat(par.nSig^2, 1, size(Y,2)) - mean((PY - Y).^2))); %Estimated Local Noise Level
    % estimation of noise variance
    if mod(ite-1, par.innerIter)==0
        par.nlsp = par.nlsp - par.nlspgap;
        % searching  non-local patches
        blk_arr = Block_Matching( Y, par );
        if ite == 1 && par.lambda2~=0
            Sigma = par.nSig * ones(size(Sigma));
        end
    end
    % Weighted Sparse Coding
    Y_hat = zeros(par.ps2ch, par.maxrc, 'double');
    W_hat = zeros(par.ps2ch, par.maxrc, 'double');
    for i = 1:par.lenrc
        index = blk_arr(:, i);
        nlY = Y( : , index );
        DC = mean(nlY, 2);
        nDCnlY = bsxfun(@minus, nlY, DC);
        W2 = 1 ./ (Sigma(index) + eps);
        
%        % Compute W1
%        nlNY = NY( : , index );
%        residual  =    mean((nlNY-nlY).^2, 2);
%        W1         =   1./exp(par.lambda1*residual);

        % update D
        [D, S, ~] = svd( full(nDCnlY), 'econ' );
        % update S
        S = sqrt(max( diag(S).^2 - length(index) * Sigma(index(1))^2, 0 ));
%         S = diag(S);
%         % Solve Sylvester equation AX + XB = E for X
%         A = S' * D' * diag(W1.^2) * D * S;
%         B = diag(1./(W2.^2));
%         E = S' * D' * diag(W1.^2) * nlY;
%         C = sylvester(A, B, E);
        % update weight for sparse coding
        Wsc = bsxfun( @rdivide, Sigma(index).^2, S + eps );
        % update C by soft thresholding
        B = D' * nDCnlY;
        C = sign(B) .* max( abs(B) - Wsc, 0 );
        % update Y
        nDCnlYhat = D * C;
        % add back DC components
        nlYhat = bsxfun(@plus, nDCnlYhat, DC);
        % aggregation
        Y_hat(:, index) = Y_hat(:, index) + bsxfun(@times, nlYhat, W2);
        W_hat(:, index) = W_hat(:, index) + repmat(W2, [par.ps2ch, 1]);
%         % aggregation
%         Y_hat(:, index) = Y_hat(:, index) + bsxfun(@times, bsxfun(@times, W1, nlYhat), W2);
%         W_hat(:, index) = W_hat(:, index) + W1 * W2;
    end
    % Reconstruction
    im_out = PGs2Image(Y_hat, W_hat, par);
    % calculate the PSNR and SSIM
    PSNR =   csnr( im_out , par.I , 0, 0 );
    SSIM      =  cal_ssim( im_out , par.I , 0, 0 );
    fprintf('Iter %d : PSNR = %2.4f, SSIM = %2.4f\n', ite, PSNR, SSIM);
    par.PSNR(ite, par.image) = PSNR;
    par.SSIM(ite, par.image) = SSIM;
end
return;

