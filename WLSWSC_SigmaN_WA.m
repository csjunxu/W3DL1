function  [im_out, par] = WLSWSC_SigmaN_WA(par)
im_out    =   par.nim;
% parameters for noisy image
[h,  w, ch]      =  size(im_out);
par.h = h;
par.w = w;
par.ch = ch;
par = SearchNeighborIndex( par );
% % original noisy image to patches
% NY = Image2PatchNew( par.nim, par );
% % noise and weights for each patch
% Sigma = ones(1, par.maxrc);
% Wls = ones(1, par.maxrc);
for ite  =  1 : par.outerIter
    % iterative regularization
    im_out = im_out + par.delta * (par.nim - im_out);
    % image to patches
    Y = Image2PatchNew( im_out, par );
    %     % estimate local noise variance, par.lambdals is put here since the MAP
    %     % and Bayesian rules
    %     Sigma = par.lambdals * sqrt(abs(repmat(par.nSig0^2, 1, size(Y,2)) - mean((NY - Y).^2))); %Estimated Local Noise Level
    % estimation of noise variance
    if mod(ite-1, par.innerIter)==0
        par.nlsp = par.nlsp - par.nlspgap;
        % searching  non-local patches
        blk_arr = Block_Matching( Y, par );
        %         if ite == 1
        %             Sigma = par.nSig0 * ones(size(Sigma));
        %         end
    end
    % Weighted Sparse Coding
    Y_hat = zeros(par.ps2ch, par.maxrc, 'double');
    W_hat = zeros(par.ps2ch, par.maxrc, 'double');
    for i = 1:par.lenrc
        index = blk_arr(:, i);
        nlY = Y( : , index );
        DC = mean(nlY, 2);
        nDCnlY = bsxfun(@minus, nlY, DC);
        Sigma = par.lambdals * sqrt(var(nDCnlY));
        Wls = 1 ./ Sigma;
        % update D
        %         [D, S, ~] = svd( full(nDCnlY), 'econ' );
        [D, S, ~] = svdecon( full(nDCnlY) );
        % update S
        S = sqrt(max( diag(S) .^ 2 - length(index) ./ Wls(1)^2, 0 ));
        % update weight for sparse coding
        Wsc = bsxfun( @rdivide, Sigma.^2, S + eps );
        % update C by soft thresholding
        B = D' * nDCnlY;
        C = sign(B) .* max( abs(B) - Wsc, 0 );
        % update Y
        nDCnlYhat = D * C;
        %         % update Weighting matrix for weighted least square
        %         Wls(index) = par.lambdals ./ sqrt(abs(Sigma(index) - mean((nDCnlY - nDCnlYhat).^2)));
        % add back DC components
        nlYhat = bsxfun(@plus, nDCnlYhat, DC);
        % aggregation
        Y_hat(:, index) = Y_hat(:, index) + bsxfun(@times, nlYhat, Wls);
        W_hat(:, index) = W_hat(:, index) + repmat(Wls, [par.ps2ch, 1]);
    end
    % Reconstruction
    im_out = PGs2Image(Y_hat, W_hat, par);
    % calculate the PSNR and SSIM
    PSNR =   csnr( im_out * 255, par.I * 255, 0, 0 );
    SSIM      =  cal_ssim( im_out * 255, par.I * 255, 0, 0 );
    fprintf('Iter %d : PSNR = %2.4f, SSIM = %2.4f\n', ite, PSNR, SSIM);
    par.PSNR(ite, par.image) = PSNR;
    par.SSIM(ite, par.image) = SSIM;
end
return;

