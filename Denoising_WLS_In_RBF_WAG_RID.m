function  [im_out, par] = Denoising_WLS_In_RBF_WAG_RID(par)
im_in = par.nim;
im_out    =   par.nim;
if ~isfield(par, 'nSig')
    par.nSig0 = 0;
end
% parameters for noisy image
[h,  w, ch]      =  size(im_out);
par.h = h;
par.w = w;
par.ch = ch;
par = SearchNeighborIndex( par );
Wls = ones( 1, par.maxrc );
for ite  =  1 : par.outerIter
    %     % iterative regularization
    %     im_out = im_out+par.delta*(par.nim - im_out);
    % image to patches and estimate local noise variance
    [Y, Sigma] = Image2Patch( im_out, im_in, par );
    % estimation of noise variance
    if mod(ite-1, par.innerIter)==0
        par.nlsp = par.nlsp - 10;
        % searching  non-local patches
        blk_arr = Block_Matching( Y, par );
%         if ite == 1
%             Sigma = par.nSig0 * ones(size(Sigma));
%         end
    end
    % Weighted Sparse Coding
    Y_hat = zeros(par.ps2ch, par.maxrc, 'single');
    W_hat = zeros(par.ps2ch, par.maxrc, 'single');
    for i = 1:par.lenrc
        index = blk_arr(:, i);
        nlY = Y( : , index );
        DC = mean(nlY, 2);
        nDCnlY = bsxfun(@minus, nlY, DC);
        % update D and S
        [D, S, ~] = svd( full(nDCnlY), 'econ' );
        %         [D, S, ~] = svd(cov(nDCnlY'), 'econ');
        S = diag(S);
        % update weight for sparse coding
        Wsc = bsxfun( @rdivide, par.lambda * Wls(index).^2, sqrt(S) + eps );
        % update C by soft thresholding
        B = D' * nDCnlY;
        C = sign(B) .* max( abs(B) - Wsc, 0 );
        % update Y
        nDCnlYhat = D * C;
        % update weight for least square
        Wls(index) = exp( - par.lambdaW * sqrt(sum((nDCnlY - nDCnlYhat) .^2, 1)) );
        %         % Recovered Estimated Patches
        %         nDCnlYhat = WODL(nDCnlY, Wls, C, par);
        % add back DC components
        nlYhat = bsxfun(@plus, nDCnlYhat, DC);
        % aggregation
        Y_hat(:, index) = Y_hat(:, index) + bsxfun(@times, nlYhat, Wls(index));
        W_hat(:, index) = W_hat(:, index) + repmat(Wls(index), [par.ps2ch, 1]);
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

