function  [im_out, Par] = WLSWSC_Sigma_W1W2GIN(Par)
im_out    =   Par.pim;
% parameters for noisy image
[h,  w, ch]      =  size(im_out);
Par.h = h;
Par.w = w;
Par.ch = ch;
Par = SearchNeighborIndex( Par );
% original noisy image to patches
NY = Image2PatchNew( Par.nim, Par );
PY = Image2PatchNew( Par.pim, Par );
for ite  =  1 : Par.Outerloop
    % iterative regularization
    im_out = im_out + Par.delta * (Par.pim - im_out);
    % image to patches
    Y = Image2PatchNew( im_out, Par );
    % estimate local noise variance
    %     for c = 1:Par.ch
    %         if (ite == 1) && (Par.Outerloop > 1)
    % %             TempSigmaRow = repmat(Par.nSig(c)^2, Par.ps2, 1);
    % %             TempSigmaRow = exp( - Par.lambda1*TempSigmaRow);
    %             SigmaRow((c-1)*Par.ps2+1:c*Par.ps2, :) = TempSigmaRow;
    %             %  TempSigmaRow = sqrt(abs(repmat(Par.nSig(c)^2, 1, size(Y, 2))));
    %         else
    %             SigmaRow((c-1)*Par.ps2+1:c*Par.ps2, :) = exp( - Par.lambda1*max(0, repmat(Par.nSig(c)^2, Par.ps2, length(Par.SelfIndex)) - ErrorRow((c-1)*Par.ps2+1:c*Par.ps2, :)));
    %             %  TempSigmaRow = Par.lambda1*sqrt(max(0, repmat(Par.nSig(c)^2, Par.ps2, length(Par.SelfIndex)) - ErrorRow((c-1)*Par.ps2+1:c*Par.ps2, :)));
    %             %  TempSigmaRow = Par.lambda1*sqrt(abs(repmat(Par.nSig(c)^2, 1, size(Y, 2)) - mean((NY((c-1)*Par.ps2+1:c*Par.ps2, :) - Y((c-1)*Par.ps2+1:c*Par.ps2, :)).^2)));
    %         end
    %     end
    % estimate local noise variance, par.lambdals is put here since the MAP and Bayesian rules
    SigmaRow = (NY-Y).^2;
    SigmaCol = Par.lambda2 * sqrt(abs(repmat(Par.nSig^2, 1, size(Y,2)) - mean((PY - Y).^2))); % Estimated Local Noise Level
    % estimation of noise variance
    if mod(ite-1, Par.Innerloop)==0
        Par.nlsp = Par.nlsp - Par.nlspgap;
        % searching  non-local patches
        [blk_arr, wei_arr] = Block_Matching( Y, Par );
        %         if ite == 1 && Par.lambda1 ~= 0
        %             SigmaRow = ones(Par.ps2ch, length(Par.SelfIndex));
        %         end
        if ite == 1 && Par.lambda2 ~= 0
            SigmaCol = Par.nSig * ones(size(SigmaCol));
        end
    end
    % Weighted Sparse Coding
    Y_hat = zeros(Par.ps2ch, Par.maxrc, 'double');
    W_hat = zeros(Par.ps2ch, Par.maxrc, 'double');
    for i = 1:Par.lenrc
        index = blk_arr(:, i);
        wei = wei_arr(:, i);
        nlY = Y( : , index );
        DC = bsxfun(@times, nlY, wei');
        %         DC = mean(nlY, 2);
        nDCnlY = bsxfun(@minus, nlY, DC);
        
        % Compute W1 and W2
        W1 = exp( - Par.lambda1*mean(SigmaRow(:, index), 2)); % max?
        W2 = 1 ./ (SigmaCol(index) + eps);
        
        % update D
        [D, S, ~] = svd( full(nDCnlY), 'econ' );
        % update S
        S = sqrt(max( diag(S).^2 - length(index) * SigmaCol(index(1))^2, 0 ));
        if Par.lambda1 == 0
            % update weight for sparse coding
            Wsc = bsxfun( @rdivide, SigmaCol(index).^2, S + eps );
            % update C by soft thresholding
            B = D' * nDCnlY;
            C = sign(B) .* max( abs(B) - Wsc, 0 );
            % update Y
            nDCnlYhat = D * C;
        elseif Par.model == 1 && Par.lambda1 ~= 0
            S = diag(S);
            % min |Z|_1 + |W1(Y-DSC)W2|_F,2  s.t.  C=Z
            C = W3SC_ADMM( nDCnlY, D, S, W1, W2, Par );
            % update Y
            nDCnlYhat = D * S * C;
        elseif Par.model == 2 && Par.lambda1 ~= 0
            S = diag(S);
            % min |C|_2 + |W1(Y-DSC)W2|_F,2
            % Solve Sylvester equation AX + XB = E for X
            A = S' * D' * diag(W1.^2) * D * S;
            B = diag(1./(W2.^2));
            E = S' * D' * diag(W1.^2) * nDCnlY;
            C = sylvester(A, B, E);
            % update Y
            nDCnlYhat = D * S * C;
        end
        
        % add back DC components
        nlYhat = bsxfun(@plus, nDCnlYhat, DC);
        
        %         % aggregation
        %         Y_hat(:, index) = Y_hat(:, index) + nlYhat;
        %         W_hat(:, index) = W_hat(:, index) + ones(Par.ps2ch, Par.nlsp);
        %         % aggregation
        %         Y_hat(:, index) = Y_hat(:, index) + bsxfun(@times, nlYhat, W2);
        %         W_hat(:, index) = W_hat(:, index) + repmat(W2, [par.ps2ch, 1]);
        % aggregation
        Y_hat(:, index) = Y_hat(:, index) + bsxfun(@times, bsxfun(@times, W1, nlYhat), W2);
        W_hat(:, index) = W_hat(:, index) + W1 * W2;
    end
    % Reconstruction
    im_out = PGs2Image(Y_hat, W_hat, Par);
    % calculate the PSNR and SSIM
    PSNR =   csnr( im_out , Par.I , 0, 0 );
    SSIM      =  cal_ssim( im_out , Par.I , 0, 0 );
    fprintf('Iter %d : PSNR = %2.4f, SSIM = %2.4f\n', ite, PSNR, SSIM);
    Par.PSNR(ite, Par.image) = PSNR;
    Par.SSIM(ite, Par.image) = SSIM;
end
return;

