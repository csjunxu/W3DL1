function  [im_out, Par] = WLSWSC_Sigma_W1W2GINcw(Par)
im_out    =   Par.pim;
% parameters for noisy image
[h,  w, ch]      =  size(im_out);
Par.h = h;
Par.w = w;
Par.ch = ch;
Par = SearchNeighborIndex( Par );
% original noisy image to patches
[NX, NY] = Image2PG( Par.nim, Par );
[PX, PY]  = Image2PG( Par.pim, Par );
for ite  =  1 : Par.Outerloop
    % iterative regularization
    im_out = im_out + Par.delta * (Par.pim - im_out);
    % image to patches
    [X, Y] = Image2PG( im_out, Par );
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
        %         [blk_arr, wei_arr] = Block_Matching( X, Par );
        blk_arr = Block_Matching0( X, Par );
        %         if ite == 1 && Par.lambda1 ~= 0
        %             SigmaRow = ones(Par.ps2ch, length(Par.SelfIndex));
        %         end
        if ite == 1 && Par.lambda2 ~= 0
            SigmaCol = Par.nSig * ones(size(SigmaCol));
        end
    end
    % Weighted Sparse Coding
    Y_hat = zeros(Par.ps2ch, Par.lenrc, 'double');
    W_hat = zeros(Par.ps2ch, Par.lenrc, 'double');
    %     Y_hat = zeros(Par.ps2ch, Par.maxrc, 'double');
    %     W_hat = zeros(Par.ps2ch, Par.maxrc, 'double');
    for i = 1:Par.lenrc
        index = blk_arr(:, i);
        nlY    = X( : , index );
        y        = Y(:, i);
        DC = mean(nlY, 2);
        nDCnlY = bsxfun(@minus, nlY, DC);
        nDCy = y - DC;
        % obtain D and S via SVD
        [D, S, ~] = svd( full(nDCnlY), 'econ' );
        S = sqrt(max( diag(S).^2 - length(index) * SigmaCol(1)^2, 0 ));
        if Par.lambda1 == 0
            % update weight for sparse coding
            Wsc = bsxfun( @rdivide, SigmaCol(index).^2, S + eps );
            % update C by soft thresholding
            B = D' * nDCy;
            C = sign(B) .* max( abs(B) - Wsc, 0 );
            % update Y
            nDCyhat = D * C;
        elseif Par.model == 1 && Par.lambda1 ~= 0
            S = diag(S);
            % Compute W1 and W2
            W1 = exp( - Par.lambda1 .* SigmaRow(:, i) );
            W2 = 1 ./ (SigmaCol(i) + eps);
            C    = W3SCcw_ADMM( nDCy, D, S, W1, W2, Par );
            nDCyhat = D * S * C;
        elseif Par.model == 2 && Par.lambda1 ~= 0
            S = diag(S);
            % Compute W1 and W2
            W1 = exp( - Par.lambda1 * SigmaRow(:, i) );
            W2 = 1 / (SigmaCol(i) + eps);
            % Solve min |C|_2 + |W1(Y-DSC)W2|_F,2
            Temp1 = S' * D' * diag(W1.^2);
            A = Temp1 * D * S;
            B = (1/W2)^2 * eye(size(S, 2));
            E = Temp1 * nDCy;
            C = (A + B) \ E;
            % update Y
            nDCyhat = D * S * C;
        end
        
        % add back DC components
        nlYhat = bsxfun(@plus, nDCyhat, DC);
        
        % aggregation
        Y_hat(:, i) = Y_hat(:, i) + nlYhat;
        W_hat(:, i) = W_hat(:, i) + ones(Par.ps2ch, 1);
        %         % aggregation
        %         Y_hat(:, index) = Y_hat(:, index) + nlYhat;
        %         W_hat(:, index) = W_hat(:, index) + ones(Par.ps2ch, Par.nlsp);
        %         % aggregation
        %         Y_hat(:, index) = Y_hat(:, index) + bsxfun(@times, nlYhat, W2);
        %         W_hat(:, index) = W_hat(:, index) + repmat(W2, [par.ps2ch, 1]);
        %         % aggregation
        %         Y_hat(:, index) = Y_hat(:, index) + bsxfun(@times, bsxfun(@times, W1, nlYhat), W2);
        %         W_hat(:, index) = W_hat(:, index) + W1 * W2;
    end
    % Reconstruction
    im_out = Patch2Image(Y_hat, W_hat, Par);
    %     im_out = PGs2Image(Y_hat, W_hat, Par);
    
    % calculate the PSNR and SSIM
    PSNR =   csnr( im_out , Par.I , 0, 0 );
    SSIM      =  cal_ssim( im_out , Par.I , 0, 0 );
    fprintf('Iter %d : PSNR = %2.4f, SSIM = %2.4f\n', ite, PSNR, SSIM);
    Par.PSNR(ite, Par.image) = PSNR;
    Par.SSIM(ite, Par.image) = SSIM;
end
return;

