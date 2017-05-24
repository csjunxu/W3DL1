clear;
Original_image_dir  =    'C:\Users\csjunxu\Desktop\Projects\WODL\20images\';
% Original_image_dir  =    'C:\Users\csjunxu\Desktop\Projects\WODL\20newimages\';
Sdir = regexp(Original_image_dir, '\', 'split');
fpath = fullfile(Original_image_dir, '*.png');
im_dir  = dir(fpath);
im_num = length(im_dir);


nSig                      =   10;         % The standard variance of the additive Gaussian noise;
sp                        =   0.3;        % salt and pepper
Type = 0;


par.ps = 8; % patch size
par.win = 30;

par.outerIter = 8;
par.innerIter = 2;
for step = [3]
    par.step = step;
    for nlspini = [90]
        par.nlspini = nlspini;
        for nlspgap = [10]
            par.nlspgap = nlspgap;
            for delta = 0
                par.delta = delta;
                for lambda1 = [0.0008]
                    par.lambda1 = lambda1;
                    for lambda2 = [1.1 1.2 1.3]
                        par.lambda2 = lambda2;
                        % record all the results in each iteration
                        par.PSNR = zeros(par.outerIter, im_num, 'double');
                        par.SSIM = zeros(par.outerIter, im_num, 'double');
                        T512 = [];
                        T256 = [];
                        for i = 1:im_num
                            par.nlsp = par.nlspini;  % number of non-local patches
                            par.image = i;
                            par.nSig = nSig;
                            par.I =  double( imread(fullfile(Original_image_dir, im_dir(i).name)) );
                            S = regexp(im_dir(i).name, '\.', 'split');
                            randn('seed',0);
                            par.nim =   par.I + par.nSig*randn(size(par.I));
                            randn('seed',0);
                            [par.nim,Narr]          =   impulsenoise(par.nim,sp, Type);
                            if Type == 0
                                imwrite(par.nim/255, ['GINimages/G' num2str(nSig) '_SPIN' num2str(sp) '_' im_dir(i).name]);
                            elseif Type == 1
                                imwrite(par.nim/255, ['GINimages/G' num2str(nSig) '_RVIN' num2str(sp) '_' im_dir(i).name]);
                            else
                                break;
                            end
                            [par.pim,ind]           =   adpmedft(par.nim,19);
                            
                            fprintf('%s :\n',im_dir(i).name);
                            PSNR =   csnr( par.nim, par.I, 0, 0 );
                            SSIM      =  cal_ssim( par.nim, par.I, 0, 0 );
                            fprintf('The initial value of PSNR = %2.4f, SSIM = %2.4f \n', PSNR,SSIM);
                            
                            PSNR =   csnr( par.pim, par.I, 0, 0 );
                            SSIM      =  cal_ssim( par.pim, par.I, 0, 0 );
                            fprintf('The initial value of PSNR = %2.4f, SSIM = %2.4f \n', PSNR,SSIM);
                            %
                            time0 = clock;
                            [im_out, par]  =  WLSWSC_Sigma_WAGIN(par);
                            if size(par.I,1) == 512
                                T512 = [T512 etime(clock,time0)];
                                fprintf('Total elapsed time = %f s\n', (etime(clock,time0)) );
                            elseif size(par.I,1) ==256
                                T256 = [T256 etime(clock,time0)];
                                fprintf('Total elapsed time = %f s\n', (etime(clock,time0)) );
                            end
                            im_out(im_out>255)=255;
                            im_out(im_out<0)=0;
                            % calculate the PSNR
                            par.PSNR(par.outerIter, par.image)  =   csnr( im_out, par.I, 0, 0 );
                            par.SSIM(par.outerIter, par.image)      =  cal_ssim( im_out, par.I, 0, 0 );
                            %             imname = sprintf('nSig%d_clsnum%d_delta%2.2f_lambda%2.2f_%s', nSig, cls_num, delta, lambda, im_dir(i).name);
                            %             imwrite(im_out/255,imname);
                            fprintf('%s : PSNR = %2.4f, SSIM = %2.4f \n',im_dir(i).name, par.PSNR(par.outerIter, par.image),par.SSIM(par.outerIter, par.image)     );
                        end
                        mPSNR=mean(par.PSNR,2);
                        [~, idx] = max(mPSNR);
                        PSNR =par.PSNR(idx,:);
                        SSIM = par.SSIM(idx,:);
                        mSSIM=mean(SSIM,2);
                        mT512 = mean(T512);
                        sT512 = std(T512);
                        mT256 = mean(T256);
                        sT256 = std(T256);
                        fprintf('The best PSNR result is at %d iteration. \n',idx);
                        fprintf('The average PSNR = %2.4f, SSIM = %2.4f. \n', mPSNR(idx),mSSIM);
                        name = sprintf(['WLSWSC_Sigma_WAGIN_' Sdir{end-1} '_nSig' num2str(nSig) '_step' num2str(par.step) '_nlspini' num2str(par.nlspini) '_nlspgap' num2str(par.nlspgap) '_delta' num2str(delta) '_l1' num2str(lambda1) '_l2' num2str(lambda2) '.mat']);
                        save(name,'nSig','PSNR','SSIM','mPSNR','mSSIM','mT512','sT512','mT256','sT256');
                    end
                end
            end
        end
    end
end