% Xian Wei, Research Group for Geometric Optimization and Machine Learning
% Muenchen, 2014. Contact: xian.wei@tum.de

function [B,D] = AVDL(Xtrain,Xtest,D0,param)
%preprocessing data
Dic_size = size(D0);
N = Dic_size(1);%512
K = Dic_size(2); % number of atoms 1024
param.targetMSE = param.targetMSE*param.batchSize;
param.sparsity = floor (N/10);
   
% initialize result struct
%% initial Dictionary, D0, D and res.D, and C 
    D = D0;
    it = 0;    % initial number of main (outer) iterations
    %%the following parameters are designed for memorizing the process of
    %%methods, they can be commanted directly
    rho_B =  zeros(param.mainIt+2,param.numberBatch);
    rho_D =  zeros(param.mainIt+2,param.numberBatch);
    estimate_Drho =  zeros(param.mainIt+2,param.numberBatch);
    estimate_Brho =  zeros(param.mainIt+2,param.numberBatch);    
    Error_main =  ones(param.mainIt+2,1);                      % error for \| Xtrain - DW\|_F^2
    Error_trans_sparse =  zeros(param.mainIt+2,1);              % error for \sum_i{\| W_{i+1} - BW_i\|_F^2}, B is learned by gradient method
    Error_trans_sparse_inverse =  zeros(param.mainIt+2,1);      % error for \sum_i{\| W_{i+1} - BW_i\|_F^2}, B is directly gotten from W.
    Norm_B =  zeros(param.mainIt+2,1);
    Eig_B  =  zeros(param.mainIt+2,1);    
    Norm_B_inverse =  zeros(param.mainIt+2,1);   %check the eigenvalues of B, judge how far it is away from 1.
    Eig_B_inverse  =  zeros(param.mainIt+2,1);   %check the eigenvalues of B, judge how far it is away from 1. 
    sparsity =  zeros(param.mainIt+2,1);         %check the sparsity of W.
   
    t_B = 0;  f0_B = 0; f0_D = 0;  t_D = 0;
    %% the main loop
    while 1
        it = it + 1
        %
        if sum(Error_main(it,:))< 1e-2;
            break; 
        end
        if (it > param.mainIt) 
            break; 
        end
        fprintf('---  iteration %4i, learning with data set   ',it);
        temp = 0;

        for i = 1:param.numberBatch %Designed for stochastic gradient descent method, or param.numberBatch = 1 when steepest gradient mathod
            temp = temp+1
            it  
            X =  Xtrain;
            DicSampMatrix = D;      
            if it ==1
                W = mexLasso(X, DicSampMatrix, param.paramLasso);
                % W = mexOMP(X, DicSampMatrix, param.paramOMP);
                %update B 
                W = full(W);         
                temp_sum = (X - D*W).^2;  
                Error_main(1) = sum(temp_sum(:)); %norm2 of error, L2 fit for residence : |X-DW|_2
                sparsity(1) = length(find(W~=0));
            end
             Current_sparsity = length(find(W~=0))
             size_X = size(X);
             Current_sparsity_rate = Current_sparsity/(size_X(1)*size_X(2))
            if (i == 1)&&(it==1)        
                B = zeros(size(W,1),size(W,1));
                B_update = B_diff(W,B,D,param);
                mode_B =2 ;%line search
                [B,rho_B(it,i),estimate_Brho(it,i),t_B,f0_B]=GetBestBrho_steepest(W,B,B_update,mode_B,it,t_B,f0_B);          
                W_SIZE = size(W);
                temp_sum = (W(:,2:W_SIZE(2))-B*W(:,1:(W_SIZE(2)-1))).^2;
                Error_trans_sparse(1) = sum(temp_sum(:)); %norm2 of error, L2 fit for residence : |X-DW|
            
                B_inverse  = W(:,2:param.batchSize)*W(:,1:(param.batchSize-1))'*pinv(W(:,1:(param.batchSize-1))*W(:,1:(param.batchSize-1))');%1024*1024
                temp_sum = (W(:,2:W_SIZE(2))-B_inverse*W(:,1:(W_SIZE(2)-1))).^2;
                Error_trans_sparse_inverse(1) = sum(temp_sum(:));
            
                Norm_B(1) = norm(B,'fro');
                temp         = eig(B);
                Eig_B(1)  = max(temp);
      
                B_inverse  = W(:,2:param.batchSize)*W(:,1:(param.batchSize-1))'*pinv(W(:,1:(param.batchSize-1))*W(:,1:(param.batchSize-1))');%1024*1024
                Norm_B_inverse(1) = norm(B_inverse,'fro');
                temp         = eig(B_inverse);
                Eig_B_inverse(1)  = max(temp);
            
            else
                B_update = B_diff(W,B,D,param);
                mode_B =2 ;%line search
                [B,rho_B(it,i),estimate_Brho(it,i),t_B,f0_B]=GetBestBrho_steepest(W,B,B_update,mode_B,it,t_B,f0_B);      
            end
            %update Dictionary
            %D_update= D_diff_no_phi_Combin_Matrix2(W,D,B,X,param);
            D_update= D_diff(W,D,B,X,param);      
            D_update = D_update - D*diag(diag(D'*D_update));
            mode_D =2;
            % [D,rho_D(it,i),estimate_Drho(it,i)]=GetBestDrho_steepest(X,B,W,D,D_update,mode,param);
            [D,rho_D(it,i),estimate_Drho(it,i),W,f0_D,t_D]=GetBestDrho_steepest_quick(X,B,W,D,D_update,mode_D,param,it,t_D,f0_D);        
            %[D,rho_D(it,i),estimate_Drho(it,i),W]=test_BD(X,B,W,D,D_update,mode_D,param);
            
            Wtrain = W;
        end   % useSett   
 
        temp_sum = (Xtrain - D*Wtrain).^2;  
        Error_main(it+1) = sum(temp_sum(:)); %norm2 of error, L2 fit for residence : |X-DW|_2
        %Error_main(it) = norm(temp_sum); %norm2 of error, L2 fit for residence : |X-DW|_2
      
        W_SIZE = size(Wtrain);
        temp_sum = (Wtrain(:,2:W_SIZE(2))-B*Wtrain(:,1:(W_SIZE(2)-1))).^2;
        Error_trans_sparse(it+1) = sum(temp_sum(:)); % norm2 of error, L2 fit for residence : |X-DW|_2
        %Error_main_sparse(it) = norm(temp_sum);    % norm2 of error, L2 fit for residence : |X-DW|_2
      
        B_inverse  = W(:,2:param.batchSize)*W(:,1:(param.batchSize-1))'*pinv(W(:,1:(param.batchSize-1))*W(:,1:(param.batchSize-1))');%1024*1024
        temp_sum = (W(:,2:W_SIZE(2))-B_inverse*W(:,1:(W_SIZE(2)-1))).^2;
        Error_trans_sparse_inverse(it+1) = sum(temp_sum(:));
      
        sparsity(it+1) = length(find(Wtrain~=0));
      
        Norm_B(it+1) = norm(B,'fro');
        temp         = eig(B);
        Eig_B(it+1)  = max(temp);
    
        Norm_B_inverse(it+1) = norm(B_inverse,'fro');
        temp         = eig(B_inverse);
        Eig_B_inverse(it+1)  = max(temp);
    end
    B_update = B_diff(W,B,D,param);
    [B,rho_B(it,i),estimate_Brho(it,i),t_B,f0_B]=GetBestBrho_steepest(W,B,B_update,mode_B,it,t_B,f0_B);
return


    
    



