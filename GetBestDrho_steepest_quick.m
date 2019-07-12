function [D,Drho,estimate_Drho,W,f0,t] = GetBestDrho_steepest_quick(X,B,W,D,D_update,mode,param,it,t,f0)
W = full(W); 
W_SIZE = size(W);

%% main cost function
f_c = f0;
alpha 		= 1e-2; %\in(0,0.5)
val = -D_update(:)'*D_update(:);
%t_init = 0;
dx = -D_update;%approximate optimal point along a way looks like "Z".
if it==1
    W1 = W(:,1:(W_SIZE(2)-1));
    W2 = W(:,2:W_SIZE(2));
    %val = 0;
    temps = ((W2-(B) * W1).^2);
    f0 = sum(temps(:));
    f_c = f0;
    
    %t_init = 0;
    t_init  = sqrt(sum(dx(:).^2));
    t = 1/t_init;
end
%% main cost function
%dx = normalize_D(dx);
iter_D = 0;
max_iter_D = 25;

D0 = D;
beta  		= [0.9,0.5,0.33,0.1,0.09];%is lamda; the most impot=rtant; speed of reduction for step t. should be lager for D

threshold = 1e-4;
    
while (iter_D == 0)||(f0 > f_c + alpha * t * val) &&(iter_D < max_iter_D)
    display1 = f_c + alpha * t * val;
    %display2 = elastic_net_coss_0 + alpha * t * val_ELA 
    temp = abs(f0-display1);
    if temp<threshold
        break;
    end
    if ((temp/abs(display1))>1000)
        %||((elastic_net_coss/display2)>1000)
     t  = t*beta(5); 
%     elseif sum(B(:))==0
%         t=1;
    elseif((temp/abs(display1))>10)
        %||((elastic_net_coss/display2)>10)
     t  = t*beta(4);
    elseif((temp/abs(display1))>0.1)
        %||((elastic_net_coss/display2)>10)
     t  = t*beta(3);
%     elseif iter_D == 1
%         t  = t;
    else
     t  = t*beta(1);
    end
    
     D =  D0+t*dx;
     D  = normalize_D(D);
     
      W = mexLasso(X, D, param.paramLasso);
      W = full(W); 
      
      W1 = W(:,1:(W_SIZE(2)-1));
      W2 = W(:,2:W_SIZE(2));
     temps = ((W2-(B) * W1).^2);
     f0 = sum(temps(:));
     %[elastic_net_coss] = getlossfucfromelasticnet(X,W,D,param);
     
    % f0 = sum(temps(:));
     iter_D = iter_D+1;
     
     sparsity = length(find(W~=0))
end

Drho = t;
t = t/(beta(3)^1);

 estimate_Drho = 0;
 
end
 
 function [elastic_net_coss] = getlossfucfromelasticnet(X,W,D,param)
    elastic_net_coss = 0;
    temp_sum = (X - D*W).^2;  
    elastic_net_coss = 0.5*sum(temp_sum(:));
    Sp_type  = 'PNormAbs';
    q = 1;
    mu       = 0;    % Multiplier in log(1+mu*x^2)
    %[fsp,q_w] = Sparsifying_functions(Sp_type, 'Evaluate', double(W), q, mu);
    fsp = sum(abs(W(:)));
    elastic_net_coss = elastic_net_coss + param.lamda1*fsp + 0.5*param.lamda2*sum(W(:).^2);
 end
 
 function [diff_d_ela] = diff_D_elastic(X,W,D)
    diff_d_ela = D*W*W'-X*W';
 end