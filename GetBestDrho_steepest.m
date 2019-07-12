function [D,Drho,estimate_Drho] = GetBestDrho_steepest(X,B,W,D,D_update,mode,param)
W = full(W); 
W_SIZE = size(W);
W1 = W(:,1:(W_SIZE(2)-1));
W2 = W(:,2:W_SIZE(2));
val = 0;

%% main cost function
temps = ((W2-(B) * W1).^2);
f0 = sum(temps(:));
f_c = f0;

%% for elastic net
% min_{alpha} 0.5||x-Dalpha||_2^2 + lambda||alpha||_1 +0.5 lambda2||alpha||_2^2
[elastic_net_coss_0] = getlossfucfromelasticnet(X,W,D,param);
elastic_net_coss = elastic_net_coss_0;


alpha 		= 1e-2; %\in(0,0.5)

%% main cost function
val = -D_update(:)'*D_update(:);
t_init = 0;
dx = -D_update;%approximate optimal point along a way looks like "Z".
%dx = normalize_D(dx);

t_init  = sqrt(sum(dx(:).^2));
t = 1/t_init;

iter_D = 0;
max_iter_D = 10;
%value_val = trace(val);

[diff_d_ela] = diff_D_elastic(X,W,D);
dx_ELA = -diff_d_ela;
%dx_ELA = normalize_D(dx_ELA);
val_ELA = diff_d_ela(:)'*dx_ELA(:);
t2 = 1/sqrt(sum(dx_ELA(:).^2));

t = min(t,t2);

D0 = D;
beta  		= [0.7,0.5,0.3,0.1,0.01];%is lamda; the most impot=rtant; speed of reduction for step t. should be lager for D
while (mode == 2)&&(f0 > f_c + alpha * t * val) || (elastic_net_coss > elastic_net_coss_0 + alpha * t * val_ELA)|| (iter_D == 0)&&(iter_D<max_iter_D)
    display1 = f_c + alpha * t * val
    display2 = elastic_net_coss_0 + alpha * t * val_ELA 
    if ((f0/display1)>100)||((elastic_net_coss/display2)>100)
     t  = t*beta(4); 
    elseif((f0/display1)>10)||((elastic_net_coss/display2)>10)
     t  = t*beta(4);
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
     
     
     [elastic_net_coss] = getlossfucfromelasticnet(X,W,D,param);
     
    % f0 = sum(temps(:));
     iter_D = iter_D+1;
end

Drho = t;
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
    fsp = sum(abs(X(:)));
    elastic_net_coss = elastic_net_coss + param.lamda1*fsp + 0.5*param.lamda2*sum(W(:).^2);
 end
 
 function [diff_d_ela] = diff_D_elastic(X,W,D)
    diff_d_ela = D*W*W'-X*W';
 end
