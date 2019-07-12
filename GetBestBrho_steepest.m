function [B,Brho,estimate_Brho,t,f0] = GetBestBrho_steepest(W,B,B_update,mode,it,t,f0)
W = full(W); 
W_SIZE = size(W);
W1 = W(:,1:(W_SIZE(2)-1));
W2 = W(:,2:W_SIZE(2));
val = 0;
f_c = f0;
if it == 1
    temps = ((W2-(B) * W1).^2);
    f0 = sum(temps(:));
    f_c = f0;
    
    t = 1.1;
end
alpha 		= 1e-2; %\in(0,0.5)
val = -B_update(:)'*B_update(:);
dx = -B_update;%approximate optimal point along a way looks like "Z".
iter_B = 0;
max_iter_B = 80;
%value_val = trace(val);
B0 = B;

beta  		= [0.9,0.7,0.33,0.1,0.01];%is lamda; the most impot=rtant; speed of reduction for step t. should be lager for D
while (iter_B == 0)||(f0 > f_c + alpha * t * val) && (iter_B<max_iter_B)
    display0 = f_c + alpha * t * val;
    if ((f0/display0)>10000)
     t  = t*beta(5); 
%     elseif sum(B(:))==0
%         t=1;
    elseif((f0/display0)>10)
     t  = t*beta(3);
    else
     t  = t*beta(1);
    end
     %t  = t*beta; 
     B =  B0+t*dx;
     temps = ((W2-(B) * W1).^2);
     f0 = sum(temps(:));
     iter_B = iter_B+1;
end

ii = 1;
Nr = 10;
ff0 = zeros(Nr,1);
tt = zeros(Nr,1);
tt(1) = t;
while ii<=Nr
     tt(ii)  = t*beta(1)^ii;
     %t  = t*beta; 
     BB =  B0+tt(ii)*dx;
     temps = ((W2-(BB) * W1).^2);
     ff0(ii) = sum(temps(:));
     %iter_B = iter_B+1;
     ii = ii+1;
end
[min_ff0,index] = min(ff0);
if min_ff0<f0
    t = tt(index);
    B =  B0+t*dx;
    f0 = min_ff0;
end


Brho = t;
t = t/(beta(3)^2);
estimate_Brho = 0;
%B = B - Brho*B_update;

