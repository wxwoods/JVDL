function D_update = D_diff(W,D,B,X,param)
Dsize = size(D);
D_update = zeros(Dsize);
%Phi = param.SamplingMatrix;
Wsize = size(W);
tic
for i = 1:(Wsize(2)-1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   derivative of W_k
%% first part
    %update_D_no = i
    delta_index1 = find(W(:,i)~=0);%find the index where element is not equal to zero
    delta_index1_size = size(delta_index1,1);% find nonzero number
    Pt = InversK(D(:,delta_index1),param.lamda2);
        %P_tPt = InvK*InvK';
    Firstpart1 = -X(:,i)*( W(:,i+1) - B(:,delta_index1)*W(delta_index1,i) )'*B(:,delta_index1)*Pt;
    for j = 1:delta_index1_size
        if W(delta_index1(j),i)>0
            S(j) = 1;
        else
            S(j) = -1;
        end
    end
     same_bracket = (D(:,delta_index1)'*X(:,i)- param.lamda1*S');%(D'X-ZS)
     Qt = same_bracket*( W(:,i+1) - B(:,delta_index1)*W(delta_index1,i) )'*B(:,delta_index1);
    Firstpart =  Firstpart1+D(:,delta_index1)*Pt*(Qt+Qt')*Pt;
    
    D_update(:,delta_index1) =  D_update(:,delta_index1) + Firstpart;
    clear Pt Qt S delta_index1 delta_index1_size same_bracket Firstpart Firstpart1;
%% second part
%%K+1  
    delta_index2 = find(W(:,i+1)~=0);%find the index where element is not equal to zero
    delta_index2_size = size(delta_index2,1);
        Pt1 = InversK(D(:,delta_index2),param.lamda2);
    secondpart1 = X(:,i+1)*( W(delta_index2,i+1) - B(delta_index2,:)*W(:,i) )'*Pt1;
    for j = 1:delta_index2_size
        if W(delta_index2(j),i+1)>0
            S(j) = 1;
        else
            S(j) = -1;
        end
    end
    same_bracket2 = (D(:,delta_index2)'*X(:,i+1)- param.lamda1*S');
    Qt1 = same_bracket2*( - W(delta_index2,i+1) + B(delta_index2,:)*W(:,i) )';
    
    Secondpart =  secondpart1 + D(:,delta_index2)*Pt1*( Qt1+Qt1' )*Pt1;
    
    D_update(:,delta_index2) =  D_update(:,delta_index2) + Secondpart;
    clear Pt1 Qt1 S delta_index2 delta_index2_size same_bracket2 Secondpart secondpart1;
end
    
    if param.Inco == 1
        B_norm = normalize_D(B);
        D_update = D_update + param.lamda4*(D*B_norm*B_norm');
    end
toc

function InversValue = InversK(D_delta,lamda2)
Dsize = size(D_delta);
identity_matrix = eye(Dsize(2));
temp = D_delta'*D_delta+lamda2*identity_matrix;
InversValue = inv(temp);