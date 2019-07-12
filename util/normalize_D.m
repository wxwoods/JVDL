function NormalizeD = normalize_D(D)
    D_size = size(D);
    diagG = 1./sqrt(sum(D.^2)); 
    NormalizeD = D .* repmat(diagG, D_size(1), 1);
    for i = 1:D_size(2)
%         if NormalizeD(1,i) == 'NaN'
%            NormalizeD(:,i) = 0;
%         end
        if isnan(NormalizeD(1,i))
           NormalizeD(:,i) = 0;
        end
    end