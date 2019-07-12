% Xian Wei, Research Group for Geometric Optimization and Machine Learning
% Muenchen, 2014. Contact: xian.wei@tum.de

function param = init_parameters(size_ydata)

%initialize dictionary
    param.channel_initial_dic = 'DCT';
    param.channel_initial_dic = 'Random';
    param.channel_initial_dic = 'From_data';

    K = fix(0.5*size_ydata(3)); %number of atoms
    bb =sqrt(size_ydata(1)*size_ydata(2)); %dim of dictionary
    k_lamda =[-3,-2,-1,0,1,2,3]; 
    % initial dictionary
    %1 DCT dic
    if strcmpi(param.channel_initial_dic, 'DCT')
        Pn=ceil(sqrt(K));
        DCT=zeros(bb,Pn);
        for k=0:1:Pn-1,
            V=cos([0:1:bb-1]'*k*pi/Pn);
            if k>0, V=V-mean(V); end;
            DCT(:,k+1)=V/norm(V);
        end;
        DCT=kron(DCT,DCT);
    
        param.InitialDic = DCT(:,1:K );
        clear DCT V Pn;
        % param.batchSize = 1;   % length of each dataset
        param.lamda1       = 0.5+0.025*k_lamda(4);    
        %param.lamda1       = 0.15;
        %param.lamda1       =  0.001; 
        param.lamda2      = 0;  
        
    elseif strcmpi(param.channel_initial_dic, 'Random')
        M = K; % one of Second of original data 
        N = size_ydata(1)*size_ydata(2); %1024
        Phi = randn(N,M);
        Phi = Phi./repmat(sqrt(sum(Phi.^2,1)),[N,1]);	
        param.InitialDic = Phi;
    
        param.lamda1       = 0.25+0.025*k_lamda(4);               
        param.lamda2       = 0;                 
    
    %3 Get dic from data
    else   
        load Initial_dic.mat;
        param.lamda1       = 0.07;               
        param.lamda2      = 0.00001;  
        param.lamda2      = 0;
        param.InitialDic = D0;
    end

%initialize other parameters
    N = size(param.InitialDic,1);        
    param.mainIt        = 102; % Maximal main iterations.
    %param.mainIt        = 20; % Maximal main iterations.

     param.method = 'AVDL';
     param.dataPeak = 1;
     param.dataLength = size_ydata(3);   % length of each dataset 
     param.batchSize = param.dataLength;   % length of each dataset 
     param.targetMSE = 1e-6;
     if ~isfield(param,'numberBatch');  param.numberBatch   = fix ( param.dataLength/param.batchSize );      end;
    
    param.paramDicUpdating = struct(...
            'eps',    param.targetMSE, ...
            'L',      floor(0.5*N)  ); 
    param.paramLasso = struct(...
            'mode',   2, ...            
            'lambda', param.lamda1, ...
            'lambda2', param.lamda2, ...
            'L',      floor(0.9*N*10)  );
    param.lamda3 = 0.1;
    param.Inco = 1; % penalty for incoherence between B and D
    param.lamda4 = 1e-5;
end