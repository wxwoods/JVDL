% @article{Wei2017Dynamical,
%   title={Dynamical Textures Modeling via Joint Video Dictionary Learning},
%   author={Wei, Xian and Li, Yuanxiang and Shen, Hao and Chen, Fang and Kleinsteuber, Martin and Wang, Zhongfeng},
%   journal={IEEE Transactions on Image Processing A Publication of the IEEE Signal Processing Society},
%   volume={PP},
%   number={99},
%   pages={1-1},
%   year={2017},
% }
% Muenchen, 2017. Contact: xian.wei@tum.de

clear all
close all

addpath('util');
% %%%%%%%%%%%%%%%%%%%%%%%
% %GET raw data
%input data
 load ydata;
 size_ydata = size(ydata);
 K = fix(0.5*size_ydata(3)); %number of atoms
%Initialize parameters
 param = init_parameters(size_ydata); 
 ydata = reshape(ydata,size_ydata(1)*size_ydata(2),size_ydata(3));
 %ydataTest = ydata(:,(param.batchSize*param.numberBatch+1):size_ydata(3));
 ydataTest = 0;
 ydataTrain = ydata (:,1:param.batchSize*param.numberBatch);
 
%%preprocessing
 ydataTrain_size = size(ydataTrain);%256*18....
 
 if max(ydataTrain(:))>1
    ydataTrain = ydataTrain/255;
 end
 vecOfMeans = mean(ydataTrain,2);
 ydataTrain   = ydataTrain - repmat(vecOfMeans,1,ydataTrain_size(2)); %according to the row.  

 [B,D] = AVDL(ydataTrain,ydataTest,param.InitialDic(:,1:K),param);

