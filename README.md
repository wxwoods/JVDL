# JVDL
Dynamical Textures Modeling via Joint Video Dictionary Learning  


Video representation is an important and challenging task in the computer vision community. In this paper, we
consider the problem of modeling and classifying video sequences
of dynamic scenes which could be modeled in a dynamic textures
(DT) framework. At first, we assume that image frames of a
moving scene can be modeled as a Markov random process.
We propose a sparse coding framework, named joint video
dictionary learning (JVDL), to model a video adaptively. By
treating the sparse coefficients of image frames over a learned
dictionary as the underlying “states”, we learn an efficient and
robust linear transition matrix between two adjacent frames of
sparse events in time series. Hence, a dynamic scene sequence is
represented by an appropriate transition matrix associated with
a dictionary. In order to ensure the stability of JVDL, we impose
several constraints on such transition matrix and dictionary. The
developed framework is able to capture the dynamics of a moving
scene by exploring both sparse properties and the temporal correlations of consecutive video frames. Moreover, such learned JVDL
parameters can be used for various DT applications, such as
DT synthesis and recognition. Experimental results demonstrate
the strong competitiveness of the proposed JVDL approach in
comparison with state-of-the-art video representation methods.
Especially, it performs significantly better in dealing with DT
synthesis and recognition on heavily corrupted data.  

cite：   
@article{Wei2017Dynamical,  
  title={Dynamical Textures Modeling via Joint Video Dictionary Learning},  
  author={Wei, Xian and Li, Yuanxiang and Shen, Hao and Chen, Fang and Kleinsteuber, Martin and Wang, Zhongfeng},  
  journal={IEEE Transactions on Image Processing A Publication of the IEEE Signal Processing Society},  
  volume={PP},  
  number={99},  
  pages={1-1},  
  year={2017},  
}
