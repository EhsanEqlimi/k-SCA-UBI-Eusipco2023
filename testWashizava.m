clc;clear;close all
 A=[0.5525 0.3919 0.5707 0.3934 0.6904;0.6863 -0.6066 0.5166 0.8634 -0.6007;0.4730 0.6917 -0.6383 -0.3158 0.4032];%Whashizava example
%%
A_mov=[0.5536 0.6930 -0.3908 0.3892 -0.5712;
0.6829 -0.5980 -0.8644 -0.6083 -0.5154;
0.4765 0.4027 0.3165 0.6917 0.6388];
[MixingIdentificatioError, MixingVectorerror,NMSE,NMSSum,AhatNew,Norm2Error,Error] =FnMixingIdentificationError(A,A_mov)

A_wash=[0.3865 0.3935 0.55251 0.3981 0.6907;
-0.6082 0.86341 0.6861 0.5431 -0.6004;
0.6933 -0.31581 0.4733 -0.7393 0.4031];


[MixingIdentificatioError_wash, MixingVectorerror_wash,NMSE_wash,NMSSum_wash,AhatNew_wash,Norm2Error_wash,Error_wash] =FnMixingIdentificationError(A,A_wash)