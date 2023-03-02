function [MixingIdentificatioError, MixingVectorerror,NMSE,NMSSum,AhatNew,Norm2Error,Error,BAS_Detected] = FnMixingIdentificationError(A,Ahat)
%% Error estimation
[AhatNew,M]=nearest2(Ahat,(A));
BAS=0;
n=min([size(A,2) ,size(AhatNew,2)]);
for e=1:n
    MixingVectorerror(e)=acos((AhatNew(:,e)'* A(:,e))/(norm(AhatNew(:,e))*norm(A(:,e))));
    BAS=MixingVectorerror(e)+BAS;
end
MixingIdentificatioError=abs(BAS);
MixingVectorerror=abs(MixingVectorerror);
BAS_Detected=sum(MixingVectorerror(find(MixingVectorerror<.1)));
% AhatNew
AhatNew;

%% NMSE
 [NMSE,NMSSum]=FnNMSECalc(A,AhatNew);
%% norm 2 Error
Norm2Error=norm(A(:,1:size(AhatNew,2))-AhatNew,'fro');
%% Error defined in Underdetermined mixing matrix estimation by exploiting sparsity of
%sources
Error=FnError(A,AhatNew);
