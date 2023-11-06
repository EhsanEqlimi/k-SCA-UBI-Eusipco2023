function FnMixingCalc_FMN2(sigma,L)
% Main Loop
while sigma>sigma_min
    for i=1:L
        delta = OurDelta(s,sigma);
        s = s - mu_0*delta;
        s = s - A_pinv*(A*s-x);   % Projection
    end
    
    if ShowProgress
        fprintf('     sigma=%f, SNR=%f\n',sigma,estimate_SNR(s,true_s))
    end
    
    sigma = sigma * sigma_decrease_factor;
end