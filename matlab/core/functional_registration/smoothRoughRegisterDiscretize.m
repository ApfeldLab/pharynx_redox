function [i1_reg,i2_reg] = smoothRoughRegisterDiscretize(i1, i2, smoothLambda, roughLambda, warpLambda, nSamples)
%SMOOTHROUGHREGISTERDISCRETIZE Summary of this function goes here
%   Detailed explanation goes here
[roughFD1, regRoughFD2] = smoothRoughRegister(i1, i2, smoothLambda, roughLambda, warpLambda);

xs = linspace(1,100,nSamples);
i1_reg = eval_fd(xs, roughFD1);
i2_reg = eval_fd(xs, regRoughFD2);
end