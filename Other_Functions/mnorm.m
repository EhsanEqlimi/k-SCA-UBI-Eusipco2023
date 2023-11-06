function [x,xnrm]=mnorm(x)
xnrm=(max(max(max(abs(x)))));
x=x/(xnrm+eps);