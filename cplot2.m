%
% cplot(x) - plot complex function;
%

%  written by John Pauly, 1989
%  (c) Board of Trustees, Leland Stanford Junior University

%fadilali modified (renamed to "cplot2") to include an x-axis dimension.

function  cplot2(xaxis, x)

l = length(x);
t = [1:l]/(l+1);
plot(xaxis,real(x),xaxis,imag(x), 'LineWidth',5.0);






