function [lse,sexp] = logsumexp(x,minusDeltaBeta,Q)
%LOGSUMEXP  Log-sum-exp function.
%    lse = LOGSUMEXP(x) returns the log-sum-exp function evaluated at 
%    the vector x, defined by lse = log(sum(exp(x)).
%    [lse,sm] = LOGSUMEXP(x) also returns the softmax function evaluated
%    at x, defined by sm = exp(x)/sum(exp(x)).
%    The functions are computed in a way that avoids overflow and 
%    optimizes numerical stability.   

%    Reference:
%    P. Blanchard, D. J. Higham, and N. J. Higham.  
%    Accurately computing the log-sum-exp and softmax functions. 
%    IMA J. Numer. Anal., Advance access, 2020.

if ~isvector(x), error('Input x must be a vector.'), end
if ~isvector(Q), error('Input x must be a vector.'), end

% x = x * (-deltaBeta);
n = length(x);
s = 0; e = zeros(n,1);
[xmax,k] = max(x);;
s = 0;

e = Q.*exp(minusDeltaBeta*(x-xmax));
s = sum(e);
lse = minusDeltaBeta*xmax + log(s);
if nargout > 1
   sexp = exp(lse);
end   

