function [y] = Log_Approx(x, p, epsilon)

In_Idx = (abs(x) <= epsilon);
Out_Idx = ~In_Idx;

y = x;
y(In_Idx) = (x(In_Idx)).^2 ./ (2*epsilon * (p+epsilon) * log(1+1/p));
y(Out_Idx) = ( log(1+abs(x(Out_Idx))./p) - log(1+epsilon/p) + 0.5*epsilon/(p+epsilon) ) ./ log(1+1/p);
end