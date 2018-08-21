function [y] = Log_Approx_Weights(x, p, epsilon)

In_Idx = (abs(x) <= epsilon);
Out_Idx = ~In_Idx;

y = x;
y(In_Idx) = 1 ./ (2*epsilon * (p+epsilon) * log(1+1/p));
y(Out_Idx) = 1 ./ ( 2.* abs(x(Out_Idx)) .* (abs(x(Out_Idx)) + p) .* log(1+1/p) );
end