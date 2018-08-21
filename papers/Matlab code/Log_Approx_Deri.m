function [y] = Log_Approx_Deri(x, p, epsilon)

In_Idx = (abs(x) <= epsilon);
Out_Idx = ~In_Idx;

y = x;
y(In_Idx) = (x(In_Idx)) ./ (epsilon * (p+epsilon) * log(1+1/p));
y(Out_Idx) = sign(x(Out_Idx)) ./ ( (p + abs(x(Out_Idx))) .* log(1+1/p) );
end