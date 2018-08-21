function [y] = General_Approx_Weights(x, p, epsilon, approx)

if strcmp(approx, 'Lp')
    [y] = Lp_Approx_Weights(x, p, epsilon);
elseif strcmp(approx, 'Log')
    [y] = Log_Approx_Weights(x, p, epsilon);
elseif strcmp(approx, 'Exp')
    [y] = Exp_Approx_Weights(x, p, epsilon);
else
    disp('The approximation method should be one of Lp, Log, and Exp!')
end

end