%Function to project xk_bar into the feasible set
%
%
%
function x_hat=projection(x,lower,upper)
%
%function x_hat=projection(x,lower,upper)
%
%INPUTS:
%x is a column vector containing n components
%lower is a column vector where the i-th element is the lower bound of the
%i-th dimension of the feasible set
%lower is a column vector where the i-th element is the upper bound of the
%i-th dimension of the feasible set
%
%OUTPUT:
%x_hat is equal to x if x is already in the feasible set otherwise it
%is updated to the projection of x in theboundary of the feasible set
%
x_hat=max(min(x,upper),lower);
end