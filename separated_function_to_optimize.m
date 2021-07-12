function Y = separated_function_to_optimize(x,n)
%
%INPUTS:
%x=column vector of length n;
%n=number of dimensions of x;
%OUTPUTS:
%Y=separated function R^n->R^n.
%
Y=zeros(n,1);
for i=1:n
    Y(i)=1/4*x(i)^4+1/2*x(i)^2-x(i);
end

end