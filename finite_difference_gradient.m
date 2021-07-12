%Compute the gradient with the finite difference method

function gradfx=finite_difference_gradient(F,x,n,k,type)
%
%function gradfx=finite_difference_gradient(f,x,n,h,type)
%
%INPUTS:
%F=separated function R^n->R^n, of which we want to find the gradient;
%x=column vector of dimension n;
%n=dimension of x;
%k=integer between 2 and 12 used to find h, the parameter of the finite 
%difference method;
%type='fw' or 'c' to choose the forward/centered finite difference method 
%respectively.
%
%OUTPUTS:
%gradfx=column vector of dimension n which approximates the gradient of f 
%in x.
%

gradfx=zeros(n,1);

h=10^(-k)*norm(x);

switch type
    case 'fw'
        xh=x+h;
        gradfx=(F(xh)-F(x))/h;
    case 'c'
        xh_plus=x+h;
        xh_minus=x-h;
        gradfx=(F(xh_plus)-F(xh_minus))/(2*h);
    otherwise %we do forward one
        xh=x+h;
        gradfx=(F(xh)-F(x))/h;
end
end