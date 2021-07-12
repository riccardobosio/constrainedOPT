clear all;
%
%
%Constrained optimization
%1)Projected gradient method
%
%
%Select a value for n
n=10^6; %We run C1 also with n=10^4

%Set parameters for the projected gradient method
x0=1.5*ones(n,1); %starting point
alpha0=1;
kmax=1000;
tolgrad=1e-12;
c1=1e-4;
rho=0.8;
btmax=50;
gamma=0.1;
tolx=1e-6;

%Create the function which has to be optimized
f=@(x) function_to_optimize(x,n);

%Create the separated function which has to be optimized
F=@(x) separated_function_to_optimize(x,n);

%Create the projection function
Pi_X=@(x) projection(x,ones(n,1),2*ones(n,1));

%Run the projected gradient method using the exact derivatives
%Create the gradient function
gradf=@(x) gradient(x,n);

tic
[xk, fk, gradfk_norm, deltaxk_norm, k, xseq, btseq]=...
    projected_gradient_method(x0, f, gradf, alpha0, kmax, tolgrad, c1,...
    rho, btmax, gamma, tolx, Pi_X);
toc

disp('Number of iterations with exact derivatives:')
disp(k)

disp('Value of the function in xk:')
disp(fk)

%We have to try for different values of h
for j=2:2:12
    j %print j
    
    %Create the gradient computed with the finite difference method
    tic
    fd_gradf_c=@(x) finite_difference_gradient(F,x,n,j,'c'); %centered
    [xk, fk, gradfk_norm, deltaxk_norm, k, xseq, btseq]=...
        projected_gradient_method(x0, f, fd_gradf_c, alpha0, kmax,...
        tolgrad, c1, rho, btmax, gamma, tolx, Pi_X);
    toc
    
    disp('Number of iterations with centered finite differences:')
    disp(k)
    
    disp('Value of the function in xk:')
    disp(fk)
    
    tic
    fd_gradf_fw=@(x) finite_difference_gradient(F,x,n,j,'fw'); %forward
    [xk, fk, gradfk_norm, deltaxk_norm, k, xseq, btseq]=...
        projected_gradient_method(x0, f, fd_gradf_fw, alpha0, kmax,...
        tolgrad, c1, rho, btmax, gamma, tolx, Pi_X);
    toc
    
    disp('Number of iterations with forward finite differences:')
    disp(k)
    
    disp('Value of the function in xk:')
    disp(fk)
end
