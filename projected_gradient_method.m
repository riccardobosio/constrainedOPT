%Projected gradient method function
function [xk, fk, gradfk_norm, deltaxk_norm, k, xseq, btseq]=...
    projected_gradient_method(x0, f, gradf, alpha0, kmax, tolgrad, c1,...
    rho, btmax, gamma, tolx, Pi_X)
%
%function [xk, fk, gradfk_norm, deltaxk_norm, k, xseq, btseq]=
%projected_gradient_method(x0, f, gradf, alpha0, kmax, tolgrad,
%c1, rho, btmax, gamma, tolx, Pi_X)
%
%INPUTS:
%x0=n-dimensional column vector;
%f=function R^n->R;
%gradf=function that describes the gradient of f;
%alpha0=first constant that multiplies the descent direction;
%kmax=maximum number of iterations;
%tolgrad=tolerance for the gradient (used as a stopping criterion);
%c1=Armijo condition factor that must be a scalar in (0,1);
%rho=fixed factor lesser than one used for reducing alpha0 at each
%iteration;
%btmax=maximum number of steps of the backtracking strategy to update
%alpha;
%gamma=constant that multiplies the descent direction before projection;
%tolx=tolerance for the norm of the difference between two consecutive 
%xk (used as a stopping criterion);
%Pi_X=function that handles the projection;
%
%OUTPUTS:
%xk=the last x computed by the function;
%fk=f(xk);
%gradfk_norm=norm of gradf(xk);
%k=last iteration performed;
%xseq=n-by-k matrix where the columns are the xk computed during the 
%iterations;
%btseq=1-by-k vector where elements are the number of backtracking
%iterations at each optimization step.
%

%Initializations
xseq=zeros(length(x0),kmax);
btseq=zeros(1,kmax);
xk=Pi_X(x0);
fk=f(xk);
k=0;
gradfk_norm=norm(gradf(xk));
deltaxk_norm=tolx+1; %to be sure of entering the while

%create function to handle armijo conditions
f_armijo=@(fk,alpha,xk,pk) fk+c1*alpha*gradf(xk)'*pk;

while k<kmax && gradfk_norm>=tolgrad && deltaxk_norm>=tolx
    %Find the descent direction
    pk=-gradf(xk);
    xk_bar=xk+gamma*pk;
    xk_hat=Pi_X(xk_bar);
    
    %Reset alpha
    alpha=alpha0;
    
    %Compute the new xk
    pik=xk_hat-xk;
    xnew=xk+alpha*pik;
    
    %Compute f in xnew
    fnew=f(xnew);
    
    bt=0;
    
    %Start backtracking
    while bt<btmax && fnew>f_armijo(fk,alpha,xk,pik)
        %Reduce alpha
        alpha=rho*alpha;
        
        %Update fnew and xnew
        xnew=xk+alpha*pik;
        fnew=f(xnew);
        
        bt=bt+1;
    end
    
    %Update xk, fk, gradfk_norm, deltaxk_norm
    deltaxk_norm=norm(xnew-xk);
    xk=xnew;
    fk=fnew;
    gradfk_norm=norm(gradf(xk));
    
    %Increase iteration
    k=k+1;
    
    %Store xk in k-th column of xseq
    xseq(:,k)=xk;
    
    %Store bt in btseq
    btseq(k)=bt;
end

%Resize xseq and btseq
xseq=xseq(:,1:k);
btseq=btseq(1:k);

end