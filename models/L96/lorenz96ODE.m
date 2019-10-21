function dx=lorenz96ODE(x,n,F);
% function dx=lorenz96ODE(x,n,F);
% This rountine is developed to generate Lorenz 96 attractor
%%%%By Xiaodong Luo, June 2007%%%%

for i=1:n,    
    dx(i,1) = (x(ind(i+1,n),1) - x(ind(i-2,n),1))*x(ind(i-1,n),1) - x(i,1) + F;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subfunction to sort out the true index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function trInd=ind(i,n);
        
        if i==0 | i==n,
            trInd=n;
        else
            trInd=mod(i,n);
        end
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
