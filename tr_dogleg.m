function [x,n] = tr_dogleg(f,df,H0,xn,delta,delmax,rho_ac,tol,eta)
    %tr_dogleg Implementation: returns the approximate solution xn as x
    %and the number of iterations as n.

    %The test in 6b was called as below:

    %[x,n] = tr_dogleg(@rosen,@drosen,eye(2),[-1.2;1],0.2,1,0.125,1e-5,1e-5)
    
    %then the visual function was then called:

    %visual(@rosen,@drosen,@ddrosen,x,[-1.2;1],[1;1])

    %this was repeated for the other test functions with appropriate random
    %x inputs.

    
    %initialise the first iteration
    Hn = H0;
    %Hninv = Hn\eye(length(xn)); wanted to try use \ method but turned out
    %an order of magnitude slower
    Hninv = inv(Hn); %need to make this step better particulary with harder matrices
    n = 0;
    gn = df(xn);

    %perform a while loop for when the norm of the derivative is greater
    %than some tolerance
    while norm(gn)>=tol

        %initialise the entry in the x vector (used to be able to visualise
        %tr_dogleg later)

        %x(:,n+1) = xn;
        %ignoring the warning here due to the total number of iterations
        %being unknown before the algorithm starts

        %set the new values for the function and its derivative
        gn = df(xn);
        fn = f(xn);

        %perform the dogleg function here
        [xd,md] = dogleg(xn,fn,gn,Hn,Hninv,delta);
        
        %set the values dn and yn for use in the SR1 update
        ddn = xn;
        dn = xd-xn;
        yn = df(xd)-gn;


        %calculate rho_n to be able to confirm if the dogleg will be
        %accepted
        rhon = (fn-f(xd))/(fn-md);

        %acceptance of xd as the next iteration of xn
        if rhon >= rho_ac
            xn = xd;
        end
        %section 6.2 testing
        %nrm(n+1) = norm(xn-[1;1]);

        %adjust the trust region radius according to alg 5.2
        if rhon<0.25
            delta = 0.25*delta;
        elseif ((rhon > 0.75) && (abs(norm(ddn-xn)-delta)<= 1e-12))
            delta = min(2*delta,delmax);
        end

        %update the Hn and Hninv matrices by the sr1 method
        [Hn,Hninv] = sr1(Hn,Hninv,dn,yn,eta);

        %update the interation count
        n = n+1;
    end
    %remove this and uncomment the x(:,n) section to be able to return the
    %vector for visualisation in the visual function.
    x = xn;
end