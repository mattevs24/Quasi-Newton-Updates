function [xd,md] = dogleg(xn, fn, gn, Hn, Hninv, delta)
    %Dogleg implementation: returns the candidate as xd and the model value
    %at the candidate as md.

    %create constants for ease later%
    const1 = gn'*Hn*gn;
    nrm = norm(gn);

    %define xCau as the cauchy point update of xn
    if const1 > (nrm^3)/delta
        alpha = (nrm^2/const1);
    else
        alpha = (delta/nrm);
    end

    %calculate the Cauchy point as well as the value of the model at this
    %point
    xCau = xn - alpha*gn;
    mCau = fn - alpha*nrm^2 + 0.5*alpha^2*const1;
    
    %Start the dogleg algorithm

    %const1 is less than equal zero the newton point and the unidirectional minimiser
    %are not necessarily well defined i.e. use the Cauchy point.
    if const1 <= 0
        %Set the output to be the Cauchy point
        xd = xCau;
        md = mCau;
        %disp('cau1')
        return

    else
        %else we can consider the option of finding the Newton Point and
        %the Unidirectional minimser.

        %First we check if the matrix is invertible and well-conditioned
        dt = det(Hn);
        cnd = cond(Hn);

        %check if ill-conditioned and singular
        if abs(dt) <= 1e-10 || cnd >= 1e9
            %case: Hn singular or ill conditioned -> take the Cauchy Point as output%
            xd = xCau;
            md = mCau;
            %disp('cau2 - ill')
            return
        else
            %case: Hn invertible !and! not ill cond -> continue%
            %Calculate the Newton Step: xNew
            xNew = -Hninv*gn;
            if norm(xNew) <= delta
                %This is the case when the newton step added to xn is still
                %within the radius delta of xn. We take the newton point
                %here as we know that m_n is decreasing along \Gamma with xn+xNew being
                %the limit of this in R_n
                xd = xn+xNew;
                md = fn + gn'*(xd-xn) + 0.5*(xd-xn)'*Hn*(xd-xn);
                %disp('new')
                return
            else
                %we will need to consider the unidirectional minimiser if
                %xNew is outside of the trust region
                %Calculate the unidirectional minimiser step: xUni
                xUni = -(nrm^2/const1)*gn;
                if norm(xUni) >= delta
                    %This is the case where xn + xUni lies outside R_n, in this
                    %case we take the Cauchy step

                    xd = xCau;
                    md = mCau;
                    %disp('cau3 - >=d')
                    return
                else
                    %in this case we have that the unidirectional minimiser
                    %lies within R_n and the newton point lies outside. We then
                    %use the formula from PS5 to calculate the intersection
                    %between \Gamma_n and the trust region boundary
                    gamma = (xNew-xUni)'*(-xUni);
                    
                    %test for gamma to be able to match the conditions in
                    %problem 2
                    if gamma < 0
                        nrm2 = norm(xNew-xUni);
                        const2 = delta^2 - norm(xUni)^2;
                        t_min = (gamma + sqrt(gamma^2 + const2*nrm2^2))/(nrm2^2);

                        xd = xn + xUni + t_min*(xNew-xUni);
                        md = fn + gn'*(xd-xn) + 0.5*(xd-xn)'*Hn*(xd-xn);
                        %disp('xUni')
                        return
                    else
                        xd = xCau;
                        md = mCau;
                        %disp('Cau - else end')
                        return
                    end
                end
            end
        end
    end 
end
