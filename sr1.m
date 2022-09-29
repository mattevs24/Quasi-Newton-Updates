function [Hup, Hupinv] = sr1(H, Hinv, d, y, eta)
    %SR1 update procedure: returns the SR1 updated H and Hinv matrices.

    %Using condition in the pseudocode, decides if to update the matrix
    if norm(d) <= 1e-10 || norm(d'*(y-H*d)) < eta*norm(d)*norm(y-H*d)
        Hup = H;
        Hupinv = Hinv;
    else
        %update Hn here with (3)
        Hup = H + ((y-H*d)*(y-H*d)')/(d'*(y-H*d));
        %update Hinv with (4)
        Hupinv = Hinv + ((d-Hinv*y)*(d-Hinv*y)')/((d-Hinv*y)'*y);
    end
end