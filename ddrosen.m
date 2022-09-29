function z = ddrosen(x)
    
    z(1,1) = 2+120*x(1)^2-40*x(2);
    z(1,2) = -40*x(1);
    z(2,1) = z(1,2);
    z(2,2) = 20;

end