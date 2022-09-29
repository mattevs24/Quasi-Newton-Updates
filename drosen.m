function z = drosen(x)
    z(1) = 40*x(1)^3 + (2-40*x(2))*x(1) - 2;
    z(2) = 20*(x(2)-x(1)^2);
    z = z';
end