function visual(f,df,d2f,x,x0,xex)
% Postprocessing: Plot contour lines for the objective function f and the 
% solution path of a particular sequence of iterates.

% Specify the window in parameter space that we want to plot.

x1 = [min(-1.5,min(x(1,:))):0.01:max(1.5,max(x(1,:)))];
x2 = [min(-0.5,min(x(2,:))):0.01:max(1.5,max(x(2,:)))];
%x1 = [-1.5:0.01:1.5];
%x2 = [-1.5:0.01:1.5];

% Plot the contours of the objective function

Z=zeros(length(x2),length(x1));
for i=1:length(x1)
    for j=1:length(x2)
        Z(j,i) = f([x1(i);x2(j)]);
    end
end
contour(x1,x2,Z,40,'LineWidth',2);
hold on

% Plot the exact solution and the iterates

plot(xex(1),xex(2),'b*');
plot(x(1,:),x(2,:),'r-');

% Plot the initial guess and the steepest descent direction 

plot(x0(1),x0(2),'md');

f0 = f(x0)
g0 = df(x0)
H0 = d2f(x0)

xN1 = x0 - H0\g0

% ghat = g0/norm(g0);
% plot([x0(1) x0(1)-ghat(1)],[x0(2) x0(2)-ghat(2)],'g-');

% Evaluate the quadratic model at x0 and plot it 

Z=zeros(length(x2),length(x1));
for i=1:length(x1)
    for j=1:length(x2)
        h = [x1(i);x2(j)]-x0;
        Z(j,i) = f0 + g0'*h + 0.5*h'*H0*h;
    end
end

fmin = min(min(Z(1,:)));
step=((2*f0)^0.5-fmin^0.5)/10;
zlev = [fmin^0.5:step:(2*f0)^0.5].^2;
contour(x1,x2,Z,10,'LineColor','k','LineWidth',0.5,'LevelList',zlev);

hold off

end

