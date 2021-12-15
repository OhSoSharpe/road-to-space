a = 1; % Cylinder radius
U = 1; % Free stream Velocity
N = 250; % Number of source panels
omega = (2*pi/N:2*pi/N:2*pi)*180/pi; % Angle measured clockwise from negative X-axis to mid point of each panel
alpha = 0; % Angle of Attack
beta = 2*pi/N; % Angle made at each Node

for I = 1 : N
    x(I) = round(-cos(((I-1)+1/2)*beta),8); % X component of each Node
    y(I) = round(sin(((I-1)+1/2)*beta),8); % Y component of each Node
end

for I = 1 : N-1
    xm(I) = (x(I)+x(I+1))*(1/2); % X component of each Panel except the last one
    ym(I) = (y(I)+y(I+1))*(1/2); % Y component of each Panel except the last one
    theta(I) = atan2d(y(I+1)-y(I),(x(I+1)-x(I))); % Angle from each midpoint to the next Node except the last to first
end
    xm(N) = (x(1)+x(N))*(1/2); % X component of last Panel
    ym(N) = (y(1)+y(N))*(1/2); % Y component of last Panel
    theta(N) = atan2d(y(1)-y(N),x(1)-x(N)); % Angle from last midpoint to first Node

for J = 1:N-1
    deltaS(J) = sqrt((y(J+1)-y(J))^2 + (x(J+1)-x(J))^2); % The distance from one Node to the next exept the last to the first
end
    deltaS(N) = sqrt((y(1)-y(N))^2 + (x(1)-x(N))^2); % The distance from the last Node to the first

for I = 1:N
    for J = 1:N
        if I == J
            C(I,J) = 0.5; % If I equals J then C must be equal to 0.5
        else
            R(I,J) = sqrt((xm(I) - xm(J))^2 + (ym(I)-ym(J))^2); % The distance from the mid point of the J panel to the midpoint of the I panel
            phi(I,J) = atan2d(ym(J)-ym(I),xm(J)-xm(I)); % The angle from the J panel to the I panel
            C(I,J) = ((sind(theta(I)-phi(I,J)))/(2*pi*R(I,J)))*deltaS(J); 
            Cbar(I,J) = ((cosd(phi(I,J)-theta(I)))/(2*pi*R(I,J)))*deltaS(J);
        end
    end
end

for I = 1:N-1
    n(I) = (ym(I+1)-ym(I))/deltaS(J)+(xm(I+1)-xm(I))/deltaS(J); % Unit normal vector directed outward from each panel except the last one
end
    n(N) = (ym(1)-ym(N))/deltaS(J)+(xm(1)-xm(N))/deltaS(J); % Unit normal vector directed outward from the last panel
    
    lambda = round(linsolve(C,n'),8); % The source strength of each panel

    S = Cbar.*lambda'; %Multiplying Cbar and lambda together in order to sum it and find tangential velocity
    F = sum(S,2); %Second part of Tangential Velocity equation, Summation of Cbar * lambda
    L = cosd(theta)'; % First part of Tangential Velocity equation, V0 * cosd(theta-alpha)
    Vt = L - F; % Tangential Velocity at each panel


    k1 = -2*sind(theta)'; % First part of new Tangential Velocity
    Kvalue = 0; % K value for determining source strength
    q = Kvalue * 2 * pi ; % Determining source strength from K value
    Vt = k1 + q/(2*pi); % New Tangential Velocity for plots
    theta = theta'; % Making my Theta's into a single column for easier transfer of data to excel

for I = 1:N
    Cp = 1 - Vt.^2; % Pressure Coefficient of each panel
end


% The following code is used to plot the pressure coefficients against theta

lonWrapped = wrapTo360(theta); % Adjusts the plot so it starts at 0 and ends at 260 instead of gfoing from -180 to 180
plot(lonWrapped,Cp,'.'); %Plots the Pressure Coefficient against the theta
grid on
xlabel('Theta (in degrees)');
ylabel('Pressure Coefficient');
title('Cp vs. Theta,  K = 0, -1, -2, -3')
legend('K = 0', 'K =-1', 'K = -2', 'K = -3')
theta = theta;


%Following is code used to plot Streamline Patterns

syms r(x,y) Theta(x,y) f(x,y)
r(x,y) = sqrt(x^2+y^2);
A = -1.;
Theta(x,y) = atan(y,x);
f(x,y) = (r-1/r)*sin(Theta)+ A*log(r);
f(x,y) = (r-1/r)*sin(Theta);
fc = fcontour(f,[-12 12 -12 12]);
fc.LevelList = [-2,-1, -0.75,-0.65, -0.5, -0.25, 0., 0.25, 0.5,0.65, 0.75, 1.0, 2.0];
fc.LineColor = 'k'
axis equal;
title('Streamline Pattern,  K = -3')





