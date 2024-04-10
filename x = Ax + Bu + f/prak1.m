% preparation
clc
clear

% main input
N_psi = 10;
N_xy0 = 300;
N_xy1 = 300;
N_t   = 2000;

A = [ 1 0;
     23  1];
B = [1 0;
     -1 2];
f = [3;
     0];
t0 = 0;
t1 = 10;
x_0 = [0;
       1];
   
t = linspace(t0, t1, N_t);
   
% P  set 
r = 1.4;% radius
p = [1; % center
     1];
 
% X0 set
alpha = 1;
gamma = 0.4;

% X1 set
a = 6;

% starting enumeration of psi_0
delta = linspace(0, 2*pi, N_psi)';
psi_0 = [cos(delta) sin(delta)]';

% getting vector of couples [psi(t)_1, psi(t)_2, ...]
% [ x1 x2 ...
% [ y1 y2 ...
psi_equ = @(t, psi_0) expm(-A'.*(t-t0))*psi_0;


% Start drawing 
axis equal
hold on;

% Draw X0
y1 = gamma^(1/4) + 1;
y0 = -gamma^(1/4) + 1;
x0 = -(gamma/alpha)^(1/2);
x1 = (gamma/alpha)^(1/2);

x = linspace(x0, x1, N_xy0);
y = linspace(y0, y1, N_xy0);
[X, Y] = meshgrid(x, y);
Z = alpha.*(X).^2 + (Y - 1).^4;
contour(X, Y, Z, [gamma, gamma], 'b-');

% Draw X1
x2 = 1/2*(a + sqrt(a^2 - 4));
x1 = 1/2*(a - sqrt(a^2 - 4));

xlin = linspace(x1, x2, N_xy1);
y1 = xlin.^(-1);
y2 = a - xlin;
plot(xlin, y1, 'r');
plot(xlin, y2, 'r');




for i = 1:N_psi
    for j = 1:N_t%for time
        ans = expm(-A'.*(t(j)-t0))*psi_0(:,i);
       
        psi_x(j) = ans(1);
        psi_y(j) = ans(2);
        
        u_x(j) = r * psi_x(j)/sqrt(psi_x(j)^2+psi_y(j)^2)+p(1);
        u_y(j) = r * psi_y(j)/sqrt(psi_x(j)^2+psi_y(j)^2)+p(2);
        
        
        
    
    end
    sin = psi_0(2,i);
    cos = psi_0(1,i);
    r = sqrt( -alpha*cos^2 + sqrt(alpha^2*cos^4 + 4*gamma*sin^4) )/sqrt(2*sin^4);
        
    
    x0 = [x1(i),x2(i)]'
    
    100*i/N_psi
    
    x02(i) = 1 + psi_0(2,i).*r ;
    x01(i) = 0 + psi_0(1,i).*r ;
    plot(x01(i), x02(i), '*');
    for j = 1:N_t
        x1(j+1) = h*(A(1,1)*x1(j) + A(1,2)*x2(i,j) + B(1,1)*u_x(i,j) + B(1,2)*u_y(i,j) + f(1)) + x1(i,j);
        x2(j+1) = h*(A(2,1)*x1(j) + A(2,2)*x2(i,j) + B(2,1)*u_x(i,j) + B(2,2)*u_y(i,j) + f(2)) + x2(i,j);
    end
    plot(x1(),x2());
    
    %opts = odeset('RelTol', relTol, 'AbsTol', absTol, 'Events', @(t,x)aimReached(t, x, x1, r1, Q, p));
    %[t,x] = ode45(@(t, x) odefun(t, x, t, A, B, u_x, u_y, f) , [t0 t1], x0  );
    %plot(x(:, 1), x(:, 2), 'b');
end

function dydt = myode(t,y,ft,f,gt,g)
f = interp1(ft,f,t); % Interpolate the data set (ft,f) at time t
g = interp1(gt,g,t); % Interpolate the data set (gt,g) at time t
dydt = -f.*y + g; % Evaluate ODE at time t
end

function dxdt = odefun(t, x, old_t, A, B, u_x, u_y, f)
    
    ux = interp1(old_t, u_x, t);
    uy = interp1(old_t, u_y, t);
        
    dxdt = A * x + B * [ux uy]' + f;
end


function [value,isterminal,direction] = aimReached(t, x, x1, r1, Q, p)
   value = inSet(x1, r1, Q, x, p);
   isterminal = 1;
   direction = 0;
end

function val = inSet(x1, r1, Q, x, p)
    options = optimoptions(@fminunc, 'Display', 'off');
    [x,fval,exitflag,output] = fminunc(@(y) -fun(x1, r1, Q, p, x, y), x, options);
    if exitflag ~= -3
        val = 1;
    else
        val = 0;
    end
end

function val = fun(x1, r1, Q, p, x2, x)
    val =  dot(x2, x) - dot(x1, x) - r1/2 * norm(x, p/(p - 1)) - norm(Q * x);
end
%{

% Find Candidate u
for i = 1:N_psi
    for j = 1:N_t
        u_x(i,j) = r * psi_x(i,j)/sqrt(psi_x(i,j)^2+psi_y(i,j)^2)+p(1);
        u_y(i,j) = r * psi_y(i,j)/sqrt(psi_x(i,j)^2+psi_y(i,j)^2)+p(2);
    end
    disp(i)
end

u = @(i,t) [r*psi_x(i,find(t))/sqrt(psi_x(i,find(t))^2+psi_y(i,find(t))^2)+p(1);
        r*psi_y(i,find(t))/sqrt(psi_x(i,find(t))^2+psi_y(i,find(t))^2)+p(2)]

    
% Solve equations
x1 = zeros(N_psi, N_t);
x2 = zeros(N_psi, N_t);
for i = 1:N_psi
    sin = psi_0(1,i);
    cos = psi_0(2,i);
    r = sqrt( -alpha*cos^2 + sqrt(alpha^2*cos^4 + 4*gamma*sin^4) )/sqrt(2*sin^4);
    x2(1,i) = 1 + psi_x(i,1).*r ;
    x1(1,i) = 0 + psi_y(i,1).*r ;
    plot(x1(1,i),x2(1,i),'*');
end

h = (t_1-t_0)/N_t;
for i = 1:N_psi-1
    
    for j = 1:N_t
        x1(i+1,j) = h*(A(1,1)*x1(i,j)+A(1,2)*x2(i,j) + B(1,1)*u_x(i,j) + B(1,2)*u_y(i,j) + f(1)) + x1(i,j);
        x2(i+1,j) = h*(A(2,1)*x1(i,j)+A(2,2)*x2(i,j) + B(2,1)*u_x(i,j) + B(2,2)*u_y(i,j) + f(2)) + x2(i,j);
    end
    
    disp(i)
end
for i = 1:N_psi
    plot(x1(i,:), x2(i,:))
end
%}

function res = support_f_X1(psi, a) %support_f_X1 ([psi_x,psi_y],a)
    if(a == 2)
        res = psi(1)+psi(2);
    else
        res = max( (a/2)*(psi(1)+psi(2)) + abs(psi(1)-psi(2))*sqrt(a^2-4)/2 , -2*sqrt( psi(1)*psi(2) ));
    end
end

 