%% Param Portret Sistemy (2)

mu = linspace(0, 10, 100);

nu = 4.*(mu.^2+3.*mu+3+1./mu);

plot(mu,nu);


str = 'Ã';
text(7, 270,str)

str = '\nu > 4(\mu^2 + 3 \mu + 3 + 1 / \mu)';
text(2, 370,str)


str = '\nu < 4(\mu^2 + 3 \mu + 3 + 1 / \mu)';
text(6, 170,str)


%% task 2 (1,3)
clc
clear 

% u = u(t), v = v(t)
%du == mu * nu * u^2 /(1+u) - u*v
%dv == v*(-mu+u)

f = @(mu,nu,u,v) mu * nu * u^2 /(1+u) - u*v;
g = @(mu,u,v) v*(-mu+u);

% task 1
%mu = 0.1667;
%nu = 0.5;

% task 3
%mu = 1;
%nu = 33;

mu = 2;
nu = 20 ;

tspan = linspace(0,100,1000);

hold on


for i = linspace(0,10,30)
    i
    for j = linspace(0,30,30)
    
        
        % task 1
        % y0 = [ 0.1*i 0.0001 ];

        % task 3
        %y0 = [1.01 15.501];

        % task 4
        y0 = [ i  j  ];
        


        [t,y] = ode45( @(t,y) func(t,y,mu,nu) , tspan, y0);


            %du = mu .* nu .* y(:,1).^2 /(1+y(:,1)) - y(:,1).*y(:,2);
            %dv = y(:,2).*(-mu+y(:,1));


        plot(y(:,1),y(:,2), 'r' )


        %quiver
        
        for k = 1:20
            k = k*10;
            quiver(y(k,1), y(k,2), y(k+1,1) - y(k,1), y(k+1,2) - y(k,2), 1, 'k', 'MaxHeadSize', 10 );
        end
        

    
        axis([0 10 0 30])
    end
end


hold off


function dydt = func(t,y,mu,nu)

dydt = [mu * nu * y(1)^2 /(1+y(1)) - y(1)*y(2); y(2)*(-mu+y(1))];
end