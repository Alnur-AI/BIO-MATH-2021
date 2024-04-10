%% task 2 
clc
clear 

% u = u(t), v = v(t)
%du == mu * nu * u^2 /(1+u) - u*v
%dv == v*(-mu+u)

f = @(mu,nu,u,v) mu * nu * u^2 /(1+u) - u*v;
g = @(mu,u,v) v*(-mu+u);

mu = 0.1667;
nu = 0.5;

tspan = linspace(0,1,100);
y0 = [0 0.01];
[t,u] = ode45( [ f(mu,nu,u,v)  g(mu,u,v) ] , tspan, y0);

%% NS bifurcation
clc
clear

f = @(r,u,v) r.*u.*(4-v.^2);

n = 5000;
u = zeros(n,1);
v = zeros(n,1);


r = 0.45 ;
u(1) = 1.0;
v(1) = 1.0;

hold on
for i = 1:n-1
    v(i+1) = u(i);
    u(i+1) = f(r,u(i),v(i));
    
    i/(n-1)*100
end
plot(u,v,'k*');
plot(2/sqrt(3),2/sqrt(3),'r*');
%axis([0.9 1.35 0.9 1.35])
axis([0 2 0 2])
title( strcat ('r =  ' , num2str ( r ), ', u_0 = ', num2str ( u(1) ), ', v_0 = ', num2str ( v(1) ) ) );
xlabel('u');
ylabel('v');
hold off


%% bifurcation

clc
clear
f = @(r,u) r.*u.*(4-u.^2);



k_1 = 100;
k_2 = 5000;
k_3 = 4900;

u_0 = 1;
r_0 = 0.6;
r_delta = 0.001;

plot(0,2);
hold on
for i = 1:k_1
    x = zeros(k_2);
    t = u_0;
    for j = 1:k_2
        q = f(r_0,t);
        t = q;
        x(j) = t;
        
        fprintf( "%4.2f%% %4.2f%%\n" ,  (i)/k_1*100, (j)/k_2*100  )
        
    end
    y = x(k_2-k_3 : k_2);
    A = repmat(r_0,k_3+1,1);
    plot(A,y,'ko','MarkerSize',1);
    axis([0.6 0.68 0 2])
    
    
    
    
    
    r_0=r_0 + r_delta;
end
xlabel('r') 
ylabel('u_t') 

%% f(r,u) plot animation
clc
clear 

f = @(r,u) r.*u.*(4-u.^2);

% r parameters
r_count = 10;
r_left  = 0;
r_right = 1.5;

% x parameters
x_count = 200;
x_left  = -3;
x_right = 3;

% plot
r = linspace (r_left, r_right, r_count) ;
x = linspace (x_left, x_right, x_count);
y = f(r(:),x);


%// Plot starts here
figure,hold on

%// Set x and y limits of the plot
xlim([0 2])
ylim([0 2])

%// Plot point by point
for k = 1:r_count
    plot(x(:),y(k,:)) 
    pause(0.01);     
end
%% cycle 2
clc
clear


%r = 0.649519052838;
r = 0.000;
syms temp u v ff(u) dff(u)
hold on

answer = zeros(651,10);

for i = 1:650
        
    r = 0.001*i
    answer(i+1,1) = r;
    
    
    % Syms
    fu = r*u*(4-u^2);
    ffu = r*fu*(4-fu^2);
    dffu = diff(ffu);


    equ = [ffu - u-v  == 0   dffu - 1-v  == 0  ];
    S = solve( equ, [u,v]);
    
    Sdouble = double(S.u);
    Sv = double(S.v);
    for j = 1:size(Sdouble)
        if( 0.0 < Sdouble(j) && Sdouble(j) < 2.0 && imag( Sdouble(j) ) == 0 )
            fprintf( "root #%d %4.2f\n" ,j ,Sdouble(j)  )
            
            answer(i+1,j+1) = Sdouble(j);
            
        end
    end
    
    % Functions
    ff(u) = ffu;
    dff(u) = dffu;

    x = linspace(0,2,100);
    yyy = ff(x);
    dyyy = dff(x);
    
    
    plot(x,yyy-x,x,dyyy-1, x, zeros(1,100));
    pause(0.001)
    hold off


end

%% cycl 3 ideal plot see it  arrows
f = @(r,u) r.*u.*(4-u.^2);

figure
hold on
    
    
    %// Progress
    r= 0.6126
    
    %// x parameters
    x_count = 200;
    x_left  = 0;   
    x_right = 2;

    %// Plot
    x = linspace (x_left, x_right, x_count);
    y45 = x;
    y = f(r, x);
    
    %// Calm dots
    cp1 = 0;
    cp2 = sqrt( (4.*r - 1)./r );

    %// Recur 
    steps = 20;
    u_01 = 0.3;
    for i = 1:steps
        u_rec1( i, 1 ) = f_rec( i-1 , u_01, r);
    end
    
    u_02 = 1.6;
    for i = 1:steps
        u_rec2( i, 1 ) = f_rec( i-1 , u_02, r);
    end

    %// Plot starts here
    

    %// Set x and y limits of the plot
    xlim([0 2])
    ylim([0 2])

    %// label
    title( strcat ('r =  ' , num2str ( r ), ', u_0 = ', num2str ( u_01 ) ) );
    xlabel('u');
    ylabel('f(r,u)');

    %// Plot u and f(r,u) 
    plot(x, y(1,:)) 
    plot(x, y45(1,:))

    %// Points and arrows
    for k = 3:steps
        % part 1
        p0 = [u_rec1(k-2) u_rec1(k-1)];
        p1 = [u_rec1(k-1) u_rec1(k-1)];

        quiver( p0(1), p0(2), p1(1)-p0(1), p1(2)-p1(1), 0 , 'k'   );

        plot(p0(1), p0(2), 'k*');
        plot(p1(1), p1(2), 'k*');

        % part 2
        p0 = [u_rec1(k-1) u_rec1(k-1)];
        p1 = [u_rec1(k-1) u_rec1(k)];

        quiver( p0(1), p0(2), p1(1)-p0(1), p1(2)-p1(1), 0 , 'k'   );

        plot(p0(1), p0(2), 'k*');
        plot(p1(1), p1(2), 'k*');
        
        pause(0.5)
        
    end
    
    %// calm dots
    plot(cp1, f(r,cp1), 'c*') 
    plot(cp2, f(r,cp2), 'c*')
    
    hold off


%% cycle 3 plot
clc
clear


r = 0;
hold on
syms temp u v fff(u) dfff(u)
    
for i = 1:1
    
    %r = 0.6126134
    r = 0.614
    
    % Syms
    fu = r*u*(4-u^2);
    ffu = r*fu*(4-fu^2);
    fffu  = r*ffu*(4-ffu^2);
    dfffu = diff(fffu);

    
    
    % Functions
    fff(u) = fffu;
    dfff(u) = dfffu;

    x = linspace(0,2,1000);
    yyy = fff(x);
    dyyy = dfff(x);
    
    
    
    
    

    

    plot(x,yyy-x,    x,dyyy-1,   x, zeros(1,1000));
    hold off
    

    
end    

    

%% cycle 3
clc
clear




%r = 0.649519052838;
r = 0.000;
syms temp u v fff(u) dfff(u)
hold on

answer = zeros(651,10);

for i = 612:612
        
    r = 0.001*i
    answer(i+1,1) = r;
    
    
    % Syms
    fu = r*u*(4-u^2);
    ffu = r*fu*(4-fu^2);
    fffu  = r*ffu*(4-ffu^2);
    dfffu = diff(fffu);

    %{
    equ = [fffu - u-v  == 0   dfffu - 1-v  == 0  ];
    S = solve( equ, [u,v]);
    
    Sdouble = double(S.u)
    Sv = double(S.v);
    for j = 1:size(Sdouble)
        if( 0.0 < Sdouble(j) && Sdouble(j) < 2.0 && imag( Sdouble(j) ) == 0 )
            fprintf( "root #%d %4.2f\n" ,j ,Sdouble(j)  )
            
            answer(i+1,j+1) = Sdouble(j);
            
        end
    end
    %}
    
    
    % Functions
    fff(u) = fffu;
    dfff(u) = dfffu;

    x = linspace(0,2,100);
    yyy = fff(x);
    dyyy = dfff(x);
    
    
    plot(x,yyy-x,x,dyyy-1, x, zeros(1,100));
    pause(0.001)
    hold off


end
%% u_t plot arrows animation
f = @(r,u) r.*u.*(4-u.^2);


for r = 0.5


    
    
    %// Progress
    r
    
    %// x parameters
    x_count = 200;
    x_left  = 0;   
    x_right = 2;

    %// Plot
    x = linspace (x_left, x_right, x_count);
    y45 = x;
    y = f(r, x);
    
    %// Calm dots
    cp1 = 0;
    cp2 = sqrt( (4.*r - 1)./r );

    %// Recur 
    steps = 20;
    u_01 = 1.351;
    for i = 1:steps
        u_rec1( i, 1 ) = f_rec( i-1 , u_01, r);
    end
    
    u_02 = 1.6;
    for i = 1:steps
        u_rec2( i, 1 ) = f_rec( i-1 , u_02, r);
    end

    %// Plot starts here
    figure
    hold on

    %// Set x and y limits of the plot
    xlim([0 2])
    ylim([0 2])

    %// label
    title( strcat ('r =  ' , num2str ( r ), ', u_0 = ', num2str ( u_01 ) ) );
    xlabel('u');
    ylabel('f(r,u)');

    %// Plot u and f(r,u) 
    plot(x, y(1,:)) 
    plot(x, y45(1,:))

    %// Points and arrows
    for k = 3:steps
        % part 1
        p0 = [u_rec1(k-2) u_rec1(k-1)];
        p1 = [u_rec1(k-1) u_rec1(k-1)];

        quiver( p0(1), p0(2), p1(1)-p0(1), p1(2)-p1(1), 0 , 'k'   );

        plot(p0(1), p0(2), 'k*');
        plot(p1(1), p1(2), 'k*');

        % part 2
        p0 = [u_rec1(k-1) u_rec1(k-1)];
        p1 = [u_rec1(k-1) u_rec1(k)];

        quiver( p0(1), p0(2), p1(1)-p0(1), p1(2)-p1(1), 0 , 'k'   );

        plot(p0(1), p0(2), 'k*');
        plot(p1(1), p1(2), 'k*');
        
        %{
        % part 1
        p0 = [u_rec2(k-2) u_rec2(k-1)];
        p1 = [u_rec2(k-1) u_rec2(k-1)];

        quiver( p0(1), p0(2), p1(1)-p0(1), p1(2)-p1(1), 0 , 'r'   );

        plot(p0(1), p0(2), 'r*');
        plot(p1(1), p1(2), 'r*');

        % part 2
        p0 = [u_rec2(k-1) u_rec2(k-1)];
        p1 = [u_rec2(k-1) u_rec2(k)];

        quiver( p0(1), p0(2), p1(1)-p0(1), p1(2)-p1(1), 0 , 'r'   );

        plot(p0(1), p0(2), 'r*');
        plot(p1(1), p1(2), 'r*');
        %}
    end
    
    %// calm dots
    plot(cp1, f(r,cp1), 'c*') 
    plot(cp2, f(r,cp2), 'c*')

    %// save and close
    %saveas( gcf, strcat(num2str ( r , '%.4f'),'.png') )
    %close 
end

%% u_t plot arrows
clc
clear 

f = @(r,u) r.*u.*(4-u.^2);

% x parameters
x_count = 200;
x_left  = 0;   
x_right = 2;

% plot
r = 0.6157; % max 0.65
x = linspace (x_left, x_right, x_count);
y45 = x;
y = f(r, x);

%// Calm dots
cp1 = 0;
cp2 = sqrt( (4.*r - 1)./r );

% recur
steps = 50;
u_0 = 0.91;

u_rec( 1 ) = u_0;
for i = 2:steps
    u_rec( i ) = f( r , u_rec(i-1));
end





%// Plot starts here
figure
hold on

%// Set x and y limits of the plot
xlim([0 2])
ylim([0 2])

% label
title( strcat ('r =  ' , num2str ( r ), ', u_0 = ', num2str ( u_0 ) ) );
xlabel('u');
ylabel('f(r,u)');

%// Plot point by point
plot(x, y(1,:)) 
plot(x, y45(1,:))

%// Points and arrows
for k = 3:steps
    p0 = [u_rec(k-2) u_rec(k-1)];
    p1 = [u_rec(k-1) u_rec(k-1)];

    quiver( p0(1), p0(2), p1(1)-p0(1), p1(2)-p1(1), 0 , 'k'   );

    plot(p0(1), p0(2), '*');
    plot(p1(1), p1(2), '*');
    
    p0 = [u_rec(k-1) u_rec(k-1)];
    p1 = [u_rec(k-1) u_rec(k)];

    quiver( p0(1), p0(2), p1(1)-p0(1), p1(2)-p1(1), 0 , 'k'   );

    plot(p0(1), p0(2), '*');
    plot(p1(1), p1(2), '*');
end

%// calm dots
plot(cp1, f(r,cp1), 'c*') 
plot(cp2, f(r,cp2), 'c*')
    

%% Lyapunov
clc
clear

f = @(r,u) r.*u.*(4 - u.^2);
df = @(r,u) 4.*r - 3.*r.*(u.^2);
partL = @(r,u) log( abs(df(r,u)) );



% recur
n = 1000;
m = 1000;

r = linspace(0, 3*sqrt(3)/8 , m);
u_0 = 0.81;



u_res(1,:) = u_0*ones(1,m);
for j = 1:m     % every r
    for i = 2:n % every u
    
        u_res(i,j) = f(r(j), u_res(i-1,j));
    
    end
    
    h(j) = sum( partL( r(j), u_res(:,j))   ) / n;
    
end


plot(r,h,r,zeros(m,1));
axis([0 3*sqrt(3)/8 -7 1]);
%% Lyapunov test dana
clc
clear

f = @(r,u) r.*u.*(1 - sqrt(u) );
df = @(r,u) r - 3.*r.*sqrt(u)/2;
partL = @(r,u) log( abs(df(r,u)) );



% recur
n = 1000;
m = 1000;

r = linspace(0, 6.75 , m);
u_0 = 0.25;


u_res(1,:) = u_0*ones(1,m);

for j = 1:m     % every r
    for i = 2:n % every u
    
        u_res(i,j) = f(r(j), u_res(i-1,j));
    
    end
    
    h(j) = sum( partL( r(j), u_res(:,j))   ) / n;
    
end


plot(r,h, r,zeros(m,1));
%axis([0 2 0 1]);
%xlabel('r') 
%ylabel('Lyapunov Index') 