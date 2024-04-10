%% Задача линейного управления
clc
clear
%\cdot x = A(t)x + Bu + f(t)
A = {@(x) 1, @(x) 2; @(x) 0.8, @(x) 1};
B = [0.4, 1; 0, 0.4];
f = {@(x) 0 , @(x) 0}';
t0 = 0; % начальное время

% P - множество допустимых управлений
r=2;
p1=0;
p2=1;



% Psi
n_psi = 40;
str_ang = 0;
fin_ang = 2*pi;
phi_t0 = [sin(str_ang:(fin_ang-str_ang)/n_psi:fin_ang-(fin_ang-str_ang)/n_psi); cos(str_ang:(fin_ang-str_ang)/n_psi:fin_ang-(fin_ang-str_ang)/n_psi)];
%Разбиение на варианты (сопряженной переменной)
A_sopr_sus = A';
% Матрица сопряженной системы

t_max=1;
% Максимальное время просчета
n_min=1;
T = t_max;
% T- оптимальное время
size_set = 500;
t_set = linspace(t0,t_max,size_set);
% стандартизированная сетка

ALL_PSI = zeros(2*n_psi,size_set);
ALL_U_OPT = zeros(2*n_psi,size_set);
ALL_X_OPT = zeros(2*n_psi,size_set);
% матрицы с результатами
%{
if (set_X1(x0)>0)
    T=t0;
    n_min = 1;
end
%}
if abs(B(1,1)*B(2,2)-B(2,1)*B(1,2))<0.01
    B(1,1) = B(1,1)+0.1;
    B(2,2) = B(2,2)+0.1;
end

for i=1:n_psi
    tspan_0 = [t0 t_max];
    opts1 = odeset('RelTol',1e-8,'AbsTol',1e-8,'Refine',6,'MaxStep',1e-2);
    [t_sp,psi_sp] = ode45(@(t_sp,psi_sp) odefcn(t_sp,psi_sp,A_sopr_sus), tspan_0, phi_t0(:,i),opts1);
    % Решение сопряженной системы
    psi = interp1(t_sp,psi_sp,t_set);
    ALL_PSI(2*i-1,:)=psi(:,1);
    ALL_PSI(2*i,:)=psi(:,2);
    % интерполяция решения на нормальную сетку
    u_opt = zeros(2,size_set);
    u_1 = zeros(1,size_set);
    u_2 = zeros(1,size_set);
    u_1(1:size_set) = B(1,1)*psi(1:size_set,1)+B(2,1)*psi(1:size_set,2);
    u_2(1:size_set) = B(1,2)*psi(1:size_set,1)+B(2,2)*psi(1:size_set,2);
    % u_1 u_2 вектор и P множество. Дальше ищем опорный вектор в напр u1 u2
    % множеству P в каждый момент времени
    for j=1:size_set        
        %[pho,supp_vec]=Ellipse_Lebesgue([u_1(j),u_2(j)]); вариант через
        %численное нахождение опорной функции
        u_opt(1,j) = r*u_1(j)/norm([u_1(j),u_2(j)])+p1;
        u_opt(2,j) = r*u_2(j)/norm([u_1(j),u_2(j)])+p2;
        % TODO раздуть
    end% Находим оптимальное управление в момент времени t_i
    ALL_U_OPT(2*i-1,:)=u_opt(1,:);
    ALL_U_OPT(2*i,:)=u_opt(2,:);
    %Заполняем матрицы управления информацией о полученом управлении
    
    tspan_1 = [t0 t_max];
    opts = odeset('Events',@StopEvents,'RelTol',1e-8,'AbsTol',1e-8,'Refine',6,'MaxStep',1e-2);
    [valss,ppoit]= FIRST_Lebesgue( [phi_t0(1,i),phi_t0(2,i)]);
    x0 =  ppoit;
    %x0=[0, 1]';
    [t_opt_,psi_opt_,te,ye,ie] = ode45(@(t_opt_,psi_opt_) odefcn_OC(t_opt_,psi_opt_,A,B,f,t_set,u_opt), tspan_1, x0,opts);
    % решаем диффур и получаем кандидата на оптимальную траекторию
    if length(te)==1
        if te<T
            T = te;
            n_min = i;
        end
    end
    % Если время минимальное то учитываем это
    
    psi_opt = interp1(t_opt_,psi_opt_,t_set);
    ALL_X_OPT(2*i-1,:)=psi_opt(:,1);
    ALL_X_OPT(2*i,:)=psi_opt(:,2);
    %Заполняем матрицы траектории информацией о полученной траектории
end% Решаем диффур %\cdot x = A(t)x + Bu + f(t) и получаем кандидатов на оптимальную траекторию
%Блок счета
%% Графики компонент оптимального управления
figure(1)
hold on
for i=1:n_psi
    plot(ALL_U_OPT(2*i-1,:),ALL_U_OPT(2*i,:),'g');
end
normaly_time = find(t_set<T);
u_last = zeros(2,length(normaly_time));
for i=1:length(normaly_time)
    u_last(1,i) = ALL_U_OPT(2*n_min-1,i);
    u_last(2,i) = ALL_U_OPT(2*n_min,i);
end
plot(u_last(1,1),u_last(2,1),'b.','MarkerSize',15);
plot(u_last(1,:),u_last(2,:),'b','LineWidth',2);
plot(ALL_U_OPT(2*n_min-1,:),ALL_U_OPT(2*n_min,:),'b','LineWidth',2);
title('Графики компонент оптимального управления u2(u1)')
%lgd = legend('u2(u1)','NumColumns',2);
xlabel('u opt_1');
ylabel('u opt_2');
axis equal;
%% Компонент оптимальной траектории
figure(2)
hold on
normaly_time = find(t_set<T);
ALL_X_OPT_last = ALL_X_OPT(:,1:1:length(normaly_time));
%ALL_X_OPT_last = ALL_X_OPT;
plot([ALL_X_OPT_last(2*n_min-1,length(normaly_time)),ALL_X_OPT_last(2*n_min-1,length(normaly_time))+1],[ALL_X_OPT_last(2*n_min,length(normaly_time)),ALL_X_OPT_last(2*n_min,length(normaly_time))],'k','LineWidth',2);
plot(ALL_X_OPT_last(2*n_min-1,:),ALL_X_OPT_last(2*n_min,:),'b','LineWidth',2);
for i=1:n_psi
    plot(ALL_X_OPT_last(2*i-1,:),ALL_X_OPT_last(2*i,:),'g');
end
drawSet(@FIRST_Lebesgue, 15);
drawSet(@SECOND_Lebesgue, 15);
plot(ALL_X_OPT_last(2*n_min-1,1),ALL_X_OPT_last(2*n_min,1),'r.','MarkerSize',30);
title('Графики компонент оптимального траектории x2(x1)')
xlabel('x1');
ylabel('x2');
axis equal;
hold off
%% Графики сопряженных переменных
figure(3)
hold on
normaly_time = find(t_set<T);
ALL_PSI_last = ALL_PSI(:,1:1:length(normaly_time));
%ALL_PSI_last = ALL_PSI;
plot(ALL_PSI_last(2*n_min-1,:),ALL_PSI_last(2*n_min,:),'b','LineWidth',2);
for i=1:n_psi
    plot(ALL_PSI_last(2*i-1,:),ALL_PSI_last(2*i,:),'g');
end
plot(ALL_PSI_last(2*n_min-1,1),ALL_PSI_last(2*n_min,1),'b.','MarkerSize',30);
title('Графики сопряженных переменных psi2(psi1)')

xlabel('psi1');
ylabel('psi2');
lgd = legend('psi^* opt','NumColumns',2);
axis equal;
hold off

%% Графики компонент оптимального управления от времени
figure(4)
normaly_time = find(t_set<T);
subplot(2,1,1);
hold on
ALL_U_OPT_last = ALL_U_OPT(:,1:1:length(normaly_time));
plot(t_set(1:1:length(normaly_time)),ALL_U_OPT_last(2*n_min-1,:),'b','LineWidth',2);
%ALL_U_OPT = ALL_PSI;
for i=1:n_psi
    plot(t_set(1:1:length(normaly_time)),ALL_U_OPT_last(2*i-1,:),'g');
end
title('График компоненты оптимального управления u1(t)')
xlabel('t');
ylabel('u1(t)');

subplot(2,1,2);
hold on

ALL_U_OPT_last = ALL_U_OPT(:,1:1:length(normaly_time));
plot(t_set(1:1:length(normaly_time)),ALL_U_OPT_last(2*n_min,:),'b','LineWidth',2);
%ALL_U_OPT = ALL_PSI;
for i=1:n_psi
    plot(t_set(1:1:length(normaly_time)),ALL_U_OPT_last(2*i,:),'g');
end
title('График компоненты оптимального управления u2(t)')
legend('u2(t)');
xlabel('t');
ylabel('u2(t)');
hold off


%% Графики сопряженных переменных от времени
figure(6)
normaly_time = find(t_set<T);
subplot(2,1,1);
hold on
ALL_PSI_last = ALL_PSI(:,1:1:length(normaly_time));
plot(t_set(1:1:length(normaly_time)),ALL_PSI_last(2*n_min-1,:),'b','LineWidth',2);
%ALL_U_OPT = ALL_PSI;
for i=1:n_psi
    plot(t_set(1:1:length(normaly_time)),ALL_PSI_last(2*i-1,:),'g');
end
title('График сопряженной переменной от времени psi1(t)')
xlabel('t');
ylabel('psi1(t)');
subplot(2,1,2);
hold on
ALL_PSI_last = ALL_PSI(:,1:1:length(normaly_time));
plot(t_set(1:1:length(normaly_time)),ALL_PSI_last(2*n_min,:),'b','LineWidth',2);
%ALL_U_OPT = ALL_PSI;
for i=1:n_psi
    plot(t_set(1:1:length(normaly_time)),ALL_PSI_last(2*i,:),'g');
end
title('График сопряженной переменной от времени psi2(t)')
legend('psi2(t)');
xlabel('t');
ylabel('psi2(t)');
hold off
%% Графики траекторий от времени
figure(7)
normaly_time = find(t_set<T);
subplot(2,1,1);
hold on
ALL_X_OPT_last = ALL_X_OPT(:,1:1:length(normaly_time));
plot(t_set(1:1:length(normaly_time)),ALL_X_OPT_last(2*n_min-1,:),'b','LineWidth',2);
%ALL_U_OPT = ALL_PSI;
for i=1:n_psi
    plot(t_set(1:1:length(normaly_time)),ALL_X_OPT_last(2*i-1,:),'g');
end
title('График траектории x1 от времени x_1(t)')
xlabel('t');
ylabel('x_1(t)');
subplot(2,1,2);
hold on
ALL_X_OPT_last = ALL_X_OPT(:,1:1:length(normaly_time));
plot(t_set(1:1:length(normaly_time)),ALL_X_OPT_last(2*n_min,:),'b','LineWidth',2);
%ALL_U_OPT = ALL_PSI;
for i=1:n_psi
    plot(t_set(1:1:length(normaly_time)),ALL_X_OPT_last(2*i,:),'g');
end
title('График траектории x2 от времени x_2(t)')
legend('x2(t)');
xlabel('t');
ylabel('x2(t)');
hold off
%% информация о выполнении программы
if T < t_max-0.001
    
    disp('Задача разрешима');
    disp('Оптимальное время');
    disp(T);
    disp('Номер угла');
    disp(n_min);
    disp('Вектор psi(t0) оптимальный');
    disp(phi_t0(:,n_min));
    if T~=t0
        normaly_time = find(t_set<T);
        sk_pr = -(ALL_PSI(2*n_min-1,length(normaly_time))*ALL_X_OPT(2*n_min-1,length(normaly_time))+ALL_PSI(2*n_min,length(normaly_time))*ALL_X_OPT(2*n_min,length(normaly_time)));
        op_func = -1.3;
        disp('Погреность из условия трансвенсальности для X1');
        disp(abs(sk_pr-op_func)/norm([ALL_PSI(2*n_min-1,length(normaly_time)),ALL_PSI(2*n_min,length(normaly_time))]));        
    end
else
    disp('Задача не разрешима');
    disp('Время просчета');
        disp(t_max);
end


%%
function dphidt = odefcn(t,phi,A)
    dphidt = zeros(2,1);
    dphidt(1) = A{1,1}(t)*phi(1)+A{1,2}(t)*phi(2);
    dphidt(2) = A{2,1}(t)*phi(1)+A{2,2}(t)*phi(2);
end
function [value,isterminal,direction] = StopEvents(t,y)
a=5;
g = @(x) max(max(-x(1)*x(2)+1,x(1)+x(2)-a), -x(1));

value = 0+(g([y(1), y(2)])>=0);     % Detect height = 0
isterminal = 1;   % Stop the integration
direction = 0;   % Negative direction only
end
function dphidt = odefcn_OC(t,phi,A,B,f,setka,u)
    u_1 = zeros(1,length(setka));
    u_2 = zeros(1,length(setka));
    u_1(1:length(setka)) = B(1,1)*u(1,1:length(setka))+B(1,2)*u(2,1:length(setka));
    u_2(1:length(setka)) = B(2,1)*u(1,1:length(setka))+B(2,2)*u(2,1:length(setka));

    U_1 = interp1(setka,u_1,t); % Interpolate the data set (ft,f) at time t
    U_2 = interp1(setka,u_2,t); % Interpolate the data set (ft,f) at time t

    dphidt = zeros(2,1);
    dphidt(1) = A{1,1}(t)*phi(1)+A{1,2}(t)*phi(2)+f{1}(t)+U_1;
    dphidt(2) = A{2,1}(t)*phi(1)+A{2,2}(t)*phi(2)+f{2}(t)+U_2;
end
%% ALL

function [val, point]  = SupportLebesgue_2(f,opts)
    l1=  opts.lx;
    l2 = opts.ly;
    A = [];
    b = [];
    fun = @(x) -(x(1)*l1+x(2)*l2);
    x0 = opts.x0;
    Aeq = [];
    beq = [];
    lb=[];
    ub=[];
    save params_2 f;
    nonlcon = @unitdisk_2;
    options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
    [x,fval] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
    val = -fval;
    point = x;
end

function [c,ceq] = unitdisk_2(x)
    load params_2 f;
    ceq = f(x);
    c = [];
end
function [val, point] = FIRST_Lebesgue(y)
    gamma = 0.5;
    alpha = 2;
    g = @(x) gamma-alpha*x(1)*x(1)-(x(2)-1).^4;
    s = struct('lx',y(1),'ly',y(2),'x0',[-2,-2]);
    [val, point] = SupportLebesgue_2(g,s);
end
function [val, point] = SECOND_Lebesgue(y)
    a = 5;
    g = @(x) max(max(-x(1)*x(2)+1,x(1)+x(2)-a), -x(1));

    s = struct('lx',y(1),'ly',y(2),'x0',[4,2]);
    [val, point] = SupportLebesgue_2(g,s);
end

function res = drawSet(rho,N)
    %t = linspace(0,2*pi,400);
    %x_r = cos(t)*2;
    %y_r = sin(t)-2;
    %plot(x_r,y_r);
    hold on
    p = linspace(0,2*pi-2*pi/N,N);
    x_t = cos(p);
    y_t = sin(p);
    [val, point] = rho([x_t(1),y_t(1)]);
    point1 = point;
    point_last = point;

    for i = 2:N
        [val, point] = rho([x_t(i),y_t(i)]);
        % строим внутр прямую
        alf = linspace(point_last(1),point(1),100);
        bet = linspace(point_last(2),point(2),100);
        plot(alf,bet,'r');
        % строим внешнюю прямую
        c1 = -(point_last(1)*x_t(i-1)+point_last(2)*y_t(i-1));
        c2 = -(point(1)*x_t(i)+point(2)*y_t(i));
        x0 = (c1*y_t(i)-c2*y_t(i-1))/(x_t(i)*y_t(i-1)-x_t(i-1)*y_t(i));
        y0 = (c2*x_t(i-1)-c1*x_t(i))/(x_t(i)*y_t(i-1)-x_t(i-1)*y_t(i));
        

        point_last = point;
    end
    % построение последней внутренней прямой
    alf = linspace(point_last(1),point1(1),100);
    bet = linspace(point_last(2),point1(2),100);
    plot(alf,bet,'r');
    % построение последней внешней прямой
    c1 = -(point_last(1)*x_t(N)+point_last(2)*y_t(N));
    c2 = -(point1(1)*x_t(1)+point1(2)*y_t(1));
    x0 = (c1*y_t(1)-c2*y_t(N))/(x_t(1)*y_t(N)-x_t(N)*y_t(1));
    y0 = (c2*x_t(N)-c1*x_t(1))/(x_t(1)*y_t(N)-x_t(N)*y_t(1));

end