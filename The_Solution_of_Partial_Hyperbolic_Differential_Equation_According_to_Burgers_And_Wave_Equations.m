clear all
clc
tic
%% The Solution of Partial Hyperbolic Diff.Eq. According to...
...Burger's and Wave Equations

%Constant Values
length = 70; %m
delta_t = [0.005 0.0025 0.00125 0.05];%s
delta_x = 1; %m
c = delta_t./delta_x; %Courant Number
n_x = length./delta_x +1;
n_t = 1.5./delta_t+ 1;
%Boundary Conditions
%% a)c=1
x_up = 2:2:20;
x_down = 18:-2:2;
u11(1:n_x,1:n_t(1)) = zeros;
u11(1:6,1) = 0;
u11(7:16,1) = x_up;
u11(17:25,1) = x_down;
u11(26:n_x,1) = 0;
n=1;
while n>=1 && n<=n_t(1)-1
    i=2;
    while i>=2 && i<=n_x
        u11(i,n+1) = u11(i,n) - c(1).*((u11(i,n).*u11(i,n))./2 - (u11(i-1,n).*u11(i-1,n))./2);
        i = i+1;
    end
    n = n+1;
end
%% a)c=0.5
x_up = 2:2:20;
x_down = 18:-2:2;
u12(1:n_x,1:n_t(2)) = zeros;
u12(1:6,1) = 0;
u12(7:16,1) = x_up;
u12(17:25,1) = x_down;
u12(26:n_x,1) = 0;

n=1;
while n>=1 && n<=n_t(2)-1
    i=2;
    while i>=2 && i<=n_x
        u12(i,n+1) = u12(i,n) - (c(2)./2).*(u12(i,n)^2 - u12(i-1,n)^2);
        i = i+1;
    end
    n = n+1;
end
%% a)c=0.25
x_up = 2:2:20;
x_down = 18:-2:2;
u13(1:n_x,1:n_t(3)) = zeros;
u13(1:6,1) = 0;
u13(7:16,1) = x_up;
u13(17:25,1) = x_down;
u13(26:n_x,1) = 0;

n=1;
while n>=1 && n<=n_t(3)-1
    i=2;
    while i>=2 && i<=n_x
        u13(i,n+1) = u13(i,n) - c(3).*((u13(i,n)^2)/2 - (u13(i-1,n)^2)/2);
        i = i+1;
    end
    n = n+1;
end

%% a)c=1 Updated for dt=0.05
x_up = 2:2:20;
x_down = 18:-2:2;
u14(1:n_x,1:n_t(4)) = zeros;
u14(1:6,1) = 0;
u14(7:16,1) = x_up;
u14(17:25,1) = x_down;
u14(26:n_x,1) = 0;

n=1;
while n>=1 && n<=n_t(4)-1
    i=2;
    while i>=2 && i<=n_x
        u14(i,n+1) = u14(i,n) - c(4).*((u14(i,n)^2)/2 - (u14(i-1,n)^2)/2);
        i = i+1;
    end
    n = n+1;
end
%% b)Lax-Wendroff c=1
x_up = 2:2:20;
x_down = 18:-2:2;
m11(1:n_x,1:n_t(1)) = zeros;
m11(1:6,1) = 0;
m11(7:16,1) = x_up;
m11(17:25,1) = x_down;
m11(26:n_x,1) = 0;

n=1;
while n>=1 && n<=n_t(1)-1
    i=2;
    while i>=2 && i<n_x
        m11(i,n+1) = m11(i,n) - (c(1)./2)*[(m11(i+1,n)^2)/2 - (m11(i-1,n)^2)/2] + ((c(1)^2)./4)*[(m11(i+1,n)+m11(i,n))*(((m11(i+1,n)^2)/2)-((m11(i,n)^2)/2)) - (m11(i,n)+m11(i-1,n))*(((m11(i,n)^2)/2)-((m11(i-1,n)^2)/2))]; 
        i = i+1;
    end
    n = n+1;
end
%% b) c=0.5
x_up = 2:2:20;
x_down = 18:-2:2;
m12(1:n_x,1:n_t(1)) = zeros;
m12(1:6,1) = 0;
m12(7:16,1) = x_up;
m12(17:25,1) = x_down;
m12(26:n_x,1) = 0;
n=1;
while n>=1 && n<=n_t(2)-1
    i=2;
    while i>=2 && i<n_x
        m12(i,n+1) = m12(i,n) - (c(2)./2)*[(m12(i+1,n)^2)/2 - (m12(i-1,n)^2)/2] + ((c(2)^2)./4)*[(m12(i+1,n)+m12(i,n))*(((m12(i+1,n)^2)/2)-((m12(i,n)^2)/2)) - (m12(i,n)+m12(i-1,n))*(((m12(i,n)^2)/2)-((m12(i-1,n)^2)/2))]; 
        i = i+1;
    end
    n = n+1;
end
%% b) c=0.25
x_up = 2:2:20;
x_down = 18:-2:2;
m13(1:n_x,1:n_t(1)) = zeros;
m13(1:6,1) = 0;
m13(7:16,1) = x_up;
m13(17:25,1) = x_down;
m13(26:n_x,1) = 0;
n=1;
while n>=1 && n<=n_t(3)-1
    i=2;
    while i>=2 && i<n_x
        m13(i,n+1) = m13(i,n) - (c(3)./2)*[(m13(i+1,n)^2)/2 - (m13(i-1,n)^2)/2] + ((c(3)^2)./4)*[(m13(i+1,n)+m13(i,n))*(((m13(i+1,n)^2)/2)-((m13(i,n)^2)/2)) - (m13(i,n)+m13(i-1,n))*(((m13(i,n)^2)/2)-((m13(i-1,n)^2)/2))]; 
        i = i+1;
    end
    n = n+1;
end

%% b) c=1 dt=0.05
x_up = 2:2:20;
x_down = 18:-2:2;
m14(1:n_x,1:n_t(4)) = zeros;
m14(1:6,1) = 0;
m14(7:16,1) = x_up;
m14(17:25,1) = x_down;
m14(26:n_x,1) = 0;
n=1;
while n>=1 && n<=n_t(4)-1
    i=2;
    while i>=2 && i<n_x
        m14(i,n+1) = m14(i,n) - (c(4)./2)*[(m14(i+1,n)^2)/2 - (m14(i-1,n)^2)/2] + ((c(4)^2)./4)*[(m14(i+1,n)+m14(i,n))*(((m14(i+1,n)^2)/2)-((m14(i,n)^2)/2)) - (m14(i,n)+m14(i-1,n))*(((m14(i,n)^2)/2)-((m14(i-1,n)^2)/2))]; 
        i = i+1;
    end
    n = n+1;
end
%% c) eps=0.1 c=1
x_up = 2:2:20;
x_down = 18:-2:2;
o11(1:n_x,1:n_t(4)) = zeros;
o11(1:6,1) = 0;
o11(7:16,1) = x_up;
o11(17:25,1) = x_down;
o11(26:n_x,1) = 0;
n=1;
while n>=1 && n<=n_t(4)-1
    i=2;
    while i>=2 && i<n_x
        o11(i,n+1) = o11(i,n) - 0.5*c(4)*(o11(i+1,n) - o11(i-1,n)) + (0.5*c(1)*c(1) + 0.1)*(o11(i+1,n) - 2*o11(i,n) + o11(i-1,n));
        i = i+1;
    end
    n = n+1;
end
%% c) eps=0.2 c=1
x_up = 2:2:20;
x_down = 18:-2:2;
o12(1:n_x,1:n_t(4)) = zeros;
o12(1:6,1) = 0;
o12(7:16,1) = x_up;
o12(17:25,1) = x_down;
o12(26:n_x,1) = 0;
n=1;
while n>=1 && n<=n_t(4)-1
    i=2;
    while i>=2 && i<n_x
        o12(i,n+1) = o12(i,n) - 0.5*c(4)*(o12(i+1,n) - o12(i-1,n)) + (0.5*c(1)*c(1) + 0.2)*(o12(i+1,n) - 2*o12(i,n) + o12(i-1,n));
        i = i+1;
    end
    n = n+1;
end
%% c) eps=0.3 c=1
x_up = 2:2:20;
x_down = 18:-2:2;
o13(1:n_x,1:n_t(4)) = zeros;
o13(1:6,1) = 0;
o13(7:16,1) = x_up;
o13(17:25,1) = x_down;
o13(26:n_x,1) = 0;
n=1;
while n>=1 && n<=n_t(4)-1
    i=2;
    while i>=2 && i<n_x
        o13(i,n+1) = o13(i,n) - 0.5*c(4)*(o13(i+1,n) - o13(i-1,n)) + (0.5*c(1)*c(1) + 0.3)*(o13(i+1,n) - 2*o13(i,n) + o13(i-1,n));
        i = i+1;
    end
    n = n+1;
end
%% d) c=1 
x_up = 2:2:20;
x_down = 18:-2:2;
d11(1:n_x,1:n_t(1)) = zeros;
d11(1:6,1) = 0;
d11(7:16,1) = x_up;
d11(17:25,1) = x_down;
d11(26:n_x,1) = 0;
n=1;
while n>=1 && n<=n_t(1)-1
    i=2;
    while i>=2 && i<=n_x
        d11(i,n+1) = d11(i,n) - c(1)*(d11(i,n) - d11(i-1,n)) - 1.5*c(1)*((d11(i,n)^2)/2 - (d11(i-1,n)^2)/2);
        i = i+1;
    end
    n = n+1;
end
%% d) c=0.5
x_up = 2:2:20;
x_down = 18:-2:2;
d12(1:n_x,1:n_t(2)) = zeros;
d12(1:6,1) = 0;
d12(7:16,1) = x_up;
d12(17:25,1) = x_down;
d12(26:n_x,1) = 0;
n=1;
while n>=1 && n<=n_t(2)-1
    i=2;
    while i>=2 && i<=n_x
        d12(i,n+1) = d12(i,n) - c(2)*(d12(i,n) - d12(i-1,n)) - 1.5*c(2)*((d12(i,n)^2)/2 - (d12(i-1,n)^2)/2);
        i = i+1;
    end
    n = n+1;
end
%% d) c=0.25
x_up = 2:2:20;
x_down = 18:-2:2;
d13(1:n_x,1:n_t(3)) = zeros;
d13(1:6,1) = 0;
d13(7:16,1) = x_up;
d13(17:25,1) = x_down;
d13(26:n_x,1) = 0;
n=1;
while n>=1 && n<=n_t(3)-1
    i=2;
    while i>=2 && i<=n_x
        d13(i,n+1) = d13(i,n) - c(3)*(d13(i,n) - d13(i-1,n)) - 1.5*c(3)*((d13(i,n)^2)/2 - (d13(i-1,n)^2)/2);
        i = i+1;
    end
    n = n+1;
end
%% d) c=1 dt=0.05
x_up = 2:2:20;
x_down = 18:-2:2;
d14(1:n_x,1:n_t(4)) = zeros;
d14(1:6,1) = 0;
d14(7:16,1) = x_up;
d14(17:25,1) = x_down;
d14(26:n_x,1) = 0;
n=1;
while n>=1 && n<=n_t(4)-1
    i=2;
    while i>=2 && i<=n_x
        d14(i,n+1) = d14(i,n) - c(4)*(d14(i,n) - d14(i-1,n)) - 1.5*c(4)*((d14(i,n)^2)/2 - (d14(i-1,n)^2)/2);
        i = i+1;
    end
    n = n+1;
end
%% Writing the values to Excel

excel_deltax={num2cell(0:1:70)'}; 
excel_deltat={num2cell(0:0.25:1.5)}; 
excel_u14=num2cell(u14(1:1:71,1:5:31)); 
ZZ14=['x/t',excel_deltat{1};excel_deltax{1},excel_u14];
xlswrite('FTBS(a)_Method_c=1_dt=0,05_guncelleme',ZZ14);

excel_deltax={num2cell(0:1:70)'}; 
excel_deltat={num2cell(0:0.25:1.5)}; 
excel_m14=num2cell(m14(1:1:71,1:5:31)); 
ZZm14=['x/t',excel_deltat{1};excel_deltax{1},excel_m14];
xlswrite('LaxWendroff_Method_c=1_dt=0,05_guncelleme',ZZm14);

excel_deltax={num2cell(0:1:70)'};
excel_deltat={num2cell(0:0.25:1.5)};
excel_o11=num2cell(o11(1:1:71,1:5:31)); 
ZZo11=['x/t',excel_deltat{1};excel_deltax{1},excel_o11];
xlswrite('IkinciMertebe_LaxWendroff_Method_c=1_eps=0,1_dt=0,05_updated',ZZo11);

excel_deltax={num2cell(0:1:70)'};
excel_deltat={num2cell(0:0.25:1.5)};
excel_o12=num2cell(o12(1:1:71,1:5:31)); 
ZZo12=['x/t',excel_deltat{1};excel_deltax{1},excel_o12];
xlswrite('IkinciMertebe_LaxWendroff_Method_c=1_eps=0,2_dt=0,05_updated',ZZo12);

excel_deltax={num2cell(0:1:70)'};
excel_deltat={num2cell(0:0.25:1.5)};
excel_o13=num2cell(o13(1:1:71,1:5:31)); 
ZZo13=['x/t',excel_deltat{1};excel_deltax{1},excel_o13];
xlswrite('IkinciMertebe_LaxWendroff_Method_c=1_eps=0,3_dt=0,05_updated',ZZo13);

excel_deltax={num2cell(0:1:70)'}; 
excel_deltat={num2cell(0:0.25:1.5)}; 
excel_d14=num2cell(d14(1:1:71,1:5:31)); 
ZZd14=['x/t',excel_deltat{1};excel_deltax{1},excel_d14];
xlswrite('Nonlinear_FTBS(d)_Method_c=1_dt=0,05_updated',ZZd14);
%% PLOTTING
figure('units','normalized','outerposition',[0 0 1 1]) 

subplot(1,2,1) %First Upwind Differencing Yöntemi c=1
[X11,Y11] = meshgrid(0:1:70, 0:0.05:1.5);
surface(X11',Y11',u14);
view(3)
grid on
title('Part (a) FTBS Method c=1')

subplot(1,2,2)
plot(0:1:70,u14(1:71,1),'-*')
hold on
plot(0:1:70,u14(1:71,6),'-o')
plot(0:1:70,u14(1:71,11),'-square')
plot(0:1:70,u14(1:71,16),'-x')
plot(0:1:70,u14(1:71,21),'-x')
plot(0:1:70,u14(1:71,26),'-x')
plot(0:1:70,u14(1:71,31),'-x')
title('Part (a) FTBS Method')

figure('units','normalized','outerposition',[0 0 1 1])

subplot(1,2,1) %Lax-Wendroff Method c=1
surface(X11',Y11',m14);
view(3)
grid on
title('Part (b) Lax-Wendroff Method c=1')

subplot(1,2,2)
plot(0:1:70,m14(1:71,1),'-*')
hold on
plot(0:1:70,m14(1:71,6),'-o')
plot(0:1:70,m14(1:71,11),'-square')
plot(0:1:70,m14(1:71,16),'-x')
plot(0:1:70,m14(1:71,21),'-x')
plot(0:1:70,m14(1:71,26),'-x')
plot(0:1:70,m14(1:71,31),'-x')
title('Part (b) Lax-Wendroff Method')

figure('units','normalized','outerposition',[0 0 1 1])

[X14,Y14] = meshgrid(0:1:70, 0:0.05:1.5);
subplot(2,2,1) %Lax-Wendroff 2.Order eps=0.1
surface(X14',Y14',o11);
view(3)
grid on
title('Part (c) Lax-Wendroff 2.Order eps=0.1')

subplot(2,2,2) %Lax-Wendroff 2.Order c=0.5
plot(0:1:70,o11(1:71,1),'-*')
hold on
plot(0:1:70,o11(1:71,6),'-o')
plot(0:1:70,o11(1:71,11),'-square')
plot(0:1:70,o11(1:71,16),'-x')
plot(0:1:70,o11(1:71,21),'-x')
plot(0:1:70,o11(1:71,26),'-x')
plot(0:1:70,o11(1:71,31),'-x')
title('Part (c) Lax-Wendroff 2.Order eps=0.3')

subplot(2,2,3) %Lax-Wendroff 2.Order c=0.25
plot(0:1:70,o12(1:71,1),'-*')
hold on
plot(0:1:70,o12(1:71,6),'-o')
plot(0:1:70,o12(1:71,11),'-square')
plot(0:1:70,o12(1:71,16),'-x')
plot(0:1:70,o12(1:71,21),'-x')
plot(0:1:70,o12(1:71,26),'-x')
plot(0:1:70,o12(1:71,31),'-x')
title('Part (c) Lax-Wendroff 2.Order eps=0.2')

subplot(2,2,4)
plot(0:1:70,o13(1:71,1),'-*')
hold on
plot(0:1:70,o13(1:71,6),'-o')
plot(0:1:70,o13(1:71,11),'-square')
plot(0:1:70,o13(1:71,16),'-x')
plot(0:1:70,o13(1:71,21),'-x')
plot(0:1:70,o13(1:71,26),'-x')
plot(0:1:70,o13(1:71,31),'-x')
title('Part (c) Lax-Wendroff 2.Order eps=0.3')

figure('units','normalized','outerposition',[0 0 1 1]) 

subplot(2,2,1) %FTBS c=1 (d)
[X11,Y11] = meshgrid(0:1:70, 0:0.05:1.5);
surface(X11',Y11',d14);
view(3)
grid on
title('Part (d) FTBS Method c=1')

subplot(2,2,4)
plot(0:1:70,d14(1:71,1),'-*')
hold on
plot(0:1:70,d14(1:71,6),'-o')
plot(0:1:70,d14(1:71,11),'-square')
plot(0:1:70,d14(1:71,16),'-x')
plot(0:1:70,d14(1:71,21),'-x')
plot(0:1:70,d14(1:71,26),'-x')
plot(0:1:70,d14(1:71,31),'-x')
title('Part (d) FTBS Method')
toc