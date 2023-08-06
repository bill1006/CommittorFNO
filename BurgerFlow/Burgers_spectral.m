function storage = Burgers_spectral(q)
% solves u_t + (0.5u^2)_x = nu*uxx, i.e.,
global nu
nu = .1;

%binary option
init_data = 2;

orszac = 1;

N = 1600;

storage = zeros(N, 6);


% x ticks 
x = linspace(-pi,pi,N+1);

x(N + 1) = [];
k = [-N/2 : N/2 - 1];
u = zeros(1,N+1);

% initial data
if init_data == 1; 
    u0 = exp(-4*x.^2);
end
if init_data == 2;
    u0 = max(0,sign(1 -abs(x)));

end

if init_data == 3;
    u0 = x
    u0(1:q) = 0
    u0(q+1: N) = 1
end 

if init_data == 4;
    u0 = normrnd(0,1, [1, N])
end

if init_data == 5;
    u0 = 1./(1+exp(-10*(x-x(q))))
    
end

if init_data == 6;
    u0 = 1./(1+exp(10*(x-x(q))))
    
end

storage(:, 1) = real(u0)

dt = 0.01; % time step


%draw 

% figure(1); clf; 
% hpic = plot(x,u0,'LineWidth',2,'color','r'); % plot of the numerical solution
% axis([-pi pi -0.01 1.01]);
% hold on;
% grid
% drawnow


tmax = 1;
t = 0;
freq = k; % frequencies
freq2 = freq.^2;
e2=exp(-nu*freq2*dt);
% figure(2); clf; grid;
vk=fftshift(fft(u0)); % v in the Fourier space
% sp = plot(k,abs(vk));
% set(gca,'Yscale','log');


while (t<tmax) 
    t=t+dt;
    vk=fftshift(fft(u0)); % v in the Fourier space
    % if orszac == 1
    %     ind = find(abs(k) > N/3);
    %     vk(ind) = 0;
    % end
    k1=rhs(0,vk);
    k2=rhs(0.5*dt,vk+0.5*dt*k1);
    k3=rhs(0.5*dt,vk+0.5*dt*k2);
    k4=rhs(dt,vk+dt*k3);
    vkn=vk+dt*(k1+2*k2+2*k3+k4)/6;
    
    un=ifft(ifftshift(e2.*vkn)); % return to u in the x-space
    
    if t > tmax; 
        d = (t - tmax)/dt;
        un = u0*d + un*(1 - d);
        fprintf('t = %d, tmax = %d, umax(tmax) = %d\n',t,tmax,max(un));
    end    

    %draw 
    % set(hpic,'xdata',x,'ydata',real(un));

    u0=un;

    % drawnow


   % if mod(t,tmax/10) <= dt
   %  figure(2); hold on;
   %  set(sp,'Ydata',abs(e2.*vkn));
   % end
    if min(isfinite(un)) == 0
        break;
    end

end

storage(:,2) = real(u0)


end
%%







function rk=rhs(dt,v);
global nu
[m N]=size(v);
k=[-N/2 : N/2 - 1];
freq =k;
freq2 = freq.^2;
e2=exp(-nu*freq2*dt);
em2=exp(nu*freq2*dt);
    vk1=v.*e2;          % e^{tL}v in the Fourier space 
    v2=ifft(ifftshift(vk1));      % exp(tL)v in the x-space
    v2=0.5*v2.^2;          % [exp(tL)v]^2 in the x-space
    rk=em2.*(-1i*freq).*fftshift(fft(v2)); % exp(-tL)[[(exp(tL)v)]_x] in the Fourier space
end
