% Initial condition
t = 50;
dt = 0.05;
t_ = 0:dt:t;
xh = 4*pi/63;
vh = 12/255;
xDomain = 0:xh:4*pi;
vDomain = -6:vh:6;
Nx = length(xDomain);
Nv = length(vDomain);
alpha = 10^-2;
f = @(v,x) 1/sqrt(2*pi) .* exp(-(abs(v)).^2*0.5) .* (1.+alpha.*cos(0.5.*x));


% Grid set-up
[V,X] = meshgrid(vDomain,xDomain);
fixedInitial = f(V,X);


% define discretization matrices for spatial and velocity
% spatial (periodic boundary)
Mx = diag(2*xh/3*ones(1,Nx)) + diag(xh/6*ones(1,Nx-1),1) ...
    + diag(xh/6*ones(1,Nx-1),-1);
Mx_tmp = Mx;
Mx(1,:) = Mx_tmp(1,:) + Mx_tmp(end,:);
Mx(:,1) = Mx_tmp(:,1) + Mx_tmp(:,end);
Mx(end,:) = [];
Mx(:,end) = [];

Tx = zeros(Nx-1,Nx-1) + diag(0.5*ones(1,Nx-2),1) ...
    + diag(-0.5*ones(1,Nx-2),-1);
Tx(1,end) = -0.5;
Tx(end,1) = 0.5;

% velocity
Mv = diag(2*vh/3*ones(1,Nv)) + diag(vh/6*ones(1,Nv-1),1) ...
    + diag(vh/6*ones(1,Nv-1),-1);
% Mv_tmp = Mv;
% Mv(1,:) = Mv_tmp(1,:) + Mv_tmp(end,:);
% Mv(:,1) = Mv_tmp(:,1) + Mv_tmp(:,end);
% Mv(end,:) = [];
% Mv(:,end) = [];

Tv = zeros(Nv,Nv) + diag(0.5*ones(1,Nv-1),1) ...
    + diag(-0.5*ones(1,Nv-1),-1);
Tv(1,1) = -0.5;
Tv(end,end) = 0.5;

vTemp = [vDomain(:,1)-vh vDomain vDomain(:,end)+vh];
Mvv_diag = zeros(1,Nv);
Mvv_up = zeros(1,Nv-1);
Mvv_down = zeros(1,Nv-1);

for j = 2:Nv+1
    Mvv_diag(1,j-1) = vh/3 * (vTemp(:,j-1) + vTemp(:,j+1));
end
Mvv_diag(1,1) = vh/3 * vTemp(:,3) - vh^2/4;
Mvv_diag(end,end) = vh/3 * vTemp(:,end-2) + vh^2/4;

for j = 1:Nv-1
    Mvv_up(1,j) = -0.25*vh^2 + vh/3*(vDomain(:,j+1)-2*vDomain(:,j)) + ...
        0.5*(vDomain(:,j)*vDomain(:,j+1)-vDomain(:,j)^2);
end

for j = 2:Nv
    Mvv_down(1,j-1) = -0.25*vh^2 + vh/3*(vDomain(:,j)-2*vDomain(:,j-1)) + ...
        0.5*(vDomain(:,j)*vDomain(:,j-1)-vDomain(:,j-1)^2);
end
Mvv = diag(Mvv_diag) + diag(Mvv_down,-1) + diag(Mvv_up,1);
% Mvv_tmp = Mvv;
% Mvv(1,:) = Mvv_tmp(1,:) + Mvv_tmp(end,:);
% Mvv(:,1) = Mvv_tmp(:,1) + Mvv_tmp(:,end);
% Mvv(end,:) = [];
% Mvv(:,end) = [];


% CIP stabilizer (haven't applied them yet)
Cx = diag(-1*ones(1,Nx-1)) + diag(ones(1,Nx-2),1) + diag(ones(1,Nx-2),-1);
Cx(1,end) = 1;
Cx(end,1) = 1;

Cv = diag(-1*ones(1,Nv)) + diag(ones(1,Nv-1),1) + diag(ones(1,Nv-1),-1);
Cv(1,end) = 1;
Cv(end,1) = 1;



% figure;
% ax = gca;

r = [20];
for k = 1:length(r)
    [X0,S0,V0] = svd(fixedInitial);     % svd(fixedInitial) 일지, svd(fixedInitial(1:end-1,1:end)) 일지?
    X0 = X0(1:end-1,1:r(k));
    S0 = S0(1:r(k),1:r(k));
    V0 = V0(:,1:r(k));


    initVal = X0 * S0 * V0';
    figure;
    s = surf(V(1:end-1,1:end),X(1:end-1,1:end),initVal);

    elecEnergy = zeros(1,length(t_));
    totalMass = zeros(1,length(t_));
    totalEnergy = zeros(1,length(t_));
    totalEntropy = zeros(1,length(t_));

    for i = 1:length(t_)
        % Compute the density from X0, S0, V0
        fVal = X0 * S0 * V0';
        rho = 1 - vh*(0.5*(fVal(:,1) + fVal(:,end)) + sum(fVal(:,2:end-1),2));
        Eapprox = eField(Nx-1,xh,rho);

%         elecEnergy(1,i) = sum(abs(Eapprox).^2,"all")/xh;
        elecEnergy(1,i) = norm(Eapprox);

%         plot(ax,xDomain(1:end-1),Eapprox)
%         ylim([-0.02 0.02])
%         drawnow
       
        
%         Eapprox = electricField(massMat,fvec,xh,interPoints);
        % Compute electric field E and update Mxe
%         approxDiag = Eapprox' + [Eapprox(end-1,1) Eapprox(1:end-1,1)'];
%         approxUD = Eapprox(1:end-1,1)';
%         Mxe = diag(approxDiag.*(2*xh/3*ones(1,Nx))) + diag(approxUD.*(xh/6*ones(1,Nx-1)),1) ...
%             + diag(approxUD.*(xh/6*ones(1,Nx-1)),-1);
%         Eapprox = ones(Nx,1);
        Mxe = Eapprox' .* Mx;  % 다시 구해야하나?

        % K step with RK3
        K0 = X0 * S0;
        K1 = kstep(K0,V0,Mx,Mvv,Tx,Tv,Mxe,Cx,dt);
%         [Q1,R1] = qr(K1);
%         S_hat = R1(1:r(k),1:r(k));
%         X1 = K1/S_hat;
%         [X1,~,~] = svd(Mx);
        X1 = dummy(Mx,r(k));
%         X1 = X1(:,1:r(k));
        S_hat = X1\K1;

        % S step
        S_tilde = sstep(X1,V0,S_hat,Mv,Tx,Tv,Mxe,dt);
        
        % L step
        L0 = V0 * S_tilde';
        L1 = lstep(X1,L0,Mv,Mvv,Tx,Tv,Mxe,Cv,dt);
%         [Q2,R2] = qr(L1);
%         S1 = R2(1:r(k),1:r(k));
%         V1 = L1/S1;
%         S1 = S1';
        V1 = dummy(Mv,r(k));
%         [V1,~,~] = svd(Mv);
%         V1 = V1(:,1:r(k));
        S1 = V1\L1;
        S1 = S1';
        
        % Compute the next time step solution
        newVal = X1 * S1 * V1';
        

        S0 = S1;
        X0 = X1;
        V0 = V1;
         
        s.XData = V(1:end-1,1:end);
        s.YData = X(1:end-1,1:end);
        s.ZData = newVal;
         
       pause(0.01)

%         % the electric field energy
%         elecEnergy(1,i) = 0.5*xh*sum(Eapprox.^2,"all");
%         initialEE = 0.5*xh*sum(Einitial.^2,"all");
% 
%         % the particle number (checked)
        totalMass(1,i) = sum(newVal,"all");
        initialMass = sum(fixedInitial,"all");
% 
%         % the total energy
%         totalEnergy(1,i) = 0.5*elecEnergy(1,i) + sum(0.5*(vDomain.^2 .* newVal),"all");
%         initialEnergy = totalEnergy(1,1);
% 
%         % the total entropy (checked)
%         totalEntropy(1,i) = norm(newVal);
%         initialEntropy = norm(initVal);
% 
%         initVal = newVal;
  
        
    end
    
%     figure;
%       hold on
%       plot(t_, abs((totalMass-initialMass)/initialMass))
    plot(t_,elecEnergy)

     

% 
%     semilogy(t_, abs((totalEnergy-initialEnergy)/initialEnergy));
    
%       semilogy(t_,totalEnergy);
%       plot(t_,totalEntropy);
%     plot(t_,log(abs((totalEnergy-initialEnergy)/initialEnergy)));
%     legend('rank 5', 'rank 10', 'rank 20');
%     plot(t_,totalEntropy)
%     semilogy(t_,abs((totalEntropy-initialEntropy)/initialEntropy));
%     legend('rank 5', 'rank 10', 'rank 20');
    
end
% legend('rank 20', 'rank 40');
% title('change in entropy over time')
% xlabel('time')
% ylim([0 1e-2]);