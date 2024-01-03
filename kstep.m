% K- step with RK3

function Knew = kstep(K,V,Mx,Mv,Tx,Tv,Mxe,Cx,dt)
    
    delta = 0;

    k1 = -Tx*K*(V'*Mv*V)' + Mxe*K*(V'*Tv*V)' - delta*Cx*K;
    Ktemp = K + 0.5*dt*k1;
    k2 = -Tx*Ktemp*(V'*Mv*V)' + Mxe*Ktemp*(V'*Tv*V)' - delta*Cx*K;
%     Ktemp = K - dt*k1 + 2*dt*k2;
    Ktemp = K + 3/4*dt*k2;
    k3 = -Tx*Ktemp*(V'*Mv*V)' + Mxe*Ktemp*(V'*Tv*V)' - delta*Cx*K;

%     K = K + dt/6 * (k1 + 4*k2 + k3);
    K = K + dt/9 * (2*k1 + 3*k2 + 4*k3);
    
    Knew = K;
end