% S- step with RK3

function Snew = sstep(X,V,S,Mv,Tx,Tv,Mxe,dt)

    k1 = (X'*Tx*X)*S*(V'*Mv*V)' - (X'*Mxe*X)*S*(V'*Tv*V)';
    Stemp = S + 0.5*dt*k1;
    k2 = (X'*Tx*X)*Stemp*(V'*Mv*V)' - (X'*Mxe*X)*Stemp*(V'*Tv*V)';
    Stemp = S + 3/4*dt*k2;
    k3 = (X'*Tx*X)*Stemp*(V'*Mv*V)' - (X'*Mxe*X)*Stemp*(V'*Tv*V)';

    S = S + dt/9 * (2*k1 + 3*k2 + 4*k3);

    Snew = S;
end