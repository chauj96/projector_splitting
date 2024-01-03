% L- step with RK3

function Lnew = lstep(X,L,Mv,Mvv,Tx,Tv,Mxe,Cv,dt)
    delta = 0;

    k1 = -Mvv*L*(X'*Tx*X)' + Tv*L*(X'*Mxe*X)' - delta*Cv*L;
    Ltemp = L + 0.5*dt*k1;
    k2 = -Mvv*Ltemp*(X'*Tx*X)' + Tv*Ltemp*(X'*Mxe*X)' - delta*Cv*L;
    Ltemp = L + 3/4*dt*k2;
    k3 = -Mvv*Ltemp*(X'*Tx*X)' + Tv*Ltemp*(X'*Mxe*X)' - delta*Cv*L;

    L = L + dt/9 * (2*k1 + 3*k2 + 4*k3);


    Lnew = L;
end