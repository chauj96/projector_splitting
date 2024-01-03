% Solving poisson equation (for electric field) with 3rd order FDM +
% periodic boundary condition

function f = eField(Nx,xh,rho)
    A = -11 * eye(Nx);
    AU1 = diag(18 * ones(1,Nx-1),1);
    AU2 = diag(-9 * ones(1,Nx-2),2);
    AU3 = diag(2 * ones(1,Nx-3),3);

    A = A + AU1 + AU2 + AU3;
    A(end-2,1) = 2;
    A(end-1,1) = -9;
    A(end-1,2) = 2;
    A(end,:) = 1;
    
    rho(end,1) = 0;
    f = -6 * xh * A\rho;

end