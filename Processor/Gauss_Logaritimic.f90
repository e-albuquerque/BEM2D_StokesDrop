SUBROUTINE Gauss_Logaritimic(npgs, eta, rho)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: npgs
    REAL,INTENT(OUT) :: eta(npgs)
    REAL,INTENT(OUT) :: rho(npgs)
    
    eta(1) = 0.90425944e-2; eta(2)=0.53971054e-1;
    eta(3)=0.13531134; eta(4)=0.24705169; 
    eta(5)=0.38021171; eta(6) = 0.52379159; 
    eta(7)=0.66577472; eta(8)=0.79419019; 
    eta(9)=0.89816102; eta(10)=0.96884798;
    
    rho(1)=0.12095474; rho(2)=0.18636310; 
    rho(3)=0.19566066; rho(4)=0.17357723; 
    rho(5)=0.13569597; rho(6)=0.93647084e-1; 
    rho(7) = 0.55787938e-1; rho(8)=0.27159893e-1; 
    rho(9)=0.95151992e-2; rho(10)=0.16381586e-2;
END SUBROUTINE