SUBROUTINE mat_prod_2x2_2x6(A,B,prod)
    IMPLICIT NONE
    REAL,INTENT(IN) :: A(2,2)
    REAL,INTENT(IN) :: B(2,6)
    
    REAL,INTENT(OUT) :: prod(2,6)
    
    INTEGER :: i, j
    
    DO i=1,2
        DO j=1,6
            prod(i,j) = A(i,1)*B(1,j) + A(i,2)*B(2,j)
        END DO
    END DO
END SUBROUTINE