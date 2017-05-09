PROGRAM SOR1

IMPLICIT NONE
REAL(KIND=8), PARAMETER ::epsilon1=1.0E-03
INTEGER, PARAMETER :: grid_points=11
INTEGER ,PARAMETER:: interval=grid_points-1
INTEGER  :: iterations=0
REAL(KIND=8) :: x2=1.0,x1=0.0,y2=1.0,y1=0.0
REAL(KIND=8) :: residuals,dx,dy
REAL(KIND=8), DIMENSION(0:interval,0:interval):: phi
REAL(KIND=8), DIMENSION(1:interval-1,1:interval-1) :: res_matrix

dx=(x2-x1)/REAL(interval)
dy=(y2-y1)/REAL(interval)


CALL sor(interval,dx,dy,epsilon1,iterations,residuals,phi,res_matrix)



END PROGRAM SOR1


SUBROUTINE sor(itrvl,hx,hy,eps,iter,res,u,r)

IMPLICIT NONE

REAL(KIND=8), INTENT(IN) :: eps,hx,hy
REAL(KIND=8), INTENT(INOUT) :: res
INTEGER, INTENT(IN):: itrvl
INTEGER, INTENT(INOUT):: iter
REAL(KIND=8), DIMENSION(0:itrvl,0:itrvl),INTENT(INOUT) :: u
REAL(KIND=8), DIMENSION(1:itrvl-1,1:itrvl-1), INTENT(INOUT):: r

INTEGER :: i,j,k ,ier
REAL(KIND=8), DIMENSION(0:itrvl,0:itrvl) :: u_new,u1
REAL(KIND=8), DIMENSION(0:itrvl) :: x,y,w
REAL(KIND=8) , PARAMETER:: w1=1.404

DO i=0,itrvl,1
x(i)=REAL(i)*hx
END DO

DO i=0,itrvl,1
y(i)=REAL(i)*hy
END DO

DO i=0,8
w(i)=(REAL(i)*0.25)
!write(*,*) w(i)
END DO


DO i=0,itrvl,1
  DO j=0,itrvl,1
    IF (i==0) THEN
    u(i,j)=-y(j)**2
    ELSE IF (j==0) THEN
    u(i,j)=x(i)**2
    ELSE IF (i==itrvl) THEN
    u(i,j)=1-y(j)**2
    ELSE IF (j==itrvl) THEN
    u(i,j)=x(i)**2-1
    ELSE
    u(i,j)=0.0
    END IF
  END DO
END DO

!u(0, 0:itrvl)=x(0:itrvl)**2
!u(0:itrvl, 0)=-y(0:itrvl)**2
!u(0:itrvl,itrvl)=1-y(0:itrvl)**2
!u(itrvl, 0:itvl)= x(0:itrvl)**2
!u(1:itrvl-1,1:itrvl-1):: 0.0

DO i=0,itrvl,1
  DO j=0,itrvl,1
    u_new(i,j)=u(i,j)
      u1(i,j)=u(i,j)
      END DO
END DO


OPEN(UNIT=101,FILE='initial_matrix.dat',STATUS='REPLACE',ACTION='WRITE',IOSTAT=ier)
DO i=0,itrvl,1
WRITE(101,*)(u(j,i),j=0,itrvl)
END DO
CLOSE(UNIT=101,STATUS='KEEP',IOSTAT=ier)



 
DO k=0,8,1
    u=u1
    u_new=u1
 DO iter=1,10
  DO i=1,itrvl-1,1
   DO j=1,itrvl-1,1
    u(i,j)=0.25*(u(i-1,j)+u(i+1,j)+u(i,j-1)+u(i,j+1))
    u(i,j)=((1-w(k))*u_new(i,j))+(w(k)*u(i,j))
    r(i,j)=(u(i,j)-u_new(i,j))**2
    u_new(i,j)=u(i,j)
   END DO
END DO



 OPEN(UNIT=111,FILE='final_matrix.dat',STATUS='REPLACE',ACCESS='APPEND',ACTION='WRITE',IOSTAT=ier)
 DO i=0,itrvl,1
 WRITE(111,*)(u(j,i),j=0,itrvl)
 END DO

  OPEN(UNIT=121,FILE='res_matrix.dat',STATUS='REPLACE',ACCESS='APPEND',ACTION='WRITE',IOSTAT=ier)
  DO i=1,itrvl-1,1
  WRITE(121,*)(r(j,i),j=1,itrvl-1)

  END DO

res=SQRT(SUM(r))

!IF(res < eps) EXIT
!WRITE(*,*) w(k),iter, res
    

  END DO
    OPEN(UNIT=131,FILE='res.dat',STATUS='REPLACE',ACCESS='APPEND',ACTION='WRITE',IOSTAT=ier)
    WRITE(131,*) w(k), iter, res
    
   WRITE(*,*) w(k),iter-1,res
END DO


CLOSE(UNIT=111,STATUS='KEEP',IOSTAT=ier)
CLOSE(UNIT=121,STATUS='KEEP',IOSTAT=ier)
CLOSE(UNIT=131,STATUS='KEEP',IOSTAT=ier)
END SUBROUTINE
