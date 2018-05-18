program elliptic
implicit none
integer::i,j,nitr,nx,ny,nom,k
real*8 ::dx,dy,x0,xL,y0,yL,omega,tol
real*8::dom,lmax,pi,optom,om1,om2
real*8,allocatable ::u(:,:),x(:),y(:),e(:,:)

!Domain
x0 = 0.0d0 !left
xL = 2.0d0 !right

y0 = 0.0d0 !bottom
yL = 1.0d0 !up

!number of points
nx = 20
ny = 10

!grid spacing (spatial)
dx = (xL-x0)/dfloat(nx)
dy = (yL-y0)/dfloat(ny)

!spatial coordinates 
allocate(x(0:nx))
do i=0,nx
x(i) = x0 + dfloat(i)*dx
end do

allocate(y(0:ny))
do j=0,ny
y(j) = y0 + dfloat(j)*dy
end do

allocate(e(0:nx,0:ny))
pi = 4.0d0*datan(1.0d0)

!Iterative schemes to solve Laplace Equation:
allocate(u(0:nx,0:ny))

open(19,file='error.plt')
write(19,*) 'variables ="n","Error"'

open(20,file='nitr.plt')
write(20,*) 'variables ="w","N-ITR"'

!Tolerance
tol = 1.0d-10


!Compute with various omega and plot number of iterations required by the specified omega
om1=1.0d0
om2=1.95d0
dom = 0.05d0
nom = 19

do k=0,nom

  	!initial guess which satisfies boundary condtions
	!Homogeneous Drichlet B.C.
	!we are not updating boundary points, they are zero
	do j=0,ny
	do i=0,nx
	u(i,j) = 0.0d0
	end do
	end do

  	omega = om1 + dfloat(k)*dom
       
	call SOR(nx,ny,y,u,omega,tol,nitr) !Successive over relaxation (SOR)

	!write error
    write(*,*)'omega and total interval:'
    write(20,*) omega, nitr
    write(*,*) omega, nitr    
end do

!call GS(nx,ny,y,u,tol,nitr)    !Gauss-Seidel (GS)

close(20)
close(19)

lmax=0.5d0*(dcos(pi/dfloat(nx))+dcos(pi/dfloat(ny)))
optom=2.0d0/(1.0d0 + dsqrt(1-lmax**2))


print*,"-------------------------"
print*,"optimum omega = ", optom
print*,"-------------------------"

open(200,file='field.plt')
write(200,*)'variables ="x","y","u"'
write(200,*)'zone f=point i=',nx+1,',j=',ny+1
do j=0,ny
do i=0,nx
write(200,*) x(i),y(j),u(i,j)
end do
end do
close(200)

end


!-----------------------------------------------------------------------------!
!Gauss-Seidel (GS)
!-----------------------------------------------------------------------------!
subroutine GS(nx,ny,y,u,tol,nitr)  
implicit none
integer::nx,ny,i,j,nitr
real*8 ::tol,err,omega
real*8 ::u(0:nx,0:ny),e(0:nx,0:ny),v(0:nx,0:ny),y(0:ny)


omega = 1.0d0

nitr = 0
err = 1.0d0
    
do while(err.gt.tol)

	nitr = nitr + 1

    !dirichlet boundary condition
 do j=0,ny
    u(0,j) = y(j)
 end do

 do i = 1,9
   u(i,ny) = 1.0d0
 end do

    do j=0,ny
	do i=0,nx
	v(i,j) = u(i,j)
	end do
	end do
  
	!update
    do j=1,ny-1
	do i=1,nx-1
	u(i,j) = v(i,j) + ((omega/4.0d0)*(v(i+1,j)+u(i-1,j)+v(i,j+1)+u(i,j-1)-(4.0d0*v(i,j))))
	end do
	end do
    
    !neumann boundary condition
    do j = 1,5
      u(nx,j) = (1.0d0/11.0d0)*(18.0d0*u(nx-1,j) - 9.0d0*u(nx-2,j) + 2.0d0*u(nx-3,j))
    end do

   do j = 5,ny
   do i = 10,nx
     u(i,j) = 1.0d0
   end do
 end do
    
    !compute error
 	do j=0,ny
	do i=0,nx
	e(i,j) = dabs(u(i,j)-v(i,j))
	end do
	end do   
    !max value of the error
        err = maxval(e)

	!write error
    write(19,*) nitr, err
    write(*,*) nitr, err
        
end do

end 

!-----------------------------------------------------------------------------!
!Successive over relaxation (SOR)
!-----------------------------------------------------------------------------!
subroutine SOR(nx,ny,y,u,omega,tol,nitr)  
implicit none
integer::nx,ny,i,j,nitr
real*8 ::tol,err,omega
real*8 ::u(0:nx,0:ny),e(0:nx,0:ny),v(0:nx,0:ny),y(0:ny)

nitr = 0
err=1.0d0
    
do while(err.gt.tol)

	nitr = nitr + 1

!dirichlet boundary condition
 do j=0,ny
    u(0,j) = y(j)
 end do

 do i = 1,9
   u(i,ny) = 1.0d0
 end do

    do j=0,ny
	do i=0,nx
	v(i,j) = u(i,j)
	end do
	end do

	!update
    do j=1,ny-1
	do i=1,nx-1
	u(i,j) = v(i,j) + (omega/4.0d0)*(v(i+1,j)+u(i-1,j)+v(i,j+1)+u(i,j-1)-(4.0d0*v(i,j)))
	end do
	end do

    ! Neumann boundary condition
    do j = 1,5
      u(nx,j) = (1.0d0/11.0d0)*(18.0d0*u(nx-1,j) - 9.0d0*u(nx-2,j) + 2.0d0*u(nx-3,j))
    end do

   do j = 5,ny
   do i = 10,nx
     u(i,j) = 1.0d0
   end do
 end do

    !compute error
 	do j=0,ny
	do i=0,nx
	e(i,j) = dabs(u(i,j)-v(i,j))
	end do
	end do   
    !max value of the error
    err = maxval(e)

end do
end 
