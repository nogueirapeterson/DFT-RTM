program test_atenuation
    use mdle_adds
    integer,parameter:: nz=101, nx=101, nb=20, order=2
    real,parameter:: dx=10., dz=10.
    real:: A(nz+2*nb,nx+2*nb), L(nz+2*nb,nx+2*nb), coef(order+1)

    associate (nzb => nx+2*nb, nxb => nz+2*nb)
        A   = 1.0
        call calc_coef(order, coef)
        call op_l(order,coef, nzb,nxb,dz,dx,A,L)
    end associate

    print*, "sum", sum(L)
    if (sum(L) /= 0.) then
        error stop "failed"
    end if
end program test_atenuation
