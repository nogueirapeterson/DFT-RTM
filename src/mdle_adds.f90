module mdle_adds

    contains

        subroutine RTM(is,nx, nz, nt, dx, dz, dt, fmax, fonte, nb, nxb, nzb, sx, sz, rsz, csuav, scg, Im)
            implicit none
            integer               :: nzb, nxb,isnap,is, i, j, ix, iz, it, sz, sx, nb, nt, nz, nx, rsz, dsx, order
            real, dimension(:,:) :: csuav(nz,nx), scg(nt,nx), Im(nzb,nxb)
            real,dimension (:)   :: fonte(nt)

            real, allocatable    :: p0(:,:), p1(:,:), p2(:,:),cext(:,:), L(:,:)
            real, allocatable    :: CF(:,:,:)
            real, allocatable    :: coef(:)
            real                 :: dt, dz, dx, fmax

            order = 8

            allocate(cext(nzb,nxb))
            allocate(p0(nzb,nxb), p1(nzb,nxb), p2(nzb,nxb),L(nzb,nxb))
            allocate(coef(order+1))
            allocate(CF(nz,nx,nt))

            call makevelextend(nz,nx,nzb,nxb,nb,csuav,cext)

            call calc_coef(order, coef)

            p0  =  0. ; p1  =  0. ; p2  =  0. ; L = 0. ; Im = 0.

            do it = 1, nt

                call atenuacao (nxb,nzb,nb,p0)
                call atenuacao (nxb,nzb,nb,p1)


                call op_l(order,coef, nzb,nxb,dz,dx,p1,L)

                p2 = 2*p1 - p0 + (dt**2)*(cext**2)*L

                p2(sz,sx) = p2(sz,sx) + fonte(it)

                CF(:,:,it) = p2(nb+1:nz+nb,nb+1:nx+nb)

                !scg(it,:) = p2(rsz,nb+1:nx+nb)

                p0 = p1
                p1 = p2

            end do

            p0  =  0. ; p1  =  0. ; p2  =  0. ; L = 0.


            do it = 1,nt

                call atenuacao (nxb,nzb,nb,p0)
                call atenuacao (nxb,nzb,nb,p1)

                call op_l(order,coef, nzb,nxb,dz,dx,p1,L)

                p2 = 2*p1 - p0 + (dt**2)*(cext**2)*L

                p2(rsz,nb+1:nx+nb) = p2(rsz,nb+1:nx+nb) + scg(nt-it+1,:)

                if (it == 100) then

                    open (unit=9, file='pr2.bin',  form='unformatted', access='direct', recl=nz*nx*4)
                    write(9,rec=1) p2(nb+1:nb+nz,nb+1:nx+nb)
                    close(9)

                end if

                p0 = p1
                p1 = p2

                Im(nb+1:nz+nb,nb+1:nx+nb) = Im(nb+1:nz+nb,nb+1:nx+nb) + CF(:,:,nt-it+1)*p2(nb+1:nz+nb,nb+1:nx+nb)

            end do

            deallocate (p0, p2, p1, L, cext, coef, CF)

        end subroutine RTM

        subroutine modelagem (is,nx, nz, nt, dx, dz, dt, fmax, fonte, nb, nxb, nzb, sx, sz, rsz, csuav, scg)
            implicit none
            integer ::  nzb, nxb,isnap,is, i, j, ix, iz, it, sz, sx, nb, nt, nz, nx, rsz, dsx, order
            real, dimension(:,:) :: csuav(nz,nx), scg(nt,nx)
            real,dimension (:)     :: fonte(nt)

            real, allocatable :: p0(:,:), p1(:,:), p2(:,:),cext(:,:), L(:,:)
            real, allocatable    ::   coef(:)
            real                 :: dt, dz, dx, fmax

            order = 8

            allocate(cext(nzb,nxb))
            allocate(p0(nzb,nxb), p1(nzb,nxb), p2(nzb,nxb),L(nzb,nxb))
            allocate(coef(order+1))

            call makevelextend(nz,nx,nzb,nxb,nb,csuav,cext)

            call calc_coef(order, coef)

            p0  =  0. ; p1  =  0. ; p2  =  0. ; L = 0.

            do it = 1, nt

                call atenuacao (nxb,nzb,nb,p0)
                call atenuacao (nxb,nzb,nb,p1)

                call op_l(order,coef, nzb,nxb,dz,dx,p1,L)

                p2 = 2*p1 - p0 + (dt**2)*(cext**2)*L

                p2(sz,sx) = p2(sz,sx) + fonte(it)


                scg(it,:) = p2(rsz,nb+1:nx+nb)

                p0 = p1
                p1 = p2

            end do

            deallocate (p0, p2, p1, L, cext, coef)

        end subroutine modelagem

        subroutine RTM_DFT(is,nx, nz, nt, dx, dz, dt, fmax, fonte, nb, nxb, nzb, sx, sz, rsz, csuav, scg, Im)
            implicit none
            integer ::  nzb, nxb,isnap,is, i, j, ix, iz, it, sz, sx, nb, nt, nz, nx, rsz, dsx, order, nw, iw, nw_i, nw_f
            real, dimension(:,:) :: csuav(nz,nx), scg(nt,nx), Im(nzb,nxb)
            real,dimension (:)     :: fonte(nt)

            real, allocatable :: p0(:,:), p1(:,:), p2(:,:),cext(:,:), L(:,:), kr(:,:), ki(:,:)
            real, allocatable :: pr0(:,:), pr1(:,:), pr2(:,:), Lr(:,:)
            !real, allocatable :: CF(:,:,:)
            real, allocatable    ::   coef(:), w(:)
            real                 :: dt, dz, dx, df, fmin, fmax

            real, allocatable :: psr(:,:,:), psi(:,:,:), prr(:,:,:), pri(:,:,:)

            df = 1./(nt*dt)

            order = 8

            fmin=0
            nw_i = int(fmin / df)
            nw_f = int(fmax / df) +1
            nw = (nw_f-nw_i) +1

            !nw = nt/2 + 1

            write(*,*) "nw", nw

            allocate(cext(nzb,nxb))
            allocate(p0(nzb,nxb), p1(nzb,nxb), p2(nzb,nxb),L(nzb,nxb))
            allocate(pr0(nzb,nxb),pr1(nzb,nxb),pr2(nzb,nxb),Lr(nzb,nxb))
            allocate(coef(order+1))
            !allocate(CF(nz,nx,nt))
            allocate(w(nw), kr(nw,nt), ki(nw,nt))
            allocate(psr(nzb,nxb, nw), psi(nzb,nxb,nw), prr(nzb,nxb, nw), pri(nzb,nxb, nw))

            call makevelextend(nz,nx,nzb,nxb,nb,csuav,cext)

            call calc_coef(order, coef)

            p0  =  0. ; p1  =  0. ; p2  =  0. ; L = 0. ; Im = 0.

            pr0  =  0. ; pr1  =  0. ; pr2  =  0. ; Lr = 0.

            psr = 0. ; psi=0. ; prr=0. ; pri = 0.

            do iw = 1, nw
                w(iw) = 2.*acos(-1.)*((iw-1)*df)*dt  ! Vetor de frequências
            end do

            do it = 1, nt
                do iw = 1, nw
                    kr(iw, it) = cos(w(iw)*(it-1)) ! Kernel parte real
                    ki(iw, it) = sin(w(iw)*(it-1)) ! Kernel parte imaginária
                end do
            end do

            do it = 1, nt

                call atenuacao (nxb,nzb,nb,p0)
                call atenuacao (nxb,nzb,nb,p1)

                call op_l(order,coef, nzb,nxb,dz,dx,p1,L)

                p2 = 2*p1 - p0 + (dt**2)*(cext**2)*L

                p2(sz,sx) = p2(sz,sx) + fonte(it)

                call atenuacao (nxb,nzb,nb,pr0)
                call atenuacao (nxb,nzb,nb,pr1)

                call op_l(order,coef, nzb,nxb,dz,dx,pr1,Lr)

                pr2 = 2*pr1 - pr0 + (dt**2)*(cext**2)*Lr

                pr2(rsz,nb+1:nx+nb) = pr2(rsz,nb+1:nx+nb) + scg(nt-it+1,:)

                !psr = 0. ; psi = 0. ; prr = 0. ; pri = 0.

                do iw = 1, nw

                    psr(:,:,iw) = psr(:,:,iw) + kr(iw, it)*p2  ! Aplica kernel parte real campo da fonte
                    psi(:,:,iw) = psi(:,:,iw) + ki(iw, it)*p2  ! Aplica kernel parte imaginária campo da fonte

                    prr(:,:,iw) = prr(:,:,iw) + kr(iw, it)*pr2 ! Aplica kernel parte real campo receptor
                    pri(:,:,iw) = pri(:,:,iw) + ki(iw, it)*pr2 ! Aplica kernel parte imaginária campo receptor

                    Im = Im +  w(iw)*w(iw)*(psr(:,:,iw)*prr(:,:,iw) - psi(:,:,iw)*pri(:,:,iw)) ! Condição de imagem

                    !Im = Im + (psr(:,:,iw)*prr(:,:,iw) - psi(:,:,iw)*pri(:,:,iw)) ! Condição de imagem

                end do

                p0 = p1
                p1 = p2

                pr0 = pr1
                pr1 = pr2

           end do

           deallocate (p0, p2, p1, L, cext, coef, w, kr, ki)
           deallocate(psr, psi, prr, pri)
           deallocate (pr0,pr1,pr2,Lr)

        end subroutine RTM_DFT

        subroutine atenuacao(nxb,nzb,nb,p2)

            real, dimension (:,:) :: p2(nzb,nxb)

            lz=nzb
            do i =1,nb
                do j=1,nxb
                    p2(i,j)=p2(i,j)*(exp(-1*(0.0005*(nb-i))))**2
                    p2(lz,j)=p2(lz,j)*(exp(-1*(0.0005*(nb-i))))**2
                enddo
                lz=lz-1
            enddo

            lx=nxb
            do i =1,nb
                do j=1,nzb
                    p2(j,i)=p2(j,i)*(exp(-1*(0.0005*(nb-i))))**2
                    p2(j,lx)=p2(j,lx)*(exp(-1*(0.0005*(nb-i))))**2
                enddo
                lx=lx-1
            enddo

            return
        end subroutine atenuacao

        subroutine makevelextend(nz,nx,nzb,nxb,nb,c,cext)

            implicit none
            real, dimension (:,:)     :: c(nz,nx), cext(nzb,nxb)
            integer         :: nz,nx,nb,nzb,nxb,iz,ix

            !Região Central

            cext(nb+1:nz+nb,nb+1:nx+nb)=c(:,:)

            !Lado Superior e Inferior
            do ix=1,nx
                do iz=1,nb
                    cext(iz,nb+ix)=c(1,ix)        !Lado Superior
                    cext(nz+iz+nb,nb+ix)=c(nz,ix)    !Lado Inferior
                enddo
            enddo

            !Lado Esquerdo e Direito
            do iz=1,nzb
                do ix=1,nb
                    cext(iz,ix)=cext(iz,nb+1)        !Lado Esquerdo
                    cext(iz,ix+nx+nb)=cext(iz,nb+nx)    !Lado Direito
                enddo
            enddo

            return
        end subroutine makevelextend


        !====================================================================================
        !                            Subrotina que faz a fonte
        !====================================================================================
        subroutine source(nt, dt, fpeak, fonte)

            implicit none
            integer            :: nt
            real               :: dt, fpeak
            real, dimension(:) :: fonte(nt)

            integer :: i
            real    :: tdelay, t, wpeak, waux, tt
            real    :: pi=3.141592653589793238462643383279502884197

            wpeak = 2.*pi*fpeak
            waux  = 0.5*wpeak
            tdelay = 6./(5.*fpeak)

            do i = 1,nt
                t = (i-1)*dt
                tt = t - tdelay

                fonte(i) = exp(-waux*waux*tt*tt/4.)*cos(wpeak*tt)
            end do

        end subroutine source

        !==================================================================
        !            Subrotina que gera o vetor do tapering
        !==================================================================
        subroutine get_sponge(ctap, lsp, frac)

            implicit none
            integer            :: lsp
            real               :: frac
            real, dimension(:) :: ctap(lsp)

            integer :: i
            real    :: dfrac

            if (frac .le. 0.0) then
                ctap = 0.0
                return
            endif

            dfrac = sqrt(-log(frac))/(1.*lsp)

            do i=1,lsp
                ctap(i) = exp(-((dfrac*(lsp-i+1.))**2))
            enddo

            return
        end subroutine get_sponge

        !==================================================================
        !            Subrotina que aplica o tapering na borda
        !==================================================================
        subroutine taper_apply(pp, nxb, nzb, nb, taper)

            implicit none
            integer              :: nxb, nzb, nb
            real, dimension(:,:) :: pp(nzb,nxb), taper(nb)

            integer :: ix, iz, jlx, jlz

            jlx = nxb
            do ix=1,nb
                do iz=1,nzb
                    pp(iz,ix) = pp(iz,ix)*taper(ix)
                    pp(iz,jlx) = pp(iz,jlx)*taper(ix)
                enddo
                jlx = jlx-1
            enddo

            do ix=1,nxb
                jlz = nzb
                do iz=1,nb
                    pp(iz,ix) = pp(iz,ix)*taper(iz)
                    pp(jlz,ix) = pp(jlz,ix)*taper(iz)
                    jlz = jlz-1
                enddo
            enddo

            return
        end subroutine taper_apply

        subroutine calc_coef(a, coef)

            integer :: a
            real, dimension (:) :: coef(a+1)

            !    n=2
            if(a.eq.(2))then
                coef(1) = 1.
                coef(2) = -2.
                coef(3) = 1.
            endif

            !    n=4
            if(a.eq.(4))then
                coef(1) = (-1.)/12.
                coef(2) = (4.)/3.
                coef(3) = (-5.)/2.
                coef(4) = (4.)/3.
                coef(5) = (-1.)/12.
            endif

            !    n=6
            if(a.eq.(6))then
                coef(1) = 1./90.
                coef(2) = (-3.)/20.
                coef(3) = 3./2.
                coef(4) = (-49.)/18.
                coef(5) = 3./2.
                coef(6) = (-3.)/20.
                coef(7) = 1./90.
            endif

            !    n=8
            if(a.eq.(8))then
                coef(1) = -1./560.
                coef(2) = 8./315.
                coef(3) = (-1.)/5.
                coef(4) = 8./5.
                coef(5) = (-205.)/72.
                coef(6) = 8./5.
                coef(7) = (-1.)/5.
                coef(8) = 8./315.
                coef(9) = -1./560.
            endif

        end subroutine calc_coef

        subroutine op_l(a,coef, nz,nx,dz,dx,P,lap_P)
            !Subrotina que calcula o laplaciano para as ordens 2, 4, 6 e 8.
            !nx=número de amostra no eixo x parâmetro de entrada
            !nz=número de amostra no eixo x parâmetro de entrada
            !lapla = Matriz de saída representando o laplaciano

            integer :: i, j, k, nx, nz, a, lim_nx, lim_nz, in_n
            real :: lap_P, dx, dz, Pxx, Pzz, P
            real, dimension(:)  :: coef(a+1)

            dimension P(nz,nx), lap_P(nz,nx)

            in_n = (a/2)+1

            lim_nx = nx-(a/2)
            lim_nz = nz-(a/2)

            do i=in_n,lim_nx

                do j=in_n,lim_nz

                    Pxx = 0.0
                    Pzz = 0.0

                    do k=1,a+1

                        Pxx = Pxx + coef(k)*P(j,i+k-in_n)
                        Pzz = Pzz + coef(k)*P(j+k-in_n,i)

                    enddo

                lap_p(j,i) = Pxx/(dx**2) + Pzz/(dz**2)

                enddo
            enddo

            return
        end subroutine op_l

end module mdle_adds
