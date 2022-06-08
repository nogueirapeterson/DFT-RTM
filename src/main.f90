
! Este é um código Fortran 90 da implementação da técnica RTM no domínio tempo-frequência utilizando um kernel de Fourier.
! Título : Migração reversa no tempo no domínio tempo-frequência utilizando kernel da transformada discreta de Fourier (DFT-RTM)
! Autor : Peterson Nogueira
! Data : 26/06/2017 (ultima atualizacao)

program main
    !use mpi
    !include 'mpif.h'

    use mdle_io_utils
    use mdle_adds
    !use mdle_mdlem

    implicit none
    integer           :: nz, nx, nt, nb, sz, ns, s0x, gz, fz, iz, ds
    real              :: dz, dx, dt, fpeak, fmax
    character(len=64) :: velfile, argfile, dadofile, vtrue_file

    integer              :: i, ix, j, is, it, nzb, nxb
    integer, allocatable :: sx(:)
    real, allocatable    :: vel(:,:), fonte(:), D(:,:), Im(:,:), Imt(:,:), Im_tmp(:,:)

    integer :: order


    !Parametros iniciais
    call inputdata(velfile, dadofile, vtrue_file, nz, nx, dz, dx, nt, dt, sz, ns, s0x, ds, gz, fpeak, nb)
    nzb = nz + 2*nb;    nxb = nx + 2*nb ; order = 8 ; fmax = 2.*fpeak

    ! Alocacao de memoria

    allocate(vel(nzb,nxb), Im(nzb,nxb), Imt(nz,nx), Im_tmp(nz,nx))

    allocate(fonte(nt), sx(ns), D(nt,nx))

    open (unit=9, file=velfile,  form='unformatted', access='direct', recl=nz*nx*4)
    read(9,rec=1) vel(nb+1:nb+nz,nb+1:nx+nb)
    close(9)

    ! Fontes e recetores

    sz = sz + nb
    gz = gz + nb

    do i = 1,ns
          sx(i) = nb + s0x + ds*(i-1)
    end do

    do is = 1, ns
        write(*,*) 'posicao fonte', is, sx(is)
    end do

    call source(nt, dt, fpeak, fonte)

    open (unit=20, file="dado.bin",  form='unformatted', access='direct', recl=nx*nt*4)
    D = 0. ; Imt = 0.

    do is = 1, ns

        write(*,*) 'Tiro: ',is
        write(*,*) 'Fonte: ', sx(is)-nb, sz-nb

        !subrotina que faz a RTM

        call modelagem (is,nx, nz, nt, dx, dz, dt, fmax, fonte, nb, nxb, nzb, sx(is), sz, gz, vel(nb+1:nz+nb,nb+1:nx+nb), D)

        call RTM_DFT(is,nx, nz, nt, dx, dz, dt, fmax, fonte, nb, nxb, nzb, sx(is), sz, gz, vel(nb+1:nz+nb,nb+1:nx+nb), D,Im)


        Imt = Imt + Im(nb+1:nb+nz,nb+1:nx+nb)

        write(20, rec=is) D

    end do

    open (unit=2, file='imagem.dir',  form='unformatted', access='direct', recl=nz*nx*4)
    write(2,rec=1) Imt
    close(2)

    deallocate(vel, sx, fonte, D, Im, Imt,Im_tmp)

    write(*,*)
    call timestamp( )

end program main
