program test
    !use mpi
    !include 'mpif.h'
    use,intrinsic :: iso_fortran_env, only: wp => real32
    use pyplot_module
    use mdle_io_utils
    use mdle_adds
    !use mdle_mdlem

    implicit none
    integer, parameter :: nz=51, nx=51, nt=1500, nb=10, sz=25, ns=51, s0x=51, gz=1, ds=1
    real,parameter     :: dz=25, dx=8, dt=0.0008, fpeak=6.5
    character(len=64)  :: velfile, dadofile, vtrue_file

    integer               :: i, is, nzb, nxb
    integer, allocatable  :: sx(:)
    real(wp), allocatable :: vel(:,:), fonte(:), D(:,:), Im(:,:), Imt(:,:), Im_tmp(:,:)
    real                  :: fmax

    integer :: order

    type(pyplot) :: plt

    allocate(vel(51,51))

    vel(1:25,:) = 1
    vel(26:51,:) = 2

    call plt%initialize(legend=.true.)
    call plt%add_imshow(vel)
    call plt%showfig()

    read (*,*)

    !!Parametros iniciais
    !call inputdata(velfile, dadofile, vtrue_file, nz, nx, dz, dx, nt, dt, sz, ns, s0x, ds, gz, fpeak, nb)
    !nzb = nz + 2*nb;    nxb = nx + 2*nb ; order = 8 ; fmax = 2.*fpeak

    !! Alocacao de memoria

    !allocate(vel(nzb,nxb), Im(nzb,nxb), Imt(nz,nx), Im_tmp(nz,nx))

    !allocate(fonte(nt), sx(ns), D(nt,nx))

    !open (unit=9, file=velfile,  form='unformatted', access='direct', recl=nz*nx*4)
    !read(9,rec=1) vel(nb+1:nb+nz,nb+1:nx+nb)
    !close(9)

    !! Fontes e recetores

    !sz = sz + nb
    !gz = gz + nb

    !do i = 1,ns
          !sx(i) = nb + s0x + ds*(i-1)
    !end do

    !do is = 1, ns
        !write(*,*) 'posicao fonte', is, sx(is)
    !end do

    !call source(nt, dt, fpeak, fonte)

    !open (unit=20, file="dado.bin",  form='unformatted', access='direct', recl=nx*nt*4)
    !D = 0. ; Imt = 0.

    !do is = 1, ns

        !write(*,*) 'Tiro: ',is
        !write(*,*) 'Fonte: ', sx(is)-nb, sz-nb

        !!subrotina que faz a RTM

        !call modelagem (nx, nz, nt, dx, dz, dt, fonte, nb, nxb, nzb, sx(is), sz, gz, vel(nb+1:nz+nb,nb+1:nx+nb), D)

        !call RTM_DFT(nx, nz, nt, dx, dz, dt, fmax, fonte, nb, nxb, nzb, sx(is), sz, gz, vel(nb+1:nz+nb,nb+1:nx+nb), D,Im)


        !Imt = Imt + Im(nb+1:nb+nz,nb+1:nx+nb)

        !write(20, rec=is) D

    !end do

    !open (unit=2, file='imagem.dir',  form='unformatted', access='direct', recl=nz*nx*4)
    !write(2,rec=1) Imt
    !close(2)

    !deallocate(vel, sx, fonte, D, Im, Imt,Im_tmp)

    !write(*,*)
    !call timestamp( )

end program test
