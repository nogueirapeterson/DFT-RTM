program test

    use pyplot_module
    use mdle_io_utils
    use mdle_adds
    !use mdle_mdlem

    implicit none
    integer           :: nz, nx, nt, nb, sz, ns, s0x, gz, ds
    real              :: dz, dx, dt, fpeak, fmax
    character(len=64) :: velfile, dadofile, vtrue_file

    integer              :: i, is, nzb, nxb
    integer, allocatable :: sx(:)
    real, allocatable    :: vel(:,:), fonte(:), D(:,:), Im(:,:), Imt(:,:), Im_tmp(:,:)

    integer :: order
    type(pyplot) :: plt


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

        call modelagem (nx, nz, nt, dx, dz, dt, fonte, nb, nxb, nzb, sx(is), sz, gz, vel(nb+1:nz+nb,nb+1:nx+nb), D)

        call RTM_DFT(nx, nz, nt, dx, dz, dt, fmax, fonte, nb, nxb, nzb, sx(is), sz, gz, vel(nb+1:nz+nb,nb+1:nx+nb), D,Im)


        Imt = Imt + Im(nb+1:nb+nz,nb+1:nx+nb)

        write(20, rec=is) D

    end do
    call plt%initialize(legend=.true.)
    call plt%add_imshow(Imt)
    call plt%savefig('img.png', pyfile = 'img.py')
    !call plt%showfig()

    open (unit=20, file='imagem-ref.dir',  form='unformatted', access='direct', recl=nz*nx*4)
    read(20,rec=1) Im_tmp
    close(20)

    print*, "Residual", sum(abs(Im_tmp - Imt))/(nx * nz)

    open (unit=2, file='imagem.dir',  form='unformatted', access='direct', recl=nz*nx*4)
    write(2,rec=1) Imt
    close(2)

    deallocate(vel, sx, fonte, D, Im, Imt,Im_tmp)

    write(*,*)
    call timestamp( )

end program test
