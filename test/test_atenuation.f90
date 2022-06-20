program test_atenuation
    use mdle_adds
    use pyplot_module
    integer,parameter:: nz=101, nx=101, nb=20
    real:: A(nz+2*nb,nx+2*nb)
    type(pyplot) :: plt

    A = 1.0


    associate (nzb => nx+2*nb, nxb => nz+2*nb)
        call atenuacao(nxb,nzb,nb,A)
    end associate
    call plt%initialize(legend=.true.)
    call plt%add_imshow(A)
    !call plt%savefig('img.png', pyfile = 'img.py')
    call plt%showfig()
end program test_atenuation

