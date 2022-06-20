module mdle_io_utils

    contains

        !Contém:
        !	subroutine inputdata
        !	subroutine report
        !	subroutine timestamp

        !==================================================================================
        !	Subrotina que faz a leitura dos parametros iniciais do programa
        !==================================================================================
        subroutine inputdata(velfile,dadofile,vtrue_file, nz, nx, dz, dx, nt, dt, sjz, nsx, s0x, ds, rjz, peakf, abcWidth)

            implicit none
            integer           :: nz, nx, nt, sjz, nsx, s0x, rjz, abcWidth, ds
            real              :: dz, dx, dt, peakf
            character(len=64) :: velfile, dadofile, vtrue_file

            character(len=64) :: argfile

            if(1 /= command_argument_count()) then
                print*,'wrong number of arguments'
                stop
            endif
            call getarg(1, argfile)

            open(30,file=argfile,status='unknown',action='read',form='formatted')

            read(30,'(t15,a64)')   velfile
            read(30,'(t15,a64)')   dadofile
            read(30,'(t15,a64)')   vtrue_file
            read(30,'(t15,i10)')   nz
            read(30,'(t15,i10)')   nx
            read(30,'(t15,f10.4)') dz
            read(30,'(t15,f10.4)') dx
            read(30,'(t15,i10)')   nt
            read(30,'(t15,f10.4)') dt
            read(30,'(t15,i10)')   sjz
            read(30,'(t15,i10)')   nsx
            read(30,'(t15,i10)')   s0x
            read(30,'(t15,i10)')   ds
            read(30,'(t15,i10)')   rjz
            read(30,'(t15,f10.4)') peakf
            read(30,'(t15,i10)')   abcWidth

            close(30)

            return
        end subroutine inputdata

        !==================================================================================
        !			Subrotina que imprime a data e hora
        !==================================================================================
        subroutine timestamp( )

            integer,dimension(8) :: values

            call date_and_time(values=values)
            write(*,'(1x,a5,i5,a1,i2,a1,i2,8x,a5,i3,a1,i2,a1,i2)') 'Data:',values(1),'/',values(2),'/',values(3), &
                                                                   'Hora:',values(5),':',values(6),':', values(7)
            write(*,*)

            return
        end subroutine timestamp

end module mdle_io_utils
