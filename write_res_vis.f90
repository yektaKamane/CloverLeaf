

SUBROUTINE write_my_energy(current_timestep)

    USE clover_module
    USE timestep_module
    USE viscosity_module
    USE PdV_module
    USE accelerate_module
    USE flux_calc_module
    USE advection_module
    USE reset_field_module
    USE mpi_interface

    integer(c_int), intent(in) :: current_timestep
    integer(c_int) :: rank, err, unit_number
    integer :: i, j, nx, ny
    character(len=30) :: filename
    character(len=*), parameter :: output_dir = './output/'

    call my_MPI_Comm_rank(MPI_COMM_WORLD, rank, err)

    ! Set the file name to include the rank and unit number
    write(filename, '(A, "array_data_", I0, "_", I0, ".txt")') &
                        trim(output_dir), rank, current_timestep
    unit_number = 10 + rank

    ! Open the file for writing
    open(unit=unit_number, file=filename, status='unknown', action='write')

    nx = 480
    ny = 480
    ! Write the array elements to the file
    do i = 1, nx
        do j = 1, ny
            write(unit_number, '(F6.2)', advance='no') chunk%tiles(1)%field%energy0(i, j)
            ! print *, chunk%tiles(1)%field%energy0(i, j)
            if (j < ny) then
                write(unit_number, '(A)', advance='no') ', '
            else
                write(unit_number, *) ! Move to the next line
            end if
        end do
    end do

    close(unit_number)


END SUBROUTINE


SUBROUTINE save_initial_state(buffer)

    REAL(KIND=8), POINTER :: buffer(:)

END SUBROUTINE