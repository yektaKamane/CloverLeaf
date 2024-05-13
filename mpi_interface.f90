module mpi_interface
    use iso_c_binding
    use mpi

    interface
        subroutine raise_sigint_c() bind(C, name="raise_sigint")
            use iso_c_binding
        end subroutine raise_sigint_c

        subroutine my_MPI_Init(Fcomm_w, Fcomm_s, ierr) bind(C, name="my_MPI_Init")
            use iso_c_binding
            integer(c_int), value :: Fcomm_w
            integer(c_int), value :: Fcomm_s
            integer(c_int) :: ierr
        end subroutine my_MPI_Init

        subroutine my_MPI_Comm_rank(Fcomm, rank, ierr) bind(C, name="my_MPI_Comm_rank")
            use iso_c_binding
            integer(c_int), value :: Fcomm
            integer(c_int), intent(out) :: rank
            integer(c_int), intent(out) :: ierr
        end subroutine my_MPI_Comm_rank

        subroutine my_MPI_Comm_size(Fcomm, size, ierr) bind(C, name="my_MPI_Comm_size")
            use iso_c_binding
            integer(c_int), value :: Fcomm
            integer(c_int), intent(out) :: size
            integer(c_int), intent(out) :: ierr
        end subroutine my_MPI_Comm_size

        subroutine my_MPI_Barrier(Fcomm, ierr) bind(C, name="my_MPI_Barrier")
            use iso_c_binding
            integer(c_int), value :: Fcomm
            integer(c_int), intent(out) :: ierr
        end subroutine my_MPI_Barrier

        subroutine my_MPI_Abort(Fcomm, errorcode, ierr) bind(C, name="my_MPI_Abort")
            use iso_c_binding
            integer(c_int), value :: Fcomm
            integer(c_int), value :: errorcode
            integer(c_int), intent(out) :: ierr
        end subroutine my_MPI_Abort

        subroutine my_MPI_Waitall(count, array_of_requests, array_of_statuses, ierr) bind(C, name="my_MPI_Waitall")
            use iso_c_binding
            integer(c_int), value :: count
            integer(c_int), dimension(*), intent(in) :: array_of_requests
            integer(c_int), dimension(*), intent(in) :: array_of_statuses
            integer(c_int), intent(out) :: ierr
        end subroutine my_MPI_Waitall

        subroutine my_MPI_Wait(request, status, ierr) bind(C, name="my_MPI_Wait")
            use iso_c_binding
            use mpi
            integer(c_int), value :: request
            integer(c_int), dimension(MPI_STATUS_SIZE) :: status
            integer(c_int), intent(out) :: ierr
        end subroutine my_MPI_Wait

        subroutine my_MPI_Isend(buf, count, datatype, dest, tag, Fcomm, request, ierr) bind(C, name="my_MPI_Isend")
            use iso_c_binding
            type(c_ptr), value :: buf
            integer(c_int), value :: count
            integer(c_int), value :: datatype
            integer(c_int), intent(in) :: dest
            integer(c_int), intent(in) :: tag
            integer(c_int), value :: Fcomm
            integer(c_int), intent(out) :: request
            integer(c_int), intent(out) :: ierr
        end subroutine my_MPI_Isend

        subroutine my_MPI_Irecv(buf, count, datatype, source, tag, Fcomm, request, ierr) bind(C, name="my_MPI_Irecv")
            use iso_c_binding
            type(c_ptr), value :: buf
            integer(c_int), value :: count
            integer(c_int), value :: datatype
            integer(c_int), intent(in) :: source
            integer(c_int), intent(in) :: tag
            integer(c_int), value :: Fcomm
            integer(c_int), intent(out) :: request
            integer(c_int), intent(out) :: ierr
        end subroutine my_MPI_Irecv

        subroutine my_MPI_Gather(sendbuf, sendcount, sendtype, recvbufer, recvcount, recvtype, &
             root, Fcomm, ierr) bind(C, name="my_MPI_Gather")
            use iso_c_binding
            type(c_ptr), value :: sendbuf
            integer(c_int), value :: sendcount
            integer(c_int), value :: sendtype
            type(c_ptr), value :: recvbufer
            integer(c_int), value :: recvcount
            integer(c_int), value :: recvtype
            integer(c_int), value :: root
            integer(c_int), value :: Fcomm
            integer(c_int), intent(out) :: ierr
        end subroutine my_MPI_Gather

        subroutine my_MPI_Allgather(sendbuf, sendcount, sendtype, recvbufer, recvcount, recvtype, &
             Fcomm, ierr) bind(C, name="my_MPI_Allgather")
            use iso_c_binding
            type(c_ptr), value :: sendbuf
            integer(c_int), value :: sendcount
            integer(c_int), value :: sendtype
            type(c_ptr), value :: recvbufer
            integer(c_int), value :: recvcount
            integer(c_int), value :: recvtype
            integer(c_int), value :: Fcomm
            integer(c_int), intent(out) :: ierr
        end subroutine my_MPI_Allgather

        subroutine my_MPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm, ierr) bind(c, name="my_MPI_Reduce")
            use iso_c_binding
            type(c_ptr), value :: sendbuf
            type(c_ptr), value :: recvbuf
            integer(C_INT), value :: count, root
            integer(C_INT), value :: datatype, op, comm
            integer(c_int), intent(out) :: ierr
        end subroutine my_MPI_Reduce

        subroutine my_MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm, ierr) bind(c, name="my_MPI_Allreduce")
            use iso_c_binding
            type(c_ptr), value :: sendbuf
            type(c_ptr), value :: recvbuf
            integer(C_INT), value :: count
            integer(C_INT), value :: datatype, op, comm
            integer(c_int), intent(out) :: ierr
        end subroutine my_MPI_Allreduce

    end interface
end module mpi_interface
