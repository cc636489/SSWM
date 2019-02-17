module convert

   !> \author Ali Samii - 2016
   !! \brief This object stores the data required for construction of a parallel or serial
   !! ESMF_Mesh from <tt>fort.14, fort.18, partmesh.txt</tt> files.
   !!
   type meshdata
      !> \details This array contains the node coordinates of the mesh. For
      !! example, in a 2D mesh, the \c jth coordinate of the \c nth node
      !! is stored in location <tt> 2*(n-1)+j</tt> of this array.
      real(kind=8), allocatable    :: NdCoords(:)
      !> \details This array contains the elevation of different nodes of the mesh
      real(kind=8), allocatable    :: bathymetry(:)
      real(kind=8), allocatable    :: elevation(:, :)
      real(kind=8), allocatable    :: uvelocity(:, :, :)
      real(kind=8), allocatable    :: timestamp(:)
      !> \details Number of nodes present in the current PE. This is different from the
      !! number of nodes owned by this PE (cf. NumOwnedNd)
      integer(kind=4)              :: NumNd
      integer(kind=4)              :: NumTime
      !> \details Number of elements in the current PE. This includes ghost elements and
      !! owned elements. However, we do not bother to distinguish between owned
      !! element and present element (as we did for the nodes).
      integer(kind=4)              :: NumEl
      !> \details Number of nodes of each element, which is simply three in 2D ADCIRC.
      integer(kind=4)              :: NumND_per_El
      !> \details Global node numbers of the nodes which are present in the current PE.
      integer(kind=4), allocatable :: NdIDs(:)
      !> \details Global element numbers which are present in the current PE.
      integer(kind=4), allocatable :: ElIDs(:)
      !> \details The element connectivity array, for the present elements in the current PE.
      !! The node numbers are the local numbers of the present nodes. All the element
      !! connectivities are arranged in this one-dimensional array.
      integer(kind=4), allocatable :: ElConnect(:)
      !> \details An array containing the element types, which are all triangles in our
      !! application.
      integer(kind=4), allocatable :: ElTypes(:)
      !> \details This flag tells if the meshdata has been initialized
      logical                      :: is_initialized = .false.
   end type meshdata

   contains

   !> \details This function creates an object of type meshdata from the fort14 file given by
   !!fort14_filename. Unlike extract_parallel_data_from_mesh(), this function does not
   !! create a parallel meshdata, so it can be called by only one PE and the created meshdata
   !! object can later be used to create an ESMF_Mesh object.
   subroutine extract_global_data_from_fort14(fort14_filename, the_data)
      implicit none
      type(meshdata), intent(inout)   :: the_data
      character(len=*), intent(in)    :: fort14_filename
      integer(kind=4)                 :: i1, i_num, rc
      integer(kind=4), parameter      :: dim1 = 2, spacedim = 2, NumND_per_El = 3
      logical                         :: iOpen, iExist

      open (unit=14, file=fort14_filename, form='FORMATTED', status='OLD', action='READ')
      read (unit=14, fmt=*)
      read (unit=14, fmt=*) the_data%NumEl, the_data%NumNd
      allocate (the_data%NdIDs(the_data%NumNd))
      allocate (the_data%ElIDs(the_data%NumEl))
      allocate (the_data%NdCoords(dim1*the_data%NumNd))
      allocate (the_data%bathymetry(the_data%NumNd))
      allocate (the_data%ElConnect(NumND_per_El*the_data%NumEl))
      allocate (the_data%ElTypes(the_data%NumEl))
      do i1 = 1, the_data%NumNd, 1
         read (unit=14, fmt=*) the_data%NdIDs(i1), &
            the_data%NdCoords((i1 - 1)*dim1 + 1), &
            the_data%NdCoords((i1 - 1)*dim1 + 2), &
            the_data%bathymetry(i1)
      end do
      do i1 = 1, the_data%NumEl, 1
         read (unit=14, fmt=*) the_data%ElIDs(i1), i_num, &
            the_data%ElConnect((i1 - 1)*NumND_per_El + 1), &
            the_data%ElConnect((i1 - 1)*NumND_per_El + 2), &
            the_data%ElConnect((i1 - 1)*NumND_per_El + 3)
      end do
      the_data%ElTypes = 3
      close (14)
   end subroutine

   !> read data from fort.63.
   subroutine extract_global_data_from_fort63(fort63_filename, the_data)
      implicit none
      type(meshdata), intent(inout)   :: the_data
      character(len=*), intent(in)    :: fort63_filename
      integer(kind=4)                 :: i1, j1, i_num, rc, a
      integer(kind=4), parameter      :: dim1 = 2, spacedim = 2, NumND_per_El = 3
      logical                         :: iOpen, iExist

      open (unit=63, file=fort63_filename, form='FORMATTED', status='OLD', action='READ')
      read (unit=63, fmt=*)
      read (unit=63, fmt=*) the_data%NumTime
      allocate (the_data%elevation(the_data%NumNd, the_data%NumTime))
      allocate (the_data%timestamp(the_data%NumTime))

      do j1 = 1, the_data%NumTime, 1
         read (unit=63, fmt=*) the_data%timestamp(j1)
      do i1 = 1, the_data%NumNd, 1
         read (unit=63, fmt=*) a, the_data%elevation(i1, j1)
      end do
      end do

      close (63)
   end subroutine

   !> read data from fort.64.
   subroutine extract_global_data_from_fort64(fort64_filename, the_data)
      implicit none
      type(meshdata), intent(inout)   :: the_data
      character(len=*), intent(in)    :: fort64_filename
      integer(kind=4)                 :: i1, j1, i_num, rc, a
      integer(kind=4), parameter      :: dim1 = 2, spacedim = 2, NumND_per_El = 3
      logical                         :: iOpen, iExist

      open (unit=64, file=fort64_filename, form='FORMATTED', status='OLD', action='READ')
      read (unit=64, fmt=*)
      read (unit=64, fmt=*)
      allocate (the_data%uvelocity(the_data%NumNd, spacedim, the_data%NumTime))

      do j1 = 1, the_data%NumTime, 1
         read (unit=64, fmt=*)
      do i1 = 1, the_data%NumNd, 1
         read (unit=64, fmt=*) a, the_data%uvelocity(i1, 1, j1), the_data%uvelocity(i1, 2, j1)
      end do
      end do

      close (64)
   end subroutine
   !> \details This function writes the input meshdata object to a \c vtu file.
   !! The \c vtu file is in \c XML format. This function can be used for both parallel
   !! and serial mesh writing. If one uses this function for parallel write, the
   !! processing element with \c localPE=0 should also enter this function, otherwise
   !! the \c pvtu file will not be written. This function assumes that the \c vtu file
   !! which we want to write does not exist. If we want to add fields to the files which
   !! are created before, we have to call write_node_field_to_vtu() function. If the user
   !! wants to add more data fields to the created \c vtu file, the \c last_write parameter
   !! should be passed <tt>.false.</tt> so that the program do not close the file.
   !! By closing we mean writing the last three closing lines in the XML \c vtu files.
   !! However, if this is the last time we want to write on the same \c vtu file, we
   !! have to pass \c last_write equal to <tt>.true.</tt>
   !! \param the_data This is the input data for which we create the vtu file
   !! \param vtu_filename This is the name of the vtu file
   !! \param last_write This parameter indicates if this is the last time we want to
   !! write something to this \c vtu file.
   subroutine write_meshdata_to_vtu(the_data, vtu_filename, last_write)
      implicit none
      type(meshdata), intent(in)      :: the_data
      character(len=*), intent(in)    :: vtu_filename
      logical, intent(in)             :: last_write
      integer(kind=4)                 :: i1, indent, offset_counter, rc, indent2
      integer(kind=4), parameter      :: dim1 = 2, spacedim = 2, NumND_per_El = 3, vtk_triangle = 5
      indent = 0

      open (unit=1014, file=vtu_filename, form="FORMATTED", &
            status="UNKNOWN", action="WRITE")
      write (unit=1014, fmt="(A)") '<?xml version="1.0"?>'
      write (unit=1014, fmt="(A,A)") '<VTKFile type="UnstructuredGrid"', &
         ' version="0.1" byte_order="LittleEndian">'
      indent = indent + 2
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         "<UnstructuredGrid>"
      indent = indent + 2
      write (unit=1014, fmt="(A,A,I0,A,I0,A)") repeat(" ", indent), &
         '<Piece NumberOfPoints="', the_data%NumNd, &
         '" NumberOfCells="', the_data%NumEl, '">'
      indent = indent + 2

      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '<Points>'
      indent = indent + 2
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '<DataArray type="Float32" NumberOfComponents="3" Format="ascii">'
      indent = indent + 2
      do i1 = 1, the_data%NumNd, 1
         indent2 = 1
         if (i1 == 1) then; indent2 = indent; endif
         write (unit=1014, fmt="(A,F0.4,' ',F0.4,' ',F0.4,' ')", advance='no') &
            repeat(" ", indent2), &
            the_data%NdCoords((i1 - 1)*dim1 + 1), &
            the_data%NdCoords((i1 - 1)*dim1 + 2), 0.0
      end do
      write (unit=1014, fmt="(A)") ""
      indent = indent - 2
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '</DataArray>'
      indent = indent - 2
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '</Points>'

      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '<Cells>'
      indent = indent + 2
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '<DataArray type="Int32" Name="connectivity" Format="ascii">'
      indent = indent + 2
      do i1 = 1, the_data%NumEl, 1
         indent2 = 1
         if (i1 == 1) then; indent2 = indent; endif
         write (unit=1014, fmt="(A,I0,' ',I0,' ',I0,' ')", advance='no') &
            repeat(" ", indent2), &
            the_data%ElConnect((i1 - 1)*NumND_per_El + 1) - 1, &
            the_data%ElConnect((i1 - 1)*NumND_per_El + 2) - 1, &
            the_data%ElConnect((i1 - 1)*NumND_per_El + 3) - 1
      end do
      write (unit=1014, fmt="(A)") ""
      indent = indent - 2
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '</DataArray>'
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '<DataArray type="Int32" Name="offsets" Format="ascii">'
      indent = indent + 2
      offset_counter = 0
      do i1 = 1, the_data%NumEl, 1
         indent2 = 1
         if (i1 == 1) then; indent2 = indent; endif
         offset_counter = offset_counter + 3
         write (unit=1014, fmt="(A,I0)", advance='no') repeat(" ", indent2), &
            offset_counter
      end do
      write (unit=1014, fmt="(A)") ""
      indent = indent - 2
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '</DataArray>'
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '<DataArray type="Int32" Name="types" Format="ascii">'
      indent = indent + 2
      do i1 = 1, the_data%NumEl, 1
         indent2 = 1
         if (i1 == 1) then; indent2 = indent; endif
         write (unit=1014, fmt="(A,I2)", advance='no') &
            repeat(" ", indent2), vtk_triangle
      end do
      write (unit=1014, fmt="(A)") ""
      indent = indent - 2
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '</DataArray>'
      indent = indent - 2
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '</Cells>'

      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '<CellData Scalars="scalars">'
      indent = indent + 2
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '<DataArray type="Int32" Name="subdomain_id" NumberOfComponents="1" Format="ascii">'
      indent = indent + 2
      do i1 = 1, the_data%NumEl, 1
         indent2 = 1
         if (i1 == 1) then; indent2 = indent; endif
         write (unit=1014, fmt="(A,I0)", advance='no') &
            repeat(" ", indent2), 0
      end do
      write (unit=1014, fmt="(A)") ""
      indent = indent - 2
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '</DataArray>'
      !
      indent = indent - 2
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '</CellData>'


      if (last_write) then
         indent = indent - 2
         write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
            '</Piece>'
         indent = indent - 2
         write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
            '</UnstructuredGrid>'
         indent = indent - 2
         write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
            '</VTKFile>'
      end if
      close (1014)
   end subroutine

   !> \details This function writes the input array (\c field_array) and its name (\c field_name)
   !! to the vtu file (which should already exist and not closed). Refer to write_meshdata_to_vtu()
   !! to know more about opening vtu file and closing them. If the parameter \c last_write is true
   !! then we close this file and as such we should not write anything else on this file.
   subroutine write_node_1D_field_to_vtu(field_array, field_name, vtu_filename, last_write)
      implicit none
      character(len=*), intent(in)   :: vtu_filename, field_name
      logical, intent(in)            :: last_write
      real(kind=8), intent(in)       :: field_array(:)
      integer(kind=4)                :: i1, indent, indent2, num_recs
      open (unit=1014, file=vtu_filename, form='FORMATTED', &
            position='APPEND', status='OLD', action='WRITE')

      indent = 8
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '<PointData Scalars="scalars">'
      indent = indent + 2
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '<DataArray type="Float32" Name="'//field_name//'" NumberOfComponents="1" Format="ascii">'
      indent = indent + 2
      num_recs = size(field_array)
      do i1 = 1, num_recs, 1
         indent2 = 1
         if (i1 == 1) then; indent2 = indent; endif
         write (unit=1014, fmt="(A,F0.4)", advance='no') &
            repeat(" ", indent2), field_array(i1)
      end do
      write (unit=1014, fmt="(A)") ""
      indent = indent - 2
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '</DataArray>'
      indent = indent - 2
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '</PointData>'

      if (last_write) then
         indent = indent - 2
         write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
            '</Piece>'
         indent = indent - 2
         write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
            '</UnstructuredGrid>'
         indent = indent - 2
         write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
            '</VTKFile>'
      end if
      close (1014)
   end subroutine

   subroutine write_node_2D_field_to_vtu(field_array, field_name, vtu_filename, last_write)
      implicit none
      character(len=*), intent(in)   :: vtu_filename, field_name
      logical, intent(in)            :: last_write
      real(kind=8), intent(in)       :: field_array(:, :)
      integer(kind=4)                :: i1, indent, indent2, num_recs
      open (unit=1014, file=vtu_filename, form='FORMATTED', &
            position='APPEND', status='OLD', action='WRITE')

      indent = 8
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '<PointData Vectors="vectors">'
      indent = indent + 2
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '<DataArray type="Float64" Name="'//field_name//'" NumberOfComponents="3" Format="ascii">'
      indent = indent + 2
      num_recs = size(field_array, 1)
      do i1 = 1, num_recs, 1
         indent2 = 1
         if (i1 == 1) then; indent2 = indent; endif
         write (unit=1014, fmt="(A,F0.4,A,F0.4,A,F0.4)", advance='no') &
            repeat(" ", indent2), field_array(i1, 1), repeat(" ", indent2), field_array(i1, 2), repeat(" ", indent2), 0.0
      end do
      write (unit=1014, fmt="(A)") ""
      indent = indent - 2
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '</DataArray>'
      indent = indent - 2
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '</PointData>'

      if (last_write) then
         indent = indent - 2
         write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
            '</Piece>'
         indent = indent - 2
         write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
            '</UnstructuredGrid>'
         indent = indent - 2
         write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
            '</VTKFile>'
      end if
      close (1014)
   end subroutine

   !> \details This function writes the input array (\c field_array) and its name (\c field_name)
   !! to the vtu file (which should already exist and not closed). Refer to write_meshdata_to_vtu()
   !! to know more about opening vtu file and closing them. If the parameter \c last_write is true
   !! then we close this file and as such we should not write anything else on this file.
   subroutine write_int_node_field_to_vtu(field_array, field_name, vtu_filename, last_write)
      implicit none
      character(len=*), intent(in)      :: vtu_filename, field_name
      logical, intent(in)               :: last_write
      integer(kind=4), intent(in)       :: field_array(:)
      integer(kind=4)                   :: i1, indent, indent2, num_recs
      open (unit=1014, file=vtu_filename, form='FORMATTED', &
            position='APPEND', status='OLD', action='WRITE')

      indent = 8
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '<PointData Vectors="vectors">'
      indent = indent + 2
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '<DataArray type="Int32" Name="'//field_name//'" NumberOfComponents="1" Format="ascii">'
      indent = indent + 2
      num_recs = size(field_array)
      do i1 = 1, num_recs, 1
         indent2 = 1
         if (i1 == 1) then; indent2 = indent; endif
         write (unit=1014, fmt="(A,I4)", advance='no') &
            repeat(" ", indent2), field_array(i1)
      end do
      write (unit=1014, fmt="(A)") ""
      indent = indent - 2
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '</DataArray>'
      indent = indent - 2
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '</PointData>'

      if (last_write) then
         indent = indent - 2
         write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
            '</Piece>'
         indent = indent - 2
         write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
            '</UnstructuredGrid>'
         indent = indent - 2
         write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
            '</VTKFile>'
      end if
      close (1014)
   end subroutine

   ! write corresponding pvd collection file.
   subroutine write_series_to_pvd(pvd_filename, time_array)
   
      implicit none
      character(len=*), intent(in)    :: pvd_filename
      character(len=80)               :: vtu_filename, temp, time
      real(kind=8), intent(in)        :: time_array(:)
      integer(kind=4)                 :: k1, indent, num_timestep
    
      num_timestep = size(time_array)

      indent = 0

      open (unit=1014, file=pvd_filename, form="FORMATTED", &
            status="UNKNOWN", action="WRITE")

      write (unit=1014, fmt="(A)") '<?xml version="1.0"?>'
      write (unit=1014, fmt="(A,A,A,A)") '<VTKFile type="Collection"', ' version="0.1"', &
          ' byte_order="LittleEndian"', ' compressor="vtkZLibDataCompressor" >'
      indent = indent + 2
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), "<Collection>"
      indent = indent + 2


      write(temp, fmt='(I7.7)') 0
      write(time, fmt='(f0.8)') 0.0
      write(vtu_filename, fmt="(A)") pvd_filename(:len(pvd_filename)-4)//trim(temp)//'.vtu'
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
          '<DataSet timestep="'//trim(time)//'" part="0" file="'//trim(vtu_filename)//'" />'

      Do k1=1, num_timestep, 1
      write(temp, fmt='(I7.7)') k1
      write(time, fmt='(f0.8)') time_array(k1)
      write(vtu_filename, fmt="(A)") pvd_filename(:len(pvd_filename)-4)//trim(temp)//'.vtu'
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
          '<DataSet timestep="'//trim(time)//'" part="0" file="'//trim(vtu_filename)//'" />'
      End do

      indent = indent - 2
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), '</Collection>'
      write (unit=1014, fmt="(A)") '</VTKFile>'
      

      close (1014)
   
   end subroutine


   ! write time series 1D data to series of vtu
   subroutine write_node_series_1D_field_to_vtu_pvd(the_data, field_array, time_array, field_name, pvd_filename)

      implicit none
      type(meshdata), intent(in)      :: the_data
      character(len=*), intent(in)    :: pvd_filename, field_name
      character(len=80)               :: vtu_filename, temp, time
      real(kind=8), intent(in)        :: field_array(:, :)
      real(kind=8), intent(in)        :: time_array(:)
      real(kind=8), allocatable       :: zero_array(:)
      integer(kind=4)                 :: i1, j1, k1, indent, num_timestep

      num_timestep = size(time_array)

      allocate(zero_array(the_data%NumNd))
      zero_array(:) = 0.0

      ! write the .vtu files in all time step. 

      write(temp, fmt='(I7.7)') 0
      write(vtu_filename, fmt="(A)") pvd_filename(:len(pvd_filename)-4)//trim(temp)//'.vtu'
      print *, trim(vtu_filename), num_timestep
      call write_meshdata_to_vtu(the_data, vtu_filename, .false.)
      call write_node_1D_field_to_vtu(zero_array, field_name, vtu_filename, .true.)

      Do j1=1, num_timestep, 1
     
      write(temp, fmt='(I7.7)') j1
      write(vtu_filename, fmt="(A)") pvd_filename(:len(pvd_filename)-4)//trim(temp)//'.vtu'
      print *, trim(vtu_filename), num_timestep

      call write_meshdata_to_vtu(the_data, vtu_filename, .false.)
      call write_node_1D_field_to_vtu(field_array(:, j1), field_name, vtu_filename, .true.)
      
      end do

      ! write .pvd file for collecting vtu. 
      call write_series_to_pvd(pvd_filename, time_array)

   end subroutine

   ! write time series 2D data to series of vtu
   subroutine write_node_series_2D_field_to_vtu_pvd(the_data, field_array, time_array, field_name, pvd_filename)

      implicit none
      type(meshdata), intent(in)      :: the_data
      character(len=*), intent(in)    :: pvd_filename, field_name
      character(len=80)               :: vtu_filename, temp, time
      real(kind=8), intent(in)        :: field_array(:, :, :)
      real(kind=8), intent(in)        :: time_array(:)
      real(kind=8), allocatable       :: zero_array(:, :)
      integer(kind=4)                 :: i1, j1, k1, indent, num_timestep

      num_timestep = size(time_array)

      allocate(zero_array(the_data%NumNd, 2))
      zero_array(:, :) = 0.0

      ! write the .vtu files in all time step. 

      write(temp, fmt='(I7.7)') 0
      write(vtu_filename, fmt="(A)") pvd_filename(:len(pvd_filename)-4)//trim(temp)//'.vtu'
      print *, trim(vtu_filename), num_timestep
      call write_meshdata_to_vtu(the_data, vtu_filename, .false.)
      call write_node_2D_field_to_vtu(zero_array, field_name, vtu_filename, .true.)

      Do j1=1, num_timestep, 1
     
      write(temp, fmt='(I7.7)') j1
      write(vtu_filename, fmt="(A)") pvd_filename(:len(pvd_filename)-4)//trim(temp)//'.vtu'
      print *, trim(vtu_filename), num_timestep

      call write_meshdata_to_vtu(the_data, vtu_filename, .false.)
      call write_node_2D_field_to_vtu(field_array(:, :, j1), field_name, vtu_filename, .true.)
      
      end do

      ! write .pvd file for collecting vtu. 
      call write_series_to_pvd(pvd_filename, time_array)

   end subroutine

end module convert

program main

  use convert

  type(meshdata) :: global_data

  call extract_global_data_from_fort14("fort.14", global_data)
  call extract_global_data_from_fort63("fort.63", global_data)
  call extract_global_data_from_fort64("fort.64", global_data)
  call write_meshdata_to_vtu(global_data, "global_mesh.vtu", .false.)
  call write_node_1D_field_to_vtu(global_data%bathymetry, &
                                   "bath52", "global_mesh.vtu", .true.)
 
  call write_node_series_1D_field_to_vtu_pvd(global_data, global_data%elevation, global_data%timestamp, &
                                   "elev", "elev.pvd")
  call write_node_series_2D_field_to_vtu_pvd(global_data, global_data%uvelocity, global_data%timestamp, &
                                   "u", "u.pvd")

end program main

