subroutine getfilename
    use global
    implicit none
    
    open(unit=1, file='file_name.dat',status='old')
    
    read(unit=1, fmt=*)inputfilename !read the name of the input file 
    
    close(unit=1, status='keep')
    
    return
end subroutine getfilename