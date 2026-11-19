subroutine output_file
    use global
    implicit none
    
    integer :: i
    
    !Writing an singular outoput data file
    OPEN(FILE = 'output.txt', UNIT = 117)
    write(117,101)nnode,nele,ndrop
    write(117,*)
    write(117,118)dt,max_pos,max_it
    write(117,*)
    write(117,*)'NODES COORDINATES (i# node, x, y)'
    write(117,*)
    
    do i=1,nnode
        write(117,102)xnode(i,:)
    end do
    
    write(117,*)
    write(117,*)'CONNECTIVITY (i# elem, 1st nodes, 2nd node, 3th node, type)'
    write(117,*)
    
    do i=1,nele
        write(117,103)conn(i,:)
    end do
    
    write(117,*)
    write(117,*)'ELEM BOUNDARY CONDITIONS (i# elem, nd type x, vlaue x, nd type y, value y)'
    write(117,*)
    
    do i=1,nele
        write(117,104)bc(i,:)
    end do
    
    write(117,*)
    write(117,*)'DROPS FIRST AND LAST NODES (i# drop, first node, last node)'
    write(117,*)
    
    do i=1,ndrop
        write(117,105)drop_first_last_node(i,:)
    end do
    
    write(117,*)
    write(117,*)'DROPS PROPERTIES (Lambda, Ca, aH)'
    write(117,*)
    write(117,106)lambda, Ca, aH
    write(117,*)
    
    CLOSE(UNIT = 117)
    
    
112 FORMAT(3(F9.2,2X)) 
110 FORMAT(6(F9.2,2X))
101 format(3(i10),'        #number of nodes, elements and drops')
102 format(3(f18.7))
103 format(10(i10))
104 format(13(f18.7))
105 format(3(i7))
106 format(3(f18.7))
118 format((f18.7,f18.7,i10),'   #time step, max drop position and max iteration')
end subroutine