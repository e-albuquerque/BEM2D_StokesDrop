subroutine meshdisc
    use global
    
    call calc_maxdim !calculates the max dimension of the problem
    call create_nodes_conn !crates the connective and the nodes coordinates
    
end subroutine meshdisc