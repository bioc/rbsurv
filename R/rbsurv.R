##########################################################################
#
#        
#    Robust Likelihood-based Survival Modeling for High-throughput Data
#
#             by
#
#      HyungJun Cho, Sukwoo Kim, and Soo-heang Eo
#      Deparments of Statistics and Biostatistics
#      Korea University
#
#        June 2008
#
##########################################################################

.First.lib <- function(lib, pkg) {  
     invisible()
     if(.Platform$OS.type=="windows" && require(Biobase) && interactive() 
        && .Platform$GUI=="Rgui") { addVigs2WinMenu("rbsurv") }
}

'rbsurv' <- function(time, ...)
{
    UseMethod("rbsurv")
}

