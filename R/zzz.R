# $Id: zzz.R,v 1.3 2006/12/29 21:29:57 edd Exp $
.First.lib <- function(lib, pkg) {
  library.dynam("rvmbinary", pkg, lib )
}

.Last.lib <- function(libpath){
	library.dynam.unload("rvmbinary",libpath)
}