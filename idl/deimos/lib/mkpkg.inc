# I don't know what this is all about....
# Global MKPKG definitions for the NEWXXXX package.

$set XFLAGS	= "$(XFLAGS) -p deimos"
$set XVFLAGS	= "$(XVFLAGS) -p deimos"
$set LFLAGS	= "$(LFLAGS) -p deimos"

# Uncomment and modify the following to add special file list entries for
# various machine architectures and Fortran compilers.

# $ifeq (MACH, sparc) then
# $include "deimos$lib/mkpkg.sf.sun4"
# $else $ifeq (MACH, vms) then
# $include "deimos$lib/mkpkg.sf.vms"
# $end
