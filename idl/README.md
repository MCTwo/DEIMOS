# Get Dependencies

This pipeline is heavily dependent on IDL and Xwindows. It has been run on IDL 7.06 and 7.1 but may well work with later versions.

# Download Code

The code in this repo is currently setup to run on a Mac running OSX Snow Leapord or later.  If you have a different operating system then you will need to follow the 
DEEP2 pipeline installation instructions:

http://www2.keck.hawaii.edu/inst/deimos/pipeline.html

# Setup IDL Environment

## Edit your .bash_profile file (or similar)

Copy the following lines to your .bash_profile file:

```
# Setup idl
#setup cvs
export CVS_RSH=ssh
export CVS_DIR=$HOME/Git/DEIMOS/idl
#setup idl and related utilities
source ~/.idlenv
```

## Create/edit your .idlenv file

If you don’t alread have a .idlenv file then copy the one in this repository to your $HOME folder location, or where ever you have decided to source it in your .bash_profile.

If you already have one then copy the relevant information from the .idlenv file in this repo to your’s.

## Create/edit your .idlstartup file

If you don’t alread have a . idlstartup file then copy the one in this repository to your $HOME folder location, or where ever you have decided to source it in your .idlenv.

If you already have one then copy the relevant information from the . idlstartup file in this repo to your’s.