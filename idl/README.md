# Install the Pipeline
## Get Dependencies

This pipeline is heavily dependent on IDL and Xwindows. It has been run on IDL 7.06 and 7.1 but may well work with later versions.

## Download Code

The code in this repo is currently setup to run on a Mac running OSX Snow Leapord or later.  If you have a different operating system then you will need to follow the 
DEEP2 pipeline installation instructions:

http://www2.keck.hawaii.edu/inst/deimos/pipeline.html

## Setup IDL Environment

### Edit your .bash_profile file (or similar)

Copy the following lines to your .bash_profile file:

```
# Setup idl
#setup cvs
export CVS_RSH=ssh
export CVS_DIR=$HOME/Git/DEIMOS/idl
#setup idl and related utilities
source ~/.idlenv
```

### Create/edit your .idlenv file

If you don’t alread have a .idlenv file then copy the one in this repository to your $HOME folder location, or where ever you have decided to source it in your .bash_profile.

If you already have one then copy the relevant information from the .idlenv file in this repo to your’s.

### Create/edit your .idlstartup file

If you don’t alread have a . idlstartup file then copy the one in this repository to your $HOME folder location, or where ever you have decided to source it in your .idlenv.

If you already have one then copy the relevant information from the . idlstartup file in this repo to your’s.

# Reduce the Data
1. Place the raw data in a folder with the date of the observation (e.g. 2015feb16) in the $DEIMOS_DATA directory specified in the .idlenv file.
2. For a given mask create a folder with that mask short name and a subfolder with the date of the observation in the $D2_RESULTS directory (e.g. $D2_RESULTS/a5233A/2015feb16)
3. cd to that directory and create a plan file (see e.g. http://deep.ps.uci.edu/spec2d/xxxx.plan.html)
5. start idl
6. run the following command:
> domask_2013a,nlsky=0
Note that this version of domask is good for masks observed during or after the 2013a semester. If the pipeline fails during the reduction proceedure it may be do to a bad slit. See the Skipping Bad Slits section below.

After the this portion of the reduction is complete then run the spec1d code in IDL to generate the zspec files.
> reduce1d,nlsky=0,doplot=0 

## Skipping Bad Slits
If bad slits are encountered during the reduction process then it might be necessary to create a slitcat file (listing all of the valid slits, i.e. deleting the bad slit number from the list).
So, to make the initial slitcat, open the "bintabs" file in the reduction directory, which is created after you start the reduction process (so I've assumed here that you've already tried to run the reduction once and it's failed), by doing in IDL:

> m = mrdfits('*bintabs*',3)

Then create the slitcat:

> openw, lun, 'XXXslits.cat', /get_lun
> for i=0, n_elements(m.slitname)-1 do printf, lun, m[i].slitname
> free_lun, lun

Then go into the slitcat and remove the lines that correspond to the slits that have failed. After that you add to your normal IDL command line:

> slitcat = 'XXXslits.cat'
> domask... slitcat=slitcat

# Run zspec
> cd $D2_RESULTS
> zspec,'maskname',/orelse

where 'maskname' is the short deimos mask name (e.g. a5233A)
