# Miscellaneous useful routines go here.

# int procedure	line_count (fd)	-- returns count
# int procedure get_valr (str, fmt, curr, new, vmin, vmax) -- returns stat
# int procedure get_vali (str, fmt, curr, new, vmin, vmax) -- returns stat


#
# LINE_COUNT: (I know I have another version of this around somewhere ...)
#

int procedure	line_count (fd)

pointer	fd			# file descriptor of open file

char	tchar
int	ndx

int	fscan(), nscan()
begin

# Count the entries
	ndx = 0
	while (fscan (fd) != EOF) {
		call gargwrd (tchar, 1)
		if (tchar == '#' || nscan() == 0) {
			next
		}
		ndx = ndx + 1
	}
	call seek (fd, BOF)

	return (ndx)
end

#
# GET_VALR: get a real value, allowing a default
#

int	procedure get_valr (str, fmt, curr, new, vmin, vmax)

char	str[ARB], fmt[ARB]			# value description, format
real	curr					# current value
real	new					# new value
real	vmin, vmax				# acceptible range

real	val

char	fmtline[80]
int	fscan(), nscan()

begin
	call sprintf (fmtline, 80, "%s (%s): ")
		call pargstr (str)
		call pargstr (fmt)
	call printf (fmtline)
		call pargr (curr)
	call flush (STDOUT)
	if (fscan(STDIN) != EOF) {
		call gargr (val)
		if (nscan() < 1) {
			new = curr
			return (OK)
		} else if (val <= vmax && val >= vmin) {
			new = val
			return (OK)
		} else {
			call eprintf ("Value out of range!! [%f,%f]\007")
				call pargr (vmin)
				call pargr (vmax)
			return (ERR)
		}
	}
	return (ERR)
end

#
# GET_VALI: get an integer  value, allowing a default
#


int	procedure get_vali (str, fmt, curr, new, vmin, vmax)

char	str[ARB], fmt[ARB]			# value description, format
int	curr					# current value
int	new					# new value
int	vmin, vmax				# acceptible range

int	val

char	fmtline[80]
int	fscan(), nscan()

begin
	call sprintf (fmtline, 80, "%s (%s): ")
		call pargstr (str)
		call pargstr (fmt)
	call printf (fmtline)
		call pargi (curr)
	call flush (STDOUT)
	if (fscan(STDIN) != EOF) {
		call gargi (val)
		if (nscan() < 1) {
			new = curr
			return (OK)
		} else if (val <= vmax && val >= vmin) {
			new = val
			return (OK)
		} else {
			call eprintf ("Value out of range!! [%d,%d]\007")
				call pargr (vmin)
				call pargr (vmax)
			return (ERR)
		}
	}
	return (ERR)
end

