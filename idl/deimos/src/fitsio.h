#  Version Info: This file is distributed with version 2.401 of CFITSIO   */


# global variables */
 
define FLEN_FILENAME 1025 # max length of a filename  */
define FLEN_KEYWORD   72  # max length of a keyword (HIERARCH convention) */
define FLEN_CARD      81  # length of a FITS header card */
define FLEN_VALUE     71  # max length of a keyword value string */
define FLEN_COMMENT   73  # max length of a keyword comment string */
define FLEN_ERRMSG    81  # max length of a FITSIO error message */
define FLEN_STATUS    31  # max length of a FITSIO status text string */
 
define TBIT          1  # codes for FITS table data types */
define TBYTE        11
define TLOGICAL     14
define TSTRING      16
define TUSHORT      20
define TSHORT       21
define TUINT        30
define TINT         31
define TULONG       40
define TLONG        41
define TINT32BIT    41  # used when returning datatype of a column */
define TFLOAT       42
define TLONGLONG    81
define TDOUBLE      82
define TCOMPLEX     83
define TDBLCOMPLEX 163

define TYP_STRUC_KEY 10
define TYP_CMPRS_KEY 20
define TYP_SCAL_KEY  30
define TYP_NULL_KEY  40
define TYP_DIM_KEY   50
define TYP_RANG_KEY  60
define TYP_UNIT_KEY  70
define TYP_DISP_KEY  80
define TYP_HDUID_KEY 90
define TYP_CKSUM_KEY 100
define TYP_WCS_KEY   110
define TYP_REFSYS_KEY 120
define TYP_COMM_KEY  130
define TYP_CONT_KEY  140
define TYP_USER_KEY  150



define BYTE_IMG      8  # BITPIX code values for FITS image types */
define SHORT_IMG    16
define LONG_IMG     32
define LONGLONG_IMG 64
define FLOAT_IMG   -32
define DOUBLE_IMG  -64
                         # The following 2 codes are not true FITS         */
                         # datatypes; these codes are only used internally */
                         # within cfitsio to make it easier for users      */
                         # to deal with unsigned integers.                 */
define USHORT_IMG   20
define ULONG_IMG    40

define IMAGE_HDU  0  # Primary Array or IMAGE HDU */
define ASCII_TBL  1  # ASCII table HDU  */
define BINARY_TBL 2  # Binary table HDU */
define ANY_HDU   -1  # matches any HDU type */

define READONLY  0    # options when opening a file */
define READWRITE 1

