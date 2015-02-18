pro uberdattofits

	;readcol, 'full.fromcat.merged.ACS+LFC+RADIO+XRAY+SPECT+WIRC+SPITZER.dat', format='I,I,I,I,I,I,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D', mask_ACS_LFC, mask_AL_RX, mask_spect, mask_wirc_opt, mask_irac_mips, mask_spitz_opt, RA, dec,   RA_acs, dec_acs, RA_lfc, dec_lfc, RA_radio, dec_radio, RA_xray,  dec_xray, RA_wirc, dec_wirc, RA_spitz, dec_spitz, xs0, xs1,numline=3 		;readcol will only do this many columns, using the .sav file instead
	
	restore, 'ubercatV2.sav'
	

	help, Q, /struct
end
