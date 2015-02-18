pro DEIMOS_raw_mesh, offset, file

   foo = intarr(8192, 8560-4*offset)
   aa = mrdfits(file,1)
   bb = mrdfits(file,2)
   cc = mrdfits(file,3)
   dd = mrdfits(file,4)
   ee = mrdfits(file,5)
   ff = mrdfits(file,6)
   gg = mrdfits(file,7)
   hh = mrdfits(file,8)
   AA = transpose(aa)
   BB = transpose(bb)
   CC = transpose(cc)
   DD = transpose(dd)
   EE = transpose(ee)
   FF = transpose(ff)
   GG = transpose(gg)
   HH = transpose(hh)
   AAflip = rotate(AA,2)
   BBflip = rotate(BB,2)
   CCflip = rotate(CC,2)
   DDflip = rotate(DD,2)
   EEflip = rotate(EE,2)
   FFflip = rotate(FF,2)
   GGflip = rotate(GG,2)
   HHflip = rotate(HH,2)
   foo[0:4095, 0:2139-offset] = AA[0:4095, 0:2139-offset]
   foo[4096:8191, 0:2139-offset] = EEflip[0:4095, offset:2139]
   foo[0:4095, 2139-offset+1:4279-2*offset] = BB[0:4095, 0:2139-offset]
   foo[4096:8191, 2139-offset+1:4279-2*offset] = FFflip[0:4095, offset:2139]
   foo[0:4095, 4279-2*offset+1:6419-3*offset] = CC[0:4095, 0:2139-offset]
   foo[4096:8191, 4279-2*offset+1:6419-3*offset] = GGflip[0:4095, offset:2139]
   foo[0:4095, 6419-3*offset+1:8559-4*offset] = DD[0:4095, 0:2139-offset]
   foo[4096:8191, 6419-3*offset+1:8559-4*offset] = HHflip[0:4095, offset:2139]

   atv, foo
   mwrfits, foo, 'DEIMOSraw4compare.fits'
end
