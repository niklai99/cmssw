# example usage
# python scoutingToRaw_bmtf.py /eos/user/m/matteomi/scoutingData/run362653/scout_GMT_362653_000506.dat SDSID ramdisk/run362653/ 362653 1000

import struct
import os,sys
import json
import shutil

os.umask(0)

#need to read this to find orbit ("event") boundary and calculate size per orbit
class header_shifts:
  bxmatch    = 32
  mAcount    = 16
  orbitmatch = 8
  mBcount    = 0

class header_masks:
  bxmatch    = 0xff << header_shifts.bxmatch
  mAcount    = 0xf  << header_shifts.mAcount
  orbitmatch = 0xff << header_shifts.orbitmatch
  mBcount    = 0xf


#new V2 FRD file header (32 bytes)
class frd_file_header_v2:
  ver_id = "RAW_0002".encode() # 64 (offset 0B)
  header_size = 32 #16 (offset 8B)
  data_type = 20 #16 (offset 10)
  event_count = 0 #32 (offset 12B)
  run_number = 0 #32 (offset 16B)
  lumisection = 0 #32 (offset 20B)
  file_size = 0 #64 (offset 24B)


def parseBmtfScoutingRawFile(infilepath, sdsid, outdir, rn_override, maxorbits):

  if infilepath != 'stdin':
    fin = open(infilepath,'rb')
  else:
    fin = sys.stdin.buffer

  #orbit count per file
  orbitcount=0
  #total
  orbitcount_total=0

  last_ls = 0

  orbit_data = bytes()
  orbit_nr = 0
  orbit_size = 0
  flags = 0
  c_crc32c = 0

  #ls = 1
  #event header (FRD format) const
  version = 6

  #files
  fout = None
  if infilepath != 'stdin':
    fin = open(infilepath,'rb')
  else:
    fin = sys.stdin.buffer


  #write header before closing the file
  def update_header():
    nonlocal orbitcount
    nonlocal last_ls
    h = frd_file_header_v2()
    h.event_count = orbitcount
    h.run_number = rn_override
    h.lumisection = last_ls
    h.file_size = fout.tell()
    fout.seek(0, 0)
    fout.write(frd_file_header_v2.ver_id)
    fout.write(struct.pack('H',h.header_size))
    fout.write(struct.pack('H',h.data_type))
    fout.write(struct.pack('I',h.event_count))
    fout.write(struct.pack('I',h.run_number))
    fout.write(struct.pack('I',h.lumisection))
    fout.write(struct.pack('Q',h.file_size))

    orbitcount = 0
    print(h.ver_id, h.header_size, h.data_type, h.event_count, h.lumisection, h.file_size)


  #write orbit when next one is detected or file is closed
  def write_orbit():
    nonlocal orbit_size
    nonlocal orbit_data
    if not orbit_size:
        return

    #print(fout.tell(), struct.pack('H',version))
    fout.write(struct.pack('H',version)) #could be 8 bytes
    fout.write(struct.pack('H',flags)) #could be 8 bytes
    fout.write(struct.pack('I',rn_override)) #run
    #fout.write(struct.pack('I',ls)) #ls
    fout.write(struct.pack('I',last_ls)) #ls
    fout.write(struct.pack('I',orbit_nr)) #eid (orbit number, 32-bit)
    fout.write(struct.pack('I',orbit_size+4)) #payload size +4 for source id
    fout.write(struct.pack('I',c_crc32c)) #payload checksum (not used)

    fout.write(struct.pack('I', sdsid)) # write the BMTF source ID

    fout.write(orbit_data)

    orbit_data = bytes()
    orbit_size = 0

  def writeout_close():
    write_orbit()
    update_header()
    fout.close()
    orbit_nr = 0

  #read loop
  endBmtfOrbit = False
  while True:

    #check if exceeded max orbits specified
    if (orbitcount_total > maxorbits) | (endBmtfOrbit):
      print(f"finish: {orbitcount_total}/{maxorbits} orbits")
      writeout_close()

      if infilepath != 'stdin':
        fin.close()
      sys.exit(0)

    try:
      h_raw = fin.read(4)

      if (h_raw==''):
        print('Reached end of file')
        writeout_close()
        sys.exit(0)

      bx     = struct.unpack('B', h_raw[2:4])[0]
      sCount = struct.unpack('B', h_raw[0:1])[0]

      #print("bx", bx, "sCount", sCount)

      orbit_raw = fin.read(4)
      orbit = struct.unpack('i', orbit_raw)[0] - 1 ## TEST

      if orbitcount_total >= maxorbits:
        if orbit != orbit_nr:
          endBmtfOrbit = True
          continue

      new_ls = orbit >> 18

      if new_ls > last_ls:
        #open a new output file if crossing LS boundary or on first orbit
        if last_ls:
          write_orbit()
          update_header()
          fout.close()
          orbitcount = 0

        last_ls = new_ls
        fout = open(os.path.join(outdir, f'run{rn_override}_ls{str(new_ls).zfill(4)}_index000000.raw') ,'wb')
        #empty file header, will be updated later
        fout.write(frd_file_header_v2.ver_id)
  #      fout.write(bytes(16))
        fout.write(bytes(24))

      read_len = 8*sCount
      stub_blk = fin.read(8*sCount)
      if len(stub_blk) != read_len:
        print('incomplete read')
        sys.exit(1)

      if not orbit_nr or orbit != orbit_nr:
        #received new orbit, write previous one
        if orbit_nr:
          write_orbit()

        if abs(orbit-orbit_nr)>1:
          print(f"\n\n{orbit_nr} -> {orbit}\n\n")

        #should not decrease:
        if orbit < orbit_nr:
          orbit_count = -1
          print("got smaller orbit than earlier!")
          sys.exit(1)

        # new valid orbit
        #print("new orbit", orbit)
        orbit_nr = orbit

        if (orbitcount_total%1000==0):
          print(f"Processed {orbitcount_total} orbits")

        #per LS file counter:
        orbitcount += 1
        #total counter:
        orbitcount_total += 1

      #update orbit size and data variables
      orbit_size += 8 + read_len
      orbit_data += (h_raw + bx_raw + orbit_raw) + stub_blk

    except Exception as ex:
      #reached premature end of the file?
      print(f"exception: {ex}")
      #writeout_close()
      #if infilepath != 'stdin':
      #    fin.close()
      sys.exit(1)

    #print count," : ",version,run,lumi,eid,esize,crc32c,"override id/ls/run:",count,1,rn_override
    #lumi=1

if len(sys.argv) < 6:
  print("parameters: input file (or stdin), output directory, run number (use same as input file), orbits to write")
else:
  parseBmtfScoutingRawFile(sys.argv[1], int(sys.argv[2]), sys.argv[3], int(sys.argv[4]), int(sys.argv[5]))