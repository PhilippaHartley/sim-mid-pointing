
from simutil import simutil
sim=simutil()
d2r=pi/180.0

cx, cy, cz = simutil.long2xyz(sim, 116.4525771*d2r, -26.60055525*d2r, 300.0, "WGS84")

config_file = "/Users/timcornwell/Code/algorithm-reference-library/data/configurations/LOW_SKA-TEL-SKO-0000422_Rev3.txt"
tel_name = "low"
f = open(config_file, "r")
fo = open("/Users/timcornwell/Code/algorithm-reference-library/data/configurations/ska1low_local.cfg", "w")
stationid=0
while True:
    line = f.readline()
    if not line: break
    items = line.split()
    if (items[0] != "#"):
        dd = float(items[0])
        lo = float(items[1])
        la = float(items[2])
        ss = "SKA1LOW%03d" % stationid
        stationid += 1
        xx, yy, zz = simutil.long2xyz(sim, lo * d2r, la * d2r, 0.0, "WGS84")
        lx, ly, lz = simutil.itrf2loc(sim, x=xx, y=yy, z=zz, cx=cx, cy=cy, cz=cz)
        line = "%.5f %.5f %.5f %.2f %s\n" % (lx[0], ly[0], lz[0], dd, ss)
        fo.write(line)
f.close()
fo.close()
