
from simutil import simutil

sim=simutil()

cx=5108395.07973599
cy=2006464.69354945
cz=-3238571.37787314

posskamid = me.position("wgs84", "21.443803deg", "-30.712925deg", "1053.000000m")

config_file = "shared/ska1mid.cfg"
tel_name = "MID"
f = open(config_file, "r")
fo = open("shared/ska1mid_local.cfg", "w")
while True:
    line = f.readline()
    if not line: break
    items = line.split()
    if (items[0] != "#"):
        xx = float(items[0])
        yy = float(items[1])
        zz = float(items[2])
        ss = str(items[4])
        if ss[0] == "M":
            dd = float(13.5)
        else:
            dd = float(15.0)
        lx, ly, lz = simutil.itrf2loc(sim, x=xx, y=yy, z=zz, cx=cx, cy=cy, cz=cz)
        line = "%.5f %.5f %.5f %.2f %s\n" % (lx[0], ly[0], lz[0], dd, ss)
        fo.write(line)
f.close()
fo.close()
