import sys

fname = sys.argv[1]
fname2 = sys.argv[2]

f = open(fname)
f2 = open(fname2)

#create dictionary for f2
f2dict = {}
for l2 in f2:
    l2 = l2.strip()
    f2dict[l2] = 1
f2.close()

#check f
out = open(fname[:-4]+fname2[:-4]+"_overlap.cmp","w")
for l in f:
    l =l.strip()
    if l in f2dict and l.find('jet_') != -1:
        out.write(l+"\n")

out.close()
f.close()
