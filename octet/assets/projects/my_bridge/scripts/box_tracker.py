#Blender has to be run "by administrator" for this to work
import bpy
import csv
import os, sys
#this is unfortunately not working for me...
#blendfilepath = bpy.data.filepath
#directory = os.path.dirname(blendfilepath)
file_ = open("box_list.csv", 'w')
wr = csv.writer(file_, quoting=csv.QUOTE_NONE)
for b in list(bpy.data.objects):
	if "Cube" in b.name:
		print(b.location)
		wr.writerow([b.location[0], b.location[1], b.location[2] ])
		
file_.close()
