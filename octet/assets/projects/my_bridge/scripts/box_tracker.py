#Blender has to be run "by administrator" for this to work.
#Output file is saved in the Blender directory.
import bpy
import csv
import os, sys
#this is unfortunately not working for me...
#blendfilepath = bpy.data.filepath
#directory = os.path.dirname(blendfilepath)
file_ = open("box_list.csv", 'w')
wr = csv.writer(file_, quoting=csv.QUOTE_NONE, delimiter = ' ')
objs = [item for item in list(bpy.data.objects) if "Cube" in item.name]
wr.writerow([len(objs)])
for b in objs:
		wr.writerow([b.location[0], b.location[1], b.location[2], b.rotation_quaterion[0], b.rotation_quaterion[1], b.rotation_quaterion[2], b.rotation_quaterion[3], b.scale[0], b.scale[1], b.scale[2]])
		
file_.close()
