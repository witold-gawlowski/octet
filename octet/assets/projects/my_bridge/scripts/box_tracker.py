#Blender has to be run "by administrator" for this to work.
#Output file is saved in the Blender directory.
import bpy
import csv
import os, sys
#this is unfortunately not working for me...
#blendfilepath = bpy.data.filepath
#directory = os.path.dirname(blendfilepath)

#Cubes that make colider are called "cc.001" etc, whereas
#bridge cubes are called "bc.001" etc.
file_name = "box_list.csv"
file_ = open(file_name, 'w')
wr = csv.writer(file_, quoting=csv.QUOTE_NONE, delimiter = ' ')
colliders = [item for item in list(bpy.data.objects) if "box_collider" in item.name]
bridge_elems = [item for item in list(bpy.data.objects) if "bridge_box" in item.name]
#colliders first
wr.writerow([len(colliders)])
for b in colliders:
  wr.writerow([b.location[0], b.location[1], b.location[2], b.rotation_quaterion[0], b.rotation_quaterion[1], b.rotation_quaterion[2], b.rotation_quaterion[3], b.scale[0], b.scale[1], b.scale[2]])
#now bridge elements

wr.writerow([len(bridge_elems)])

for b in bridge_elems:
  wr.writerow([b.location[0], b.location[1], b.location[2], b.rotation_quaterion[0], b.rotation_quaterion[1], b.rotation_quaterion[2], b.rotation_quaterion[3], b.scale[0], b.scale[1], b.scale[2]])

file_.close()
