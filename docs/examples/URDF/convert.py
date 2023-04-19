import pymeshlab
import os
os.chdir('examples\\URDF\\anymal_b_simple_description\\meshes')
os.chdir('..\\..\\..\\..')
os.chdir('examples\\URDF\\unitree_ros\\robots\\a1_description\\meshes')
os.chdir('examples\\URDF\\example-robot-data\\robots\\anymal_b_simple_description\\meshes')
os.chdir('examples\\URDF\\example-robot-data\\robots\\panda_description\\meshes\\visual')
filenames = os.listdir()
ms = pymeshlab.MeshSet()
def convert_dae(filenames):
    for f in filenames:
        if os.path.isdir(f):
            os.chdir(f)
            subdir_filenames = os.listdir()
            convert_dae(subdir_filenames)
            os.chdir('..')
        else:
            fn,ext = os.path.splitext(f)
            if ext == '.dae':
                ms.load_new_mesh(f)
                ms.save_current_mesh(fn+'.obj')

convert_dae(filenames)


