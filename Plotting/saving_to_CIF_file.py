import datetime
from pathlib import Path
import numpy as np

PROJECT_DIR = Path(__file__).parent

def save_ice_structure_to_CIF(water_particles,
                              net_size,
                              ice_type='1c',
                              cell_size=3.566986):
    if ice_type=='1c':
        file_name='CIF\\'+'Ice_1c_'+str(net_size**3*8)+datetime.datetime.now().strftime('_oxygens_%d.%m.%y_%H_%M_%S.')+str(np.random.randint(100,999))+'.cif'
        path1 = PROJECT_DIR / 'Ice_1c_template.txt'
    elif ice_type=='6':
        file_name='CIF\\'+'Ice_6_'+str(net_size**3*5*2)+datetime.datetime.now().strftime('_oxygens_%d.%m.%y_%H_%M_%S.')+str(np.random.randint(100,999))+'.cif'
        path1 = PROJECT_DIR / 'Ice_6_template.txt'
    elif ice_type=='7':
        file_name='CIF\\'+'Ice_7_'+str(net_size**3*8*2)+datetime.datetime.now().strftime('_oxygens_%d.%m.%y_%H_%M_%S.')+str(np.random.randint(100,999))+'.cif'
        path1 = PROJECT_DIR / 'Ice_7_template.txt'
    else:
        print('Ice type unknown. Try 1c, 6 or 7.')
        
    divider=4*net_size
    cif_file_content=[]
    with open(path1, 'r') as template:
        lines = template.readlines()
        for line in lines:
            if '_cell_length_a' in line:
                cif_file_content.append('_cell_length_a '+str(float(line.strip('_cell_length_a '))*net_size)+'\n')
            elif '_cell_length_b' in line:
                cif_file_content.append('_cell_length_b '+str(float(line.strip('_cell_length_b '))*net_size)+'\n')
            elif '_cell_length_c' in line:
                cif_file_content.append('_cell_length_c '+str(float(line.strip('_cell_length_c '))*net_size)+'\n')
            elif '_cell_volume' in line:
                cif_file_content.append('_cell_volume '+str(float(line.strip('_cell_volume '))*net_size**3)+'\n')
            else:
                cif_file_content.append(line)
                
        hydrogen_atoms=[]
        for i, water_particle in enumerate(water_particles):
            hydrogen_atoms.append(water_particle[-3])

        j=0
        for i, water_particle in enumerate(water_particles):
            oxygen_atom_x=water_particle[0][0]
            oxygen_atom_y=water_particle[0][1]
            oxygen_atom_z=water_particle[0][2]

            oxygen_coords=np.array(water_particles[i][0])/divider
            cif_file_content.append(str('O{}'.format(i)+'\t'+str(oxygen_coords[0])+'\t'+str(oxygen_coords[1])+'\t'+str(oxygen_coords[2]))+'\n')

            H_O_d_1=0.25        #x,y,z coord lenght of OH bond in ice 1c
            H_O_d_2=1-H_O_d_1   #x,y,z coord lenght of H2O to H2O particle bond in ice 1c
           
            if water_particle[0][2]%2==0:
                if hydrogen_atoms[i][0]:
                    hydrogen_coords=(oxygen_atom_x-H_O_d_1, oxygen_atom_y-H_O_d_1, oxygen_atom_z+H_O_d_1)
                else:
                    hydrogen_coords=(oxygen_atom_x-H_O_d_2, oxygen_atom_y-H_O_d_2, oxygen_atom_z+H_O_d_2)
                cif_file_content.append(str('H{}'.format(j)+'\t'+str(hydrogen_coords[0]/divider)+'\t'+str(hydrogen_coords[1]/divider)+'\t'+str(hydrogen_coords[2]/divider))+'\n')
                j+=1
                if hydrogen_atoms[i][1]:
                    hydrogen_coords=(oxygen_atom_x+H_O_d_1, oxygen_atom_y+H_O_d_1, oxygen_atom_z+H_O_d_1)
                else:
                    hydrogen_coords=(oxygen_atom_x+H_O_d_2, oxygen_atom_y+H_O_d_2, oxygen_atom_z+H_O_d_2)
                cif_file_content.append(str('H{}'.format(j)+'\t'+str(hydrogen_coords[0]/divider)+'\t'+str(hydrogen_coords[1]/divider)+'\t'+str(hydrogen_coords[2]/divider))+'\n')
                j+=1
            if water_particle[0][2]%2==1:
                if hydrogen_atoms[i][0]:
                    hydrogen_coords=(oxygen_atom_x-H_O_d_1, oxygen_atom_y+H_O_d_1, oxygen_atom_z+H_O_d_1)
                else:
                    hydrogen_coords=(oxygen_atom_x-H_O_d_2, oxygen_atom_y+H_O_d_2, oxygen_atom_z+H_O_d_2)
                cif_file_content.append(str('H{}'.format(j)+'\t'+str(hydrogen_coords[0]/divider)+'\t'+str(hydrogen_coords[1]/divider)+'\t'+str(hydrogen_coords[2]/divider))+'\n')
                j+=1
                if hydrogen_atoms[i][1]:
                    hydrogen_coords=(oxygen_atom_x+H_O_d_1, oxygen_atom_y-H_O_d_1, oxygen_atom_z+H_O_d_1)
                else:
                    hydrogen_coords=(oxygen_atom_x+H_O_d_2, oxygen_atom_y-H_O_d_2, oxygen_atom_z+H_O_d_2)
                cif_file_content.append(str('H{}'.format(j)+'\t'+str(hydrogen_coords[0]/divider)+'\t'+str(hydrogen_coords[1]/divider)+'\t'+str(hydrogen_coords[2]/divider))+'\n')
                j+=1
            
    path2= PROJECT_DIR / file_name
    (PROJECT_DIR / 'CIF').mkdir(exist_ok=True)    #git nie przechowuje pustych katalogów, więc po sklonowaniu CIF/ może nie istnieć
    with open(path2, 'w') as cif_file:
        for cif_line in cif_file_content:
            cif_file.write(cif_line)
    print(file_name, 'saved')
    print(path2)
    return file_name

def save_ice_7_structure_to_CIF(water_particles1, waterparticles2, net_size,
                                 file_name='Ice_7_'+datetime.datetime.now().strftime('%d.%m.%y_%H:%M:%S')):
    path = PROJECT_DIR / 'Ice_7_template.txt'
    template=open(path, 'r')
    print('oho')
    return True
