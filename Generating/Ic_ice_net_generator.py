import numpy as np
import copy
import secrets
#import energy_calculations #import calculate_particle_energy

def calculate_particle_energy(potential_net, water_particle, new_hydrogen_position=None, new_hydrogen_value=None):
    own_hydrogens_sum=0
    neighbors_hydrogen_sum=0
    x_coordinate_first_neighbor=water_particle[1][0][0]
    y_coordinate_first_neighbor=water_particle[1][0][1]
    z_coordinate_first_neighbor=water_particle[1][0][2]

    x_coordinate_second_neighbor=water_particle[1][1][0]
    y_coordinate_second_neighbor=water_particle[1][1][1]
    z_coordinate_second_neighbor=water_particle[1][1][2]
    #print('somsiad', x_coordinate_first_neighbor, y_coordinate_first_neighbor, z_coordinate_first_neighbor,
    #      potential_net[x_coordinate_first_neighbor][y_coordinate_first_neighbor][z_coordinate_first_neighbor],
    #      x_coordinate_second_neighbor, y_coordinate_second_neighbor, z_coordinate_second_neighbor,
    #      potential_net[x_coordinate_second_neighbor][y_coordinate_second_neighbor][z_coordinate_second_neighbor])
    if new_hydrogen_position==0:
        first_neighbor_hydrogens=new_hydrogen_value
        second_neighbor_hydrogens=potential_net[x_coordinate_first_neighbor][y_coordinate_first_neighbor][z_coordinate_first_neighbor][2][1]
    elif new_hydrogen_value==1:
        first_neighbor_hydrogens=potential_net[x_coordinate_second_neighbor][y_coordinate_second_neighbor][z_coordinate_second_neighbor][2][0]             
        second_neighbor_hydrogens=new_hydrogen_value

    #print(water_particle,
    #      first_neighbor_hydrogens,
    #      second_neighbor_hydrogens)

    own_hydrogens_sum+=sum(water_particle[2])
    neighbors_hydrogen_sum+=(first_neighbor_hydrogens+second_neighbor_hydrogens)
    energy_penalty=abs(own_hydrogens_sum+neighbors_hydrogen_sum-2)

    return energy_penalty


def generate_ice_1c_structure(net_size=1):

    
    oxygen_atoms=[]
    hydrogen_atoms=[]

    potential_net = np.zeros((net_size*4,net_size*4,net_size*4), dtype=list)    

    for x_iterator, net_layer in enumerate(potential_net):
      for y_iterator, net_row in enumerate(net_layer):
        for z_iterator, net_node in enumerate(net_row):
        
          if x_iterator%4==0 and y_iterator%4==0 and z_iterator%4==0:
            hydrogen_presence=[secrets.randbelow(2),secrets.randbelow(2)]
            oxygen_atoms.append([[x_iterator, y_iterator, z_iterator],
                                 [[x_iterator+1, y_iterator-1, z_iterator-1], 
                                  [x_iterator-1, y_iterator+1, z_iterator-1], 
                                  [x_iterator-1, y_iterator-1, z_iterator+1],
                                  [x_iterator+1, y_iterator+1, z_iterator+1]
                                  ],
                                  hydrogen_presence])
            hydrogen_atoms.append(hydrogen_presence)

          if x_iterator%4==2 and y_iterator%4==2 and z_iterator%4==0:
            hydrogen_presence=[secrets.randbelow(2),secrets.randbelow(2)]
            oxygen_atoms.append([[x_iterator, y_iterator, z_iterator],
                                 [[x_iterator+1, y_iterator-1, z_iterator-1], 
                                  [x_iterator-1, y_iterator+1, z_iterator-1], 
                                  [x_iterator-1, y_iterator-1, z_iterator+1],
                                  [x_iterator+1, y_iterator+1, z_iterator+1]
                                  ],
                                  hydrogen_presence])
            hydrogen_atoms.append(hydrogen_presence)

          if x_iterator%4==1 and y_iterator%4==1 and z_iterator%4==1:
            hydrogen_presence=[secrets.randbelow(2),secrets.randbelow(2)]
            oxygen_atoms.append([[x_iterator, y_iterator, z_iterator],
                                 [[x_iterator+1, y_iterator+1, z_iterator-1], 
                                  [x_iterator-1, y_iterator-1, z_iterator-1], 
                                  [x_iterator-1, y_iterator+1, z_iterator+1],
                                  [x_iterator+1, y_iterator-1, z_iterator+1]
                                  ],
                                  hydrogen_presence])
            hydrogen_atoms.append(hydrogen_presence)

          if x_iterator%4==0 and y_iterator%4==2 and z_iterator%4==2:
            hydrogen_presence=[secrets.randbelow(2),secrets.randbelow(2)]
            oxygen_atoms.append([[x_iterator, y_iterator, z_iterator],
                                 [[x_iterator+1, y_iterator-1, z_iterator-1], 
                                  [x_iterator-1, y_iterator+1, z_iterator-1], 
                                  [x_iterator-1, y_iterator-1, z_iterator+1],
                                  [x_iterator+1, y_iterator+1, z_iterator+1]
                                  ],
                                  hydrogen_presence])
            hydrogen_atoms.append(hydrogen_presence)

          if x_iterator%4==2 and y_iterator%4==0 and z_iterator%4==2:
            hydrogen_presence=[secrets.randbelow(2),secrets.randbelow(2)]
            oxygen_atoms.append([[x_iterator, y_iterator, z_iterator],
                                 [[x_iterator+1, y_iterator-1, z_iterator-1], 
                                  [x_iterator-1, y_iterator+1, z_iterator-1], 
                                  [x_iterator-1, y_iterator-1, z_iterator+1],
                                  [x_iterator+1, y_iterator+1, z_iterator+1]
                                  ],
                                  hydrogen_presence])
            hydrogen_atoms.append(hydrogen_presence)

          if x_iterator%4==3 and y_iterator%4==1 and z_iterator%4==3:
            hydrogen_presence=[secrets.randbelow(2),secrets.randbelow(2)]
            oxygen_atoms.append([[x_iterator, y_iterator, z_iterator],
                                 [[x_iterator+1, y_iterator+1, z_iterator-1], 
                                  [x_iterator-1, y_iterator-1, z_iterator-1], 
                                  [x_iterator-1, y_iterator+1, z_iterator+1],
                                  [x_iterator+1, y_iterator-1, z_iterator+1]
                                  ],
                                  hydrogen_presence])
            hydrogen_atoms.append(hydrogen_presence)

          if x_iterator%4==1 and y_iterator%4==3 and z_iterator%4==3:
            hydrogen_presence=[secrets.randbelow(2),secrets.randbelow(2)]
            oxygen_atoms.append([[x_iterator, y_iterator, z_iterator],
                                 [[x_iterator+1, y_iterator+1, z_iterator-1], 
                                  [x_iterator-1, y_iterator-1, z_iterator-1], 
                                  [x_iterator-1, y_iterator+1, z_iterator+1],
                                  [x_iterator+1, y_iterator-1, z_iterator+1]
                                  ],
                                  hydrogen_presence])
            hydrogen_atoms.append(hydrogen_presence)

          if x_iterator%4==3 and y_iterator%4==3 and z_iterator%4==1:
            hydrogen_presence=[secrets.randbelow(2),secrets.randbelow(2)]
            oxygen_atoms.append([[x_iterator, y_iterator, z_iterator],
                                 [[x_iterator+1, y_iterator+1, z_iterator-1], 
                                  [x_iterator-1, y_iterator-1, z_iterator-1], 
                                  [x_iterator-1, y_iterator+1, z_iterator+1],
                                  [x_iterator+1, y_iterator-1, z_iterator+1]
                                  ],
                                  hydrogen_presence])
            hydrogen_atoms.append(hydrogen_presence)


    for i, oxygen_atom in enumerate(oxygen_atoms):
        #print(oxygen_atom)
        for j, neighbor in enumerate(oxygen_atom[1]):
          #print(oxygen_atoms[i][1][j])
          oxygen_atoms[i][1][j]=[x%(net_size*4-1)-1 if x>(net_size*4-1) else x for x in oxygen_atoms[i][1][j]]
          oxygen_atoms[i][1][j]=[x%(net_size*4-1)+1 if x<0 else x for x in oxygen_atoms[i][1][j]]
          #print(oxygen_atoms[i][1][j])

    water_particles=[]
    for i in range(len(oxygen_atoms)):
        water_particles.append(oxygen_atoms[i])

    #rint(len(water_particles), len(oxygen_atoms), len(hydrogen_atoms))
    
    water_particles2=copy.deepcopy(water_particles)
    for i, water_particle in enumerate(water_particles2):
        x_coordinate=water_particle[0][0]
        y_coordinate=water_particle[0][1]
        z_coordinate=water_particle[0][2]

        potential_net[x_coordinate][y_coordinate][z_coordinate]=water_particle
        #print(water_particle)
  

    #print('\nAll water particles:')
    
    for i, water_particle in enumerate(water_particles):
        #water_particles[i].append(calculate_particle_energy(potential_net, water_particle, initial_energy_calculation=True))
        #print('ja', water_particle)

        own_hydrogens_sum=0
        neighbors_hydrogen_sum=0
        x_coordinate_first_neighbor=water_particle[1][0][0]
        y_coordinate_first_neighbor=water_particle[1][0][1]
        z_coordinate_first_neighbor=water_particle[1][0][2]

        x_coordinate_second_neighbor=water_particle[1][1][0]
        y_coordinate_second_neighbor=water_particle[1][1][1]
        z_coordinate_second_neighbor=water_particle[1][1][2]
        #print('somsiad', x_coordinate_first_neighbor, y_coordinate_first_neighbor, z_coordinate_first_neighbor,
        #      potential_net[x_coordinate_first_neighbor][y_coordinate_first_neighbor][z_coordinate_first_neighbor],
        #      x_coordinate_second_neighbor, y_coordinate_second_neighbor, z_coordinate_second_neighbor,
        #      potential_net[x_coordinate_second_neighbor][y_coordinate_second_neighbor][z_coordinate_second_neighbor])
      
        first_neighbor_hydrogens=potential_net[x_coordinate_first_neighbor][y_coordinate_first_neighbor][z_coordinate_first_neighbor][2][0]^1             #należy ustalić czy brać zerowy czy pierwszy wodór od sąsiada
        second_neighbor_hydrogens=potential_net[x_coordinate_second_neighbor][y_coordinate_second_neighbor][z_coordinate_second_neighbor][2][1]^1         #należy ustalić czy brać zerowy czy pierwszy wodór od sąsiada

        #print(water_particle,
        #      first_neighbor_hydrogens,
        #      second_neighbor_hydrogens)
        
        own_hydrogens_sum+=sum(water_particle[2])
        neighbors_hydrogen_sum+=(first_neighbor_hydrogens+second_neighbor_hydrogens)
        energy_penalty=abs(own_hydrogens_sum+neighbors_hydrogen_sum-2)

        water_particles[i].append([first_neighbor_hydrogens, second_neighbor_hydrogens])
        water_particles[i].append(energy_penalty)
    
    
    #for i, water_particle in enumerate(water_particles):
      
    #    print(water_particles[i])
        
    print('\nwhat does it mean?\n', water_particle,
    
          '\nwater particle coordinates:\t', water_particle[0],
          
          '\nfrom who it takes hydrogen:\t', first_neighbor_hydrogens, 'from', water_particle[1][0], '\t', second_neighbor_hydrogens, 'from', water_particle[1][1],
          
          '\nto whom it gives hydrogen:\t', water_particle[2][0], 'to', water_particle[1][2], '\t', water_particle[2][1], 'to', water_particle[1][3],
          
          '\nenergy penalty for incorect number of hydrogen atoms:\t', water_particle[-1],'\n')
          
    print(len(water_particles), len(oxygen_atoms), len(hydrogen_atoms))
    
  
  
    water_particles3=copy.deepcopy(water_particles)
    for i, water_particle in enumerate(water_particles3):
        x_coordinate=water_particle[0][0]
        y_coordinate=water_particle[0][1]
        z_coordinate=water_particle[0][2]

        potential_net[x_coordinate][y_coordinate][z_coordinate]=water_particle
        #print(water_particle)
    
    return water_particles, oxygen_atoms, hydrogen_atoms, potential_net

if __name__== "__main__":
   water_particles, oxygen_atoms, hydrogen_atoms, potential_net=generate_ice_1c_structure()
   for oxygen_atom in oxygen_atoms:
      print(oxygen_atom)
   for hydrogen_atom in hydrogen_atoms:
      print(hydrogen_atom)


