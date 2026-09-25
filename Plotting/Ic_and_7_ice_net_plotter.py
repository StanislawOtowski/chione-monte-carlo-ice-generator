import matplotlib.pyplot as plt

def plot_ice_1c_or_ice_7_structure(water_particles):
    print('rysujemy')
    fig = plt.figure(figsize=(15,15))
    ax = fig.add_subplot(projection='3d')
    #plt.rcParams["figure.figsize"] = (10,10)#plt.rcParamsDefault["figure.figsize"]
    oxygen_atoms=water_particles
    hydrogen_atoms=[]
    for i, water_particle in enumerate(water_particles):
        hydrogen_atoms.append(water_particle[-3])
    #print(hydrogen_atoms[:10])
    for i, oxygen_atom in enumerate(oxygen_atoms):
        oxygen_atom_x=oxygen_atom[0][0]
        oxygen_atom_y=oxygen_atom[0][1]
        oxygen_atom_z=oxygen_atom[0][2]

        #print('old i:\t', i)
        #i=i%(len(oxygen_atoms)-1)
        #print('new i:\t', i)
        
        H_O_d_1=0.25        #x,y,z coord lenght of OH bond in ice 1c
        H_O_d_2=1-H_O_d_1   #x,y,z coord lenght of H2O to H2O particle bond in ice 1c
        ax.scatter(oxygen_atom_x, oxygen_atom_y, oxygen_atom_z, marker='o', color='red', s=200)
        alpha_factor_in_net=0.5
        if oxygen_atom[0][2]%2==0:
            plt.plot([oxygen_atom_x, oxygen_atom_x+1], [oxygen_atom_y, oxygen_atom_y+1],  [oxygen_atom_z, oxygen_atom_z+1], '--', color='grey', alpha=alpha_factor_in_net)
            plt.plot([oxygen_atom_x, oxygen_atom_x-1], [oxygen_atom_y, oxygen_atom_y-1],  [oxygen_atom_z, oxygen_atom_z+1], '--', color='grey', alpha=alpha_factor_in_net)
            if hydrogen_atoms[i][0]:
                ax.scatter(oxygen_atom_x-H_O_d_1, oxygen_atom_y-H_O_d_1, oxygen_atom_z+H_O_d_1, marker='o', color='blue', s=100)
            else:
                ax.scatter(oxygen_atom_x-H_O_d_2, oxygen_atom_y-H_O_d_2, oxygen_atom_z+H_O_d_2, marker='o', color='blue', s=100)
            if hydrogen_atoms[i][1]:
                #hydrogen_atom_x=oxygen_atom_x-H_O_distance
                ax.scatter(oxygen_atom_x+H_O_d_1, oxygen_atom_y+H_O_d_1, oxygen_atom_z+H_O_d_1, marker='o', color='blue', s=100)
            else:
                ax.scatter(oxygen_atom_x+H_O_d_2, oxygen_atom_y+H_O_d_2, oxygen_atom_z+H_O_d_2, marker='o', color='blue', s=100)
        if oxygen_atom[0][2]%2==1:
            plt.plot([oxygen_atom_x, oxygen_atom_x+1], [oxygen_atom_y, oxygen_atom_y-1],  [oxygen_atom_z, oxygen_atom_z+1], '--', color='grey', alpha=alpha_factor_in_net)
            plt.plot([oxygen_atom_x, oxygen_atom_x-1], [oxygen_atom_y, oxygen_atom_y+1],  [oxygen_atom_z, oxygen_atom_z+1], '--', color='grey', alpha=alpha_factor_in_net)
            if hydrogen_atoms[i][0]:
                ax.scatter(oxygen_atom_x-H_O_d_1, oxygen_atom_y+H_O_d_1, oxygen_atom_z+H_O_d_1, marker='o', color='blue', s=100)
            else:
                ax.scatter(oxygen_atom_x-H_O_d_2, oxygen_atom_y+H_O_d_2, oxygen_atom_z+H_O_d_2, marker='o', color='blue', s=100)
        
            if hydrogen_atoms[i][1]:
                ax.scatter(oxygen_atom_x+H_O_d_1, oxygen_atom_y-H_O_d_1, oxygen_atom_z+H_O_d_1, marker='o', color='blue', s=100)
            else:
                ax.scatter(oxygen_atom_x+H_O_d_2, oxygen_atom_y-H_O_d_2, oxygen_atom_z+H_O_d_2, marker='o', color='blue', s=100)

    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')

    plt.title('Ice Ic')

    plt.show()
