def calculate_total_energy(water_particles):
    E=0
    for water_particle in water_particles:
        E+=water_particle[-1]
    print('total energy\t', E)
    return E

def calculate_local_energy_delta(potential_net, water_particle, new_hydrogen_position, new_hydrogen_value):
    delta_E=0
    #print('Ą', water_particle)
    if new_hydrogen_position==0:
        neighbor=potential_net[water_particle[1][2][0]][water_particle[1][2][1]][water_particle[1][2][2]]
        #rint('A neighbor coordinates', water_particle[1][2])
        #rint('B calculate_total_energy function, neighbor var:', neighbor)
        #calculate delta_E in neighbor
        delta_E+=abs(sum(neighbor[-3])+(new_hydrogen_value^1)+neighbor[-2][1]-2)-neighbor[-1]                                         #abs(sum(neighbor[-3])+sum(neighbor[-2])-2)
        #calculate delta_E in water_particle
        delta_E+=abs(new_hydrogen_value+water_particle[-3][1]+sum(water_particle[-2])-2)-water_particle[-1]                         #abs(sum(water_particle[-3])+sum(water_particle[-2])-2)
    
    elif new_hydrogen_position==1:
        neighbor=potential_net[water_particle[1][3][0]][water_particle[1][3][1]][water_particle[1][3][2]]
        
        #print('C neighbor coordinates', water_particle[1][3])
        #print('D calculate_total_energy function, neighbor var:', neighbor)

        #calculate delta_E in neighbor
        delta_E+=abs(sum(neighbor[-3])+neighbor[-2][0]+(new_hydrogen_value^1)-2)-neighbor[-1]                                         #abs(sum(neighbor[-3])+sum(neighbor[-2])-2)
        #calculate delta_E in water_particle
        delta_E+=abs(water_particle[-3][0]+new_hydrogen_value+sum(water_particle[-2])-2)-water_particle[-1]                         #abs(sum(water_particle[-3])+sum(water_particle[-2])-2)

    return delta_E
