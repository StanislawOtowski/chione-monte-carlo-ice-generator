import numpy as np
from Analyzys.energy_calculations import *
import copy
import time
import secrets

randomizer=secrets.SystemRandom()

def change_acceptance(E, kT, potential_net, water_particle, new_hydrogen_position, new_hydrogen_value):
    delta_E=calculate_local_energy_delta(potential_net, water_particle, new_hydrogen_position, new_hydrogen_value)
    #print('delta_E\t\t:', delta_E)
    #print('np.exp(-delta_E/kT)\t:', np.exp(-delta_E/kT))
    if delta_E<=0:
        return E+delta_E
    elif np.exp(-delta_E/kT) > randomizer.random():
        return E+delta_E
    else:
        #print('change_rejected')
        return 'change_rejected'
    
def change_hydrogen_configuration(E,
                                  kT,
                                  potential_net,
                                  water_particles):
    
    hydrogen_atoms_count=len(water_particles)*2
    change_acceptance_rate=[0,0]
    for i in range(hydrogen_atoms_count):
        change_acceptance_rate[1]+=1
        #print('\nNew change try:')
        picked_number=secrets.randbelow(len(water_particles))
        new_hydrogen_position=secrets.randbelow(2)
        new_hydrogen_value=secrets.randbelow(2)
        if water_particles[picked_number][-3][new_hydrogen_position]==new_hydrogen_value:
            #print('No effect.')
            pass
        else:
            new_E=change_acceptance(E=E,
                                kT=kT,
                                potential_net=potential_net,
                                water_particle=water_particles[picked_number],
                                new_hydrogen_position=new_hydrogen_position,
                                new_hydrogen_value=new_hydrogen_value)
            if new_E != 'change_rejected':
                #print('changed:\t', water_particles[picked_number], new_hydrogen_position, new_hydrogen_value)
                #changing particle in water_particles
                water_particles[picked_number][-3][new_hydrogen_position] = new_hydrogen_value
                #updating particle energy in water_particles
                water_particles[picked_number][-1]=abs(sum(water_particles[picked_number][-3])+sum(water_particles[picked_number][-2])-2)
                #changing particle in potential_net
                water_particle_coordinates=copy.copy(water_particles[picked_number][0])
                x_coordinate, y_coordinate, z_coordinate = water_particle_coordinates
                potential_net[x_coordinate][y_coordinate][z_coordinate][-3][new_hydrogen_position] = new_hydrogen_value
                #updating particle energy in potential_net
                potential_net[x_coordinate][y_coordinate][z_coordinate][-1]=abs(sum(potential_net[x_coordinate][y_coordinate][z_coordinate][-3])+sum(potential_net[x_coordinate][y_coordinate][z_coordinate][-2])-2)
                
                if new_hydrogen_position==0:
                    neighbor=potential_net[water_particles[picked_number][1][2][0]][water_particles[picked_number][1][2][1]][water_particles[picked_number][1][2][2]]
                    #changing neighbor particle in potential_net
                    neighbor[-2][new_hydrogen_position]=new_hydrogen_value^1
                    #updating neighbor particle energy in potential_net
                    neighbor[-1]=abs(sum(neighbor[-3])+sum(neighbor[-2])-2)
                    
                    #poniższy kod trzeba zoptymnalizować
                    for wp in water_particles:
                        if wp[1][0]==water_particle_coordinates:
                            #changing neighbor particle in water_particles
                            wp[-2][new_hydrogen_position]=new_hydrogen_value^1
                            #updating neighbor particle energy in water_particles
                            wp[-1]=abs(sum(wp[-3])+sum(wp[-2])-2)
                            #print('wp\t', wp)

                elif new_hydrogen_position==1:
                    neighbor=potential_net[water_particles[picked_number][1][3][0]][water_particles[picked_number][1][3][1]][water_particles[picked_number][1][3][2]]
                    #changing neighbor particle in potential_net
                    neighbor[-2][new_hydrogen_position]=new_hydrogen_value^1
                    #updating neighbor particle energy in potential_net
                    neighbor[-1]=abs(sum(neighbor[-3])+sum(neighbor[-2])-2)
                    
                    for wp in water_particles:
                        if wp[1][1]==water_particle_coordinates:
                            #changing neighbor particle in water_particles
                            wp[-2][new_hydrogen_position]=new_hydrogen_value^1
                            #updating neighbor particle energy in water_particles
                            wp[-1]=abs(sum(wp[-3])+sum(wp[-2])-2)
                            #print('wp\t', wp)

                change_acceptance_rate[0]+=1
                E=new_E
                #print('total energy:\t', E)
                #print('Accepted changes:\t', change_acceptance_rate[0], '/', change_acceptance_rate[1])
                #for i, water_particle in enumerate(water_particles):
                #    print(water_particles[i])
    #print('New energy:\t', E)
    
    #for i, water_particle in enumerate(water_particles):
    #    print(water_particles[i])
    return E, potential_net, water_particles#, change_acceptance_rate[1]
        

def make_changes(potential_net,
                 water_particles,
                 iterations_at_one_temperature_to_get_stable_state=False,
                 iterations_at_one_temperature_at_stable_state=10,
                 temperatures=[100]):
    if not iterations_at_one_temperature_to_get_stable_state:
        iterations_at_one_temperature_to_get_stable_state=int(iterations_at_one_temperature_at_stable_state/2)
    E=calculate_total_energy(water_particles)
    samples_count=iterations_at_one_temperature_at_stable_state
    water_particles_count=len(water_particles)
    for kT in temperatures:
        start_kT_timer=time.time()
        sum_E=0
        sum_E_2=0
        sum_E_3=0
        sum_E_4=0
        for i in range(iterations_at_one_temperature_to_get_stable_state):
            E, potential_net, water_particles=change_hydrogen_configuration(E=E,
                                                                            kT=kT,
                                                                            potential_net=potential_net,
                                                                            water_particles=water_particles)
        #stable state achieved
        for i in range(samples_count):
            E, potential_net, water_particles=change_hydrogen_configuration(E=E,
                                                                            kT=kT,
                                                                            potential_net=potential_net,
                                                                            water_particles=water_particles)
            #zbieramy momenty energii całkowitej układu
            sum_E  +=E
            sum_E_2+=E**2
            sum_E_3+=E**3
            sum_E_4+=E**4
        #średnie po próbkach: <E>, <E^2>, <E^3>, <E^4>
        mean_E  =sum_E/samples_count
        mean_E_2=sum_E_2/samples_count
        mean_E_3=sum_E_3/samples_count
        mean_E_4=sum_E_4/samples_count

        #twierdzenie fluktuacyjno-dyssypacyjne: C = Var(E)/(kT)^2, na cząsteczkę dzielimy przez N
        var_E=mean_E_2-mean_E**2
        cV=var_E/(water_particles_count*kT**2)
        
        end_kT_timer=time.time()
        print('\n=====\tkT:\t',   kT,
              '\tEnergy:\t',    E,
              '\t<e>:\t',   "%.3e" % (mean_E  /water_particles_count),      #średnia kara energetyczna na cząsteczkę
              '\t<e^2>:\t', "%.3e" % (mean_E_2/water_particles_count**2),
              '\t<e^3>:\t', "%.3e" % (mean_E_3/water_particles_count**3),
              '\t<e^4>:\t', "%.3e" % (mean_E_4/water_particles_count**4),
              '\tcV:\t',    "%.3e" % cV,
              '\tcalculation time:\t', "%.2f" % (end_kT_timer-start_kT_timer),'s',
              '\t=====')
            #for i, water_particle in enumerate(water_particles):
            #    print(water_particles[i])
    print('\niterations per kT:\t', iterations_at_one_temperature_at_stable_state)
    return E, potential_net, water_particles

