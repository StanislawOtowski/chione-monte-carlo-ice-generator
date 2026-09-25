import time
from multiprocessing import Pool

from Generating.Ic_ice_net_generator import generate_ice_1c_structure
from Plotting.Ic_and_7_ice_net_plotter import plot_ice_1c_or_ice_7_structure
from Plotting.saving_to_CIF_file import save_ice_structure_to_CIF
from Ice_Structure_Changes.change_hydrogen_configuration import make_changes

def generate_and_iterate_net(net_size,
                             ice_type='1c',
                             plot=True,
                             cif_file_generation=True,
                             temperatures=[100,75,60,50,40,30,20,17.5,15,12.5,10,9,8,7,6,5,4,3,2.5,2,1.5,1.25,1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.09,0.08,0.07,0.06,0.05,0.04,0.03,0.02,0.01,0.009,0.008,0.007,0.006,0.005,0.004,0.003,0.002,0.001,0.0009,0.0008,0.0007,0.0006,0.0005,0.0004,0.0003,0.0002,0.0001],
                             iterations_at_one_temperature_at_stable_state=100):
    start_timer=time.time()
    paths_to_CIFs=[]
    #ice net generation
    if ice_type=='1c':
        water_particles, oxygen_atoms, hydrogen_atoms, potential_net=generate_ice_1c_structure(net_size=net_size)
        original_E=0
        for i, water_particle in enumerate(water_particles):
            original_E+=water_particles[i][-1]
        print('making changes\n')
        E, potential_net, water_particles_changed=make_changes(potential_net=potential_net,
                                                    water_particles=water_particles,
                                                    iterations_at_one_temperature_at_stable_state=iterations_at_one_temperature_at_stable_state,
                                                    temperatures=temperatures)
        print('\noriginal energy:\t', original_E)
        final_E=0
        for i, water_particle in enumerate(water_particles_changed):
            final_E+=water_particles_changed[i][-1]
        print('\nfinal energy:\t', final_E)
        print(len(water_particles), 'oxygens')
        end_timer=time.time()
        print("%.2f" % (end_timer-start_timer),'s')
        if cif_file_generation:
            paths_to_CIFs.append(save_ice_structure_to_CIF(water_particles_changed, net_size))
        plot_ice_1c_or_ice_7_structure(water_particles_changed)

    elif ice_type=='7':
        water_particles_net_1, oxygen_atoms_1, hydrogen_atoms_1, potential_net_1=generate_ice_1c_structure(net_size=net_size)
        water_particles_net_2, oxygen_atoms_2, hydrogen_atoms_2, potential_net_2=generate_ice_1c_structure(net_size=net_size)
        original_E=0
        for water_particle in water_particles_net_1:
            original_E+=water_particle[-1]
        for water_particle in water_particles_net_2:
            original_E+=water_particle[-1]

        print('making changes\n')
        E1, potential_net_1, water_particles_changed_net_1=make_changes(potential_net=potential_net_1,
                                                    water_particles=water_particles_net_1,
                                                    iterations_at_one_temperature_at_stable_state=iterations_at_one_temperature_at_stable_state,
                                                    temperatures=temperatures)
        E2, potential_net_2, water_particles_changed_net_2=make_changes(potential_net=potential_net_2,
                                                    water_particles=water_particles_net_2,
                                                    iterations_at_one_temperature_at_stable_state=iterations_at_one_temperature_at_stable_state,
                                                    temperatures=temperatures)
        print('\noriginal energy:\t', original_E)
        final_E=0
        combined_water_particles_changed=[]
        for i, water_particle in enumerate(water_particles_changed_net_1):
            print(water_particles_changed_net_1[i])
            final_E+=water_particles_changed_net_1[i][-1]
            combined_water_particles_changed.append(water_particles_changed_net_1[i])
        for i, water_particle in enumerate(water_particles_changed_net_2):
            water_particles_changed_net_2[i][0][0]+=2
            print(water_particles_changed_net_2[i])
            final_E+=water_particles_changed_net_2[i][-1]
            combined_water_particles_changed.append(water_particles_changed_net_2[i])
        print('\nfinal energy:\t', final_E)
        print(len(water_particles_changed_net_1)*2, 'oxygens')
        end_timer=time.time()
        print("%.2f" % (end_timer-start_timer),'s')
        if cif_file_generation:
           paths_to_CIFs.append(save_ice_structure_to_CIF(combined_water_particles_changed, net_size, ice_type='7'))
        plot_ice_1c_or_ice_7_structure(combined_water_particles_changed)
    else:
        print('Ice type unknown. Try 1c or 7.')

if __name__=='__main__':
    print('main')
    paths_to_CIFs=[]
    for i in range(1):          #tutaj podajemy ile razy chcemy użyć zestawu rdzeni
        with Pool(4) as p:     #tutaj podajemy liczbę rdzeni, których chcemy użyć
            paths_to_CIFs.append(p.map(generate_and_iterate_net, [4])) #tutaj podajemy rozmiary siatek, które chcemy generować równolegle na rdzeniach
