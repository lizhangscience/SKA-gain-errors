from definitions import *
from pyrap.tables import table
import os
import time

if __name__ == "__main__":

    t1 = time.time()

    #Setup file settings
    setup_file = 'setup_SKA_core.ini'
    vis_file = 'SKA_core' #visibility filename (prefix) to write visibilities to

    #Telescope configuration
    num_antennas = 256
    num_stations = 224

    #Telescope directory
    telescope_directory = 'SKA_core.tm'     
    #Skymodels directory
    skymodels_directory = 'skymodels_mf'
    #Images directory
    images_directory = 'results/images_mf'

    #FG type
    fg_type = 'diffuse'

    #Output directories
    plots_directory = 'results/plots_mf/%s'%fg_type
    data_directory = 'results/data_mf/%s'%fg_type

    #Imaging parameters
    imgsize = 300 #image size in pixels
    cellsize_arcmin = 2 #pixel spacing   
    fov = imgsize*cellsize_arcmin/60.
 
    #Add frequency and noise files in setup file
    set_settings (setup_file, 'interferometer/noise/enable', 'true')
   
    #Add noise to telescope
    frequencies = np.linspace(115e6, 125e6, 11)
    sefd_stn_0 = 1.25
    frequency_0 = frequencies[0]

    #Tuples of gain errors (G_std, phi_std)
    rms_gain_errors = [0.01, 0.1, 1]

    #quantities = ['fg_%s'%fg_type, 'noise', 'eor', 'fg_p_eor_p_ge_m_fg', 'fg_p_eor_p_ge_m_fg_m_eor']
    quantities = ['fg_%s'%fg_type, 'fg_p_eor_p_ge_m_fg', 'fg_p_eor_p_ge_m_fg_m_eor']

    #Clean up directories
    for qnt in quantities:
        os.system('rm %s/%s/*'%(images_directory, qnt))

    #Observation time parameters
    start_time = '2015-01-01T06:00:00' #Start time in UTC, in yyyy-M-dTh:m:s format
    timestep = 5 #Time step in minutes
    duration = 0.5 #Total observation time in hours
    #duration = 1.5
    duration = duration*60 #Convert duration into minutes

    #Timescale of correlated gains 
    timescale = 1. #timescale in hours
    timescale = timescale*60 #Convert timescale to minutes

    #Simulate visibilities and create images for FG, noise and EoR
    for qnt in ['fg', 'noise', 'eor']:
 
        for frequency in frequencies:

            sefd_stn = sefd_stn_0*(frequency_0/frequency)**0.55
            frequency=[frequency]

            if qnt=='fg': rms = [0.]
            if qnt=='noise': rms = [sefd_stn]
            if qnt=='eor': rms = [0.]

            #Skymodel filenames for FG and EoR
            #These skymodels cover a 10x10 FOV and the pixel spacing is 1.17 arcmin 

            if fg_type == 'diffuse':
                skymodel_fg_filename = '%s/fg_%s/fg_image_jy_per_px_%s_Hz.osm'%(skymodels_directory, fg_type, str(int(frequency[0])))
            if fg_type == 'point':
                skymodel_fg_filename = '%s/fg_%s/ScubeSkymodel1mJy.osm'%(skymodels_directory, fg_type)

            skymodel_eor_filename = '%s/eor/eor_image_jy_per_px_%s_Hz.osm'%(skymodels_directory, str(int(frequency[0])))
            skymodel_empty_filename = '%s/empty_sky/sky_empty_1.17arcmin.osm'%skymodels_directory

            skymodel_filename = {'fg':skymodel_fg_filename, 'noise':skymodel_empty_filename, 'eor':skymodel_eor_filename}   

            add_noise(telescope_directory, frequency, rms)

            set_settings (setup_file, 'sky/oskar_sky_model/file', skymodel_filename['%s'%qnt])   
            set_settings (setup_file, 'observation/start_frequency_hz', frequency[0]) 
            make_visibilities_ms(setup=setup_file, telescope_directory=telescope_directory, start_time=start_time, time_interval=timestep, duration=duration, output_ms='%s.ms'%vis_file)

            img_subdirectory = {'fg':'fg_%s'%fg_type, 'noise':'noise', 'eor':'eor'}

            img_filename = '%s/%s/%s_%s_Hz'%(images_directory, img_subdirectory[qnt], qnt, str(int(frequency[0])))
            make_image_wsclean(ms='%s.ms'%vis_file, imgsize=imgsize, cellsize_arcmin=cellsize_arcmin, img_filename=img_filename)

        os.system('python /home/users/modhurita/ps_eor/simple_ps_save.py %s/%s/*'%(images_directory, img_subdirectory[qnt]))
        os.system('cp power_spectra_spatial.eps %s/ps_%s_spatial.eps'%(plots_directory, qnt))
        os.system('cp power_spectra_2d.eps %s/ps_%s_2d.eps'%(plots_directory, qnt))
        os.system('cp power_spectra_3d.eps %s/ps_%s_3d.eps'%(plots_directory, qnt))

        os.system('cp power_spectra_spatial.txt %s/ps_%s_spatial.txt'%(data_directory, qnt))
        os.system('cp power_spectra_2d.txt %s/ps_%s_2d.txt'%(data_directory, qnt))
        os.system('cp power_spectra_3d.txt %s/ps_%s_3d.txt'%(data_directory, qnt))

    #Simulate visibilities and create images for FG + noise + EoR
 
    for gain_error_type in ['correlated', 'uncorrelated']:

        for rms_gain in rms_gain_errors:
            if gain_error_type == 'uncorrelated':
                gain_errors = generate_uncorrelated_gain_errors(rms_gain=rms_gain, timestep=timestep, duration=duration, num_antennas=num_antennas, num_stations=num_stations)
            if gain_error_type == 'correlated':
                gain_errors = generate_correlated_gain_errors(rms_gain=rms_gain, timescale=timescale, timestep=timestep, duration=duration, num_antennas=num_antennas, num_stations=num_stations)
            for frequency in frequencies:

                sefd_stn = sefd_stn_0*(frequency_0/frequency)**0.55
                frequency=[frequency]

                #Skymodel filename for FG+EoR
                #These skymodels cover a 10x10 FOV and the pixel spacing is 1.17 arcmin
                if fg_type == 'diffuse':
                    skymodel_fg_filename = '%s/fg_%s/fg_image_jy_per_px_%s_Hz.osm'%(skymodels_directory, fg_type, str(int(frequency[0])))
                if fg_type == 'point':
                    skymodel_fg_filename = '%s/fg_%s/ScubeSkymodel1mJy.osm'%(skymodels_directory, fg_type)
               
                skymodel_eor_filename = '%s/eor/eor_image_jy_per_px_%s_Hz.osm'%(skymodels_directory, str(int(frequency[0])))

                os.system('cat %s %s > fg_plus_eor.osm'%(skymodel_fg_filename, skymodel_eor_filename))

                skymodel_fg_plus_eor_filename = 'fg_plus_eor.osm'

                rms = [0.]
                add_noise(telescope_directory, frequency, rms)

                set_settings (setup_file, 'sky/oskar_sky_model/file', skymodel_fg_plus_eor_filename)   
                set_settings (setup_file, 'observation/start_frequency_hz', frequency[0])
 
                make_visibilities_ms(setup=setup_file, telescope_directory=telescope_directory, gain_errors=gain_errors, start_time=start_time, time_interval=timestep, duration=duration, output_ms='%s.ms'%vis_file)
                img_filename = '%s/fg_p_eor_p_ge/fg_p_eor_p_ge_%s_Hz'%(images_directory, str(int(frequency[0])))
                make_image_wsclean(ms='%s.ms'%vis_file, imgsize=imgsize, cellsize_arcmin=cellsize_arcmin, img_filename=img_filename)

                #(FG + EoR + gain errors) - FG image
                subtract_images('%s/fg_p_eor_p_ge/fg_p_eor_p_ge_%s_Hz.fits'%(images_directory, str(int(frequency[0]))), '%s/fg_%s/fg_%s_Hz.fits'%(images_directory, fg_type, str(int(frequency[0]))), '%s/fg_p_eor_p_ge_m_fg/fg_p_eor_p_ge_m_fg_%s_Hz.fits'%(images_directory, str(int(frequency[0]))))

                #(FG + EoR + gain errors) - FG image - EoR image
                subtract_images('%s/fg_p_eor_p_ge_m_fg/fg_p_eor_p_ge_m_fg_%s_Hz.fits'%(images_directory, str(int(frequency[0]))), '%s/eor/eor_%s_Hz.fits'%(images_directory, str(int(frequency[0]))), '%s/fg_p_eor_p_ge_m_fg_m_eor/fg_p_eor_p_ge_m_fg_m_eor_%s_Hz.fits'%(images_directory, str(int(frequency[0]))))

            for qnt in ['fg_p_eor_p_ge_m_fg_m_eor']:
                os.system('python /home/users/modhurita/ps_eor/simple_ps_save.py images_mf/%s/*'%qnt)

                os.system('cp power_spectra_spatial.eps %s/ps_%s_%s_%s_spatial_ge_rms_%s.eps'%(plots_directory, qnt, fg_type, gain_error_type, str(rms_gain)))
                os.system('cp power_spectra_2d.eps %s/ps_%s_%s_%s_2d_ge_rms_%s.eps'%(plots_directory, qnt, fg_type, gain_error_type, str(rms_gain)))
                os.system('cp power_spectra_3d.eps %s/ps_%s_%s_%s_3d_ge_rms_%s.eps'%(plots_directory, qnt, fg_type, gain_error_type, str(rms_gain)))

                os.system('cp power_spectra_spatial.txt %s/ps_%s_%s_%s_spatial_ge_rms_%s.txt'%(data_directory, qnt, fg_type, gain_error_type, str(rms_gain)))
                os.system('cp power_spectra_2d.txt %s/ps_%s_%s_%s_2d_ge_rms_%s.txt'%(data_directory, qnt, fg_type, gain_error_type, str(rms_gain)))
                os.system('cp power_spectra_3d.txt %s/ps_%s_%s_%s_3d_ge_rms_%s.txt'%(data_directory, qnt, fg_type, gain_error_type, str(rms_gain)))

    t2 = time.time()
    print 'Time taken: %f min'%((t2-t1)/60.)
