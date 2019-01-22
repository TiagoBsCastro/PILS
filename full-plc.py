#!/usr/bin/env python
from __future__ import print_function
import os,sys
import numpy as np
from string import Template
import subprocess
from astropy.cosmology import FlatLambdaCDM, z_at_value
from astropy.io import fits
from glob import glob
import astropy.units as u

################### Reading Args. #######################

if len(sys.argv) == 7:

        n_plc=int(sys.argv[1])
        box_size=float(sys.argv[2])
        npixels=int(sys.argv[3])                 # Number of Pixels for Lensing Maps
        do_pinocchio = int(sys.argv[4])          # Run Pinocchio
        do_wlmoka = int(sys.argv[5])             # Run wlmoka
        do_gamma = int(sys.argv[6])              # Run MapGetGamma

else:

        print( "Number of PLC not specified!")
        sys.exit(-1)

if do_gamma and do_wlmoka:

        print("Gamma Maps can NOT be created at the same time as Convergence Maps!")
        raise(InputError)


mpi_log = open("log_{0:02d}".format(n_plc),"w")

######################### Input ########################

seed = -9999                  # Seed for Pinocchio ICs and numpy
Om0 = 0.301                   # Omega Matter
hubble = 0.682                # Little h
zlast_PLC = 5.0               # PLC's depth (should be less or equal to highest zsource)
zsource = [1.0,2.0,3.0,5.0]   # Source redshift
fov = 1.0                     # Angular aperture of the FOV
resolution = 'hr'             # Resolution for Pinocchios Run
                              # mr:64**3 hr:128**3 uhr:256**3

#########################################################

nodes = 1                                                    # Number of nodes to be used
ncores ={'mr':2, 'hr':2}[resolution]                         # Total number of cores to be used
np.random.seed(np.abs(seed)*n_plc)                           # Starting numpy random values

################ Instatializing Cosmo ###################
# !!!!!! It should match the PLC cosmo values !!!!!!!!! #
# !!!!!!        I will work with Mpc/h        !!!!!!!!! #

cosmo = FlatLambdaCDM(H0=100., Om0=Om0)

############### Pinocchio Params File and exe ################

template_pin=Template(open("pinocchio_slave.ini").read())
pin_exe="/home/tcastro/Pinocchio/bin/pinocchio.x"
npart = {'mr':64,'hr':128,'uhr':256}[resolution]

############### WeakLMoka Params File and exe ################

template_wl=Template(open("pinmoka_slave.ini").read())
wl_exe="/home/tcastro/PinMoka-light/bin/PinMoka.x"

############### WeakLMoka Params File and exe ################

gg_exe="/home/tcastro/MapGetGamma/bin/MapGetGamma.x"

##############################################################

############### Randomization Variables #################
# Vertice position of the Light-Cone
rand_pos=[]
# Mirror random variables x->+/-x, y->+/-y
rand_mirror=[]
invert_axis=[]
# Ray-tracing direction
rand_fov=[]
#########################################################

zini=0
zini_tab=[]
zlast=0
zlast_tab=[]
nrep=0
plc_tab=[]
dir="Alice_PLC/"+resolution+"_"+str(n_plc).zfill(2)+"/"

while zlast<zlast_PLC:

	nrep+=1
	zini = z_at_value(cosmo.comoving_distance,nrep * box_size * u.Mpc)

	# Vertice position of the Light-Cone
	rand_pos.append(np.random.rand(3))
	# Mirror random variables x->+/-x, y->+/-y...
	rand_mirror.append(2*np.random.randint(0,2,3)-1)
	# Ray-tracing direction
	rand_fov.append(np.random.randint(0,3))

	xc,yc,zc=rand_pos[-1]*box_size
	xd,yd,zd=[ (i==rand_fov[-1])*rand_mirror[-1][i] for i in range(3)]
	invert_axis.append(np.array([rand_mirror[-1][i] for i in range(3) if rand_fov[-1]!=i]))

	if not do_pinocchio:

		if os.path.isfile(dir+"pinocchio.{0:.4f}.alice.plc.out".format(zini)):

			plc_tab.append("pinocchio.{0:.4f}.alice.plc.out".format(zini))
			zlast_tab.append(zlast)
	       		zini_tab.append(zini)
		        zlast=zini
			continue

		else:

			print("Pinocchio PLC for z={} was not created and do_pinocchio is set to {}".format(zini,do_pinocchio))
			raise(ImportError)


	print("---------------------- Zini = {0:.4f} ------------------------".format(zini),file=mpi_log)
	if rand_mirror[-1][rand_fov[-1]]>0:
		print( "The PLC will be contructed through the",{0:"x",1:"y",2:"z"}[rand_fov[-1]],"axis",file=mpi_log)
	else:
		print("The PLC will be contructed through the",{0:"-x",1:"-y",2:"-z"}[rand_fov[-1]],"axis",file=mpi_log)

	print("The other axis will be inverted acording to {} {}".format(invert_axis[-1][0],invert_axis[-1][1])\
													,file=mpi_log)
	print("The PLC will be centered in the following cordinates {} {} {}".format(xc,yc,zc),file=mpi_log)

	os.system("mkdir -p "+dir)
	os.system("mkdir -p "+dir+"/pinocchio")

	output_name="outputs_{0:02d}_{1:.4f}".format(n_plc,zini)
        outputs_file=open(dir+output_name,"w")
        outputs_file.write(str(zini))
        outputs_file.close()

	param_name="parameter_file_"+resolution+"_{0:02d}_{1:.4f}".format(n_plc,zini)
	param_file=open(dir+param_name,"w")
	params=template_pin.substitute(Outputs=output_name,seed=seed,Om0=Om0,OL0=1-Om0,h=hubble,\
					BoxSize=box_size,zini=zini,zlast=zlast,fov=fov,Xc=xc,Yc=yc,Zc=zc,Xd=xd,Yd=yd,Zd=zd,npart=npart)
	param_file.write(params)
	param_file.close()

	log=open(dir+"log_{0:02d}_{1:.4f}.txt".format(n_plc,zini),"w")
	err=open(dir+"err_{0:02d}_{1:.4f}.txt".format(n_plc,zini),"w")
	returned_value=subprocess.call(["mpirun","-n", "{}".format(ncores),pin_exe,param_name], cwd=dir, stdout=log, stderr=err)
	log.close()
	err.close()
	subprocess.call(["mv","pinocchio.alice.plc.out","pinocchio.{0:.4f}.alice.plc.out".format(zini)],cwd=dir)
	subprocess.call(["mv","pinocchio.alice.cosmology.out","pinocchio.{0:.4f}.alice.cosmology.out".format(zini)],cwd=dir)
	subprocess.call(["mv","pinocchio.alice.geometry.out","pinocchio.{0:.4f}.alice.geometry.out".format(zini)],cwd=dir)

	plc_tab.append("pinocchio.{0:.4f}.alice.plc.out".format(zini))
	zlast_tab.append(zlast)
	zini_tab.append(zini)
	zlast=zini

	if returned_value != 0:

		print("Pinocchio returned".format(returned_value),file=mpi_log)
		print( "Check the log and err files in "+ dir,file=mpi_log)
		print( "I'll Stop here",file=mpi_log)
		quit()

	else:

		print( "Pinocchio returned {}".format(returned_value),file=mpi_log)
                print( "\n",file=mpi_log)


	subprocess.call("mv *catalog* *mf* log* err* output* param* *cosmology* *geometry* pinocchio/",shell=True,cwd=dir)


print(  "\n############## The Full-PLC was created ###############",file=mpi_log)
print(    "#              Starting the WLMoka phase              #",file=mpi_log)
print(    "#######################################################\n",file=mpi_log)

if do_wlmoka:

	for zs in zsource:

		print( "-------------------- Zs = {0:.4f} --------------------\n".format(zs),file=mpi_log)
		nrep = np.sum( np.asarray(zlast_tab) < zs )

		for i in range(0,nrep,ncores):

			log=[]
			err=[]
			processes=[]

			if i+ncores>nrep:
				nprocesses=nrep-i
			else:
				nprocesses=ncores

			for p in range(nprocesses):

				param_name="weaklensingMOKA_{0:02d}_{1:.4f}".format(n_plc,zini_tab[i+p])
		                param_file=open(dir+param_name,"w")
		                params=template_wl.substitute(zsource=zs,PLC=plc_tab[i+p],fov=fov,npixels=npixels,Om0=Om0,OL0=1-Om0,h=hubble)
		                param_file.write(params)
		                param_file.close()

				print( "Open subprocess {0:d} for running WLMoka on Zini = {1:.4f}".format(p,zini_tab[i+p]),file=mpi_log)

				log.append(open(dir+"log_{0:02d}_{1:.4f}.txt".format(n_plc,zini_tab[i+p]),"w"))
		                err.append(open(dir+"err_{0:02d}_{1:.4f}.txt".format(n_plc,zini_tab[i+p]),"w"))
		                processes.append(subprocess.Popen(["mpirun","-n","1",wl_exe,param_name,"0"], cwd=dir, stdout=log[-1], stderr=err[-1]))

			for p in range(nprocesses):

				print( "Waiting subprocess {0:d} to finish".format(p),file=mpi_log)
				processes[p].wait()
				log[p].close()
	        	        err[p].close()
				print("Subprocess {0:d} is done".format(p),file=mpi_log)
				print("I'll check if I have to flip axis",file=mpi_log)

				if (invert_axis[i+p]==[1,1]).all():

					continue

				else:

					kappa_file=dir+"kappa_{0:.4f}.alice.plc_0.fits".format(zini_tab[i+p])
					kappa,header=fits.getdata(kappa_file,header=True)

					if (invert_axis[i+p]==[-1,1]).all():

						kappa=np.flip(kappa,0)

					elif (invert_axis[i+p]==[1,-1]).all():

						kappa=np.flip(kappa,1)

					elif (invert_axis[i+p]==[-1,-1]).all():

						kappa=np.flip(np.flip(kappa,1),0)

					fits.writeto(kappa_file,data=kappa,header=header,overwrite=True)

				print( "Done\n",file=mpi_log)

			print( "Moving EVERY Lens file to: {}".format(dir+str(zs)+"_"+str(npixels)),file=mpi_log)
			os.system("mkdir -p {}".format(dir+str(zs)+"_"+str(npixels)))
			subprocess.call("mv *kappa* weak* log* err*  {}".format(str(zs)+"_"+str(npixels)),shell=True, cwd=dir)

if do_gamma:

	print(  "\n#######################################################",file=mpi_log)
	print(    "#              Starting the Shear Maps                #",file=mpi_log)
	print(    "#######################################################\n",file=mpi_log)

	dir="Alice_PLC/"+resolution+"_*/"

	for zs in zsource:

	        kappa_dir = sorted(glob(dir+str(zs)+"_"+str(npixels)+"/"))
		print(kappa_dir)
		nsims = len(kappa_dir)
		kappa_files = sorted([ file_name.split("/")[-1] for file_name in glob(kappa_dir[0]+"/kappa_*.fits")])
		nplanes = len(kappa_files)
		km = np.zeros( nplanes )
	        kappa = np.zeros([npixels,npixels])

		print("\n  Computing the Mean density over all the planes \n",file=mpi_log)
	        for i,n in enumerate(kappa_files):

			for d in kappa_dir:

		                print("Reading Convergence map {}".format(d+n),file=mpi_log)
		                kappa += fits.getdata(d+n)

			km[i] = kappa.mean()/nsims
			kappa = np.zeros([npixels,npixels])

		print("\n  Computing the effective convergence for sim = {}\n".format(n_plc),file=mpi_log)
		d = dir.replace("*",str(n_plc).zfill(2))+str(zs)+"_"+str(npixels)+"/"
		for i,n in enumerate(kappa_files):

                        print("Reading Convergence map {}".format(d+n),file=mpi_log)
                        kappa += fits.getdata(d+n) - km[i]
			header = fits.getheader(d+n)

		header['HIERARCH PHYSICALSIZE'] = float(fov)
		header['HIERARCH REDSHIFTSOURCE'] = float(zs)
	        fits.writeto(open(d+"/kappa.fits","wb"), kappa, header=header,overwrite=True)
		kappa = np.zeros([npixels,npixels])

        	subprocess.call(["mpirun","-n","1",gg_exe,"kappa.fits"], cwd=d)

quit()
