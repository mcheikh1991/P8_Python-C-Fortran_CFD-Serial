# Start of the makefile
# Defining variables
objects = weno-sm.f90 weno-rc.f90 weno-wg.f90 weno.f90 compute_faces.f90
exe     = weno-sm.o   weno-rc.o   weno-wg.o   weno.o  compute_faces.o
f90comp = gfortran
switch  = -O2
# Makefile
# weno-sm.o weno-sm.mod: weno-sm.f90
# 	$(f90comp) -c $(switch) weno-sm.f90
# weno-rc.o weno-rc.mod: weno-rc.f90
# 	$(f90comp) -c $(switch) weno-rc.f90
# weno-wg.o weno-wg.mod: weno-wg.f90
# 	$(f90comp) -c $(switch) weno-wg.f90
$(exe): $(objects)
	$(f90comp) $(objects) -c $(switch)
# Cleaning everything
clean:
	rm $(exe)