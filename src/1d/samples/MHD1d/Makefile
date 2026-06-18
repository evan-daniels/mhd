APPNAME   = MHD_1d

CPROOT    = ../../../..
SAMPLEDIR = $(CPROOT)/samples
APPDIR    = $(SAMPLEDIR)/$(APPNAME)

OBJS = boundary_conditions.o flux_x.o initial_conditions.o \
       spectral_radius.o writeout.o ../sample_mains/main_1d_FD2.o

DIRS = $(SAMPLEDIR) $(APPDIR) $(APPDIR)/rho_files $(APPDIR)/u1_files \
       $(APPDIR)/u2_files $(APPDIR)/u3_files $(APPDIR)/b2_files \
       $(APPDIR)/b3_files $(APPDIR)/p_files $(APPDIR)/t_files $(APPDIR)/mesh_files

samples all: $(APPNAME)

clean:
	rm $(OBJS)

$(DIRS):
	mkdir $@

$(OBJS): %.o : %.cc
	g++ -c -fPIC -I$(CPROOT)/include $< -o $@

$(APPNAME): $(DIRS) $(OBJS)
	g++ $(OBJS) -I$(CPROOT)/include -L$(CPROOT)/lib -lcentpack_1d_FD2 -larray -o $@ -O3

	mv $(APPNAME) $(APPDIR)
	cp input $(APPDIR)
