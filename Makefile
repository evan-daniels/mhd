all:
	make libarray
	make libcp_1d_FD2
	make libcp_1d_SD2
	make libcp_2d_FD2
	make libcp_2d_SD2
	make burgers_FD2
	make burgers_SD2
	make euler1d_FD2
	make euler1d_SD2
	make MHD1d_FD2
	make MHD1d_SD2
	make scalar2d_FD2
	make scalar2d_SD2
	make euler2d_FD2
	make euler2d_SD2
	make MHD2d_FD2
	make MHD2d_SD2

all_libs:
	make libarray
	make libcp_1d_FD2
	make libcp_1d_SD2
	make libcp_2d_FD2
	make libcp_2d_SD2

all_samples:
	make burgers_FD2
	make burgers_SD2
	make euler1d_FD2
	make euler1d_SD2
	make MHD1d_FD2
	make MHD1d_SD2
	make scalar2d_FD2
	make scalar2d_SD2
	make euler2d_FD2
	make euler2d_SD2
	make MHD2d_FD2
	make MHD2d_SD2

libarray:
	g++ -c -fPIC src/common/doublearray*.cc -Iinclude
	ar cr libarray.a *.o
	rm *.o
	mv libarray.a lib/

libcp_1d_FD2:
	mpicxx -c -fPIC src/common/disclaimer.cc src/common/minmod3.cc src/1d/common/*.cc src/1d/FD2/*.cc -Ilib/ -Iinclude/
	ar cr libcp_1d_FD2.a *.o
	rm *.o
	mv libcp_1d_FD2.a lib/

libcp_1d_SD2:
	mpicxx -c -fPIC src/common/disclaimer.cc src/common/minmod3.cc src/1d/common/*.cc src/1d/SD2/*.cc -Ilib/ -Iinclude/
	ar cr libcp_1d_SD2.a *.o
	rm *.o
	mv libcp_1d_SD2.a lib/

libcp_2d_FD2:
	mpicxx -c -fPIC src/common/disclaimer.cc src/common/minmod3.cc src/2d/common/*.cc src/2d/FD2/*.cc -Ilib/ -Iinclude/
	ar cr libcp_2d_FD2.a *.o
	rm *.o
	mv libcp_2d_FD2.a lib/

libcp_2d_SD2:
	mpicxx -c -fPIC src/common/disclaimer.cc src/common/minmod3.cc src/2d/common/*.cc src/2d/SD2/*.cc -Ilib/ -Iinclude/
	ar cr libcp_2d_SD2.a *.o
	rm *.o
	mv libcp_2d_SD2.a lib/

burgers_FD2:
	mpicxx src/1d/samples/burgers/*.cc src/1d/sample_mains/main_1d_FD2.cc -Iinclude/ -Ilib/ -Linclude/ -Llib/ -lcp_1d_FD2 -larray -o burgers -O3

	mkdir simulation_burgers_FD2
	mkdir simulation_burgers_FD2/mesh_files
	mkdir simulation_burgers_FD2/u_frames
	mkdir simulation_burgers_FD2/u_files
	mkdir simulation_burgers_FD2/t_files
	mkdir simulation_burgers_FD2/odd_files

	mv burgers simulation_burgers_FD2
	cp src/1d/samples/burgers/input_fd simulation_burgers_FD2/input

burgers_SD2:
	mpicxx src/1d/samples/burgers/*.cc src/1d/sample_mains/main_1d_SD2.cc -Iinclude/ -Ilib/ -Linclude/ -Llib/ -lcp_1d_SD2 -larray -o burgers -O3

	mkdir simulation_burgers_SD2
	mkdir simulation_burgers_SD2/mesh_files
	mkdir simulation_burgers_SD2/u_frames
	mkdir simulation_burgers_SD2/u_files
	mkdir simulation_burgers_SD2/t_files

	mv burgers simulation_burgers_SD2
	cp src/1d/samples/burgers/input_sd simulation_burgers_SD2/input

euler1d_FD2:
	mpicxx src/1d/samples/euler1d/*.cc src/1d/sample_mains/main_1d_FD2.cc -Iinclude/ -Ilib/ -Linclude/ -Llib/ -lcp_1d_FD2 -larray -o euler1d -O3

	mkdir simulation_euler1d_FD2
	mkdir simulation_euler1d_FD2/mesh_files
	mkdir simulation_euler1d_FD2/rho_frames
	mkdir simulation_euler1d_FD2/rho_files
	mkdir simulation_euler1d_FD2/u1_frames
	mkdir simulation_euler1d_FD2/u1_files
	mkdir simulation_euler1d_FD2/p_frames
	mkdir simulation_euler1d_FD2/p_files
	mkdir simulation_euler1d_FD2/t_files
	mkdir simulation_euler1d_FD2/odd_files

	mv euler1d simulation_euler1d_FD2
	cp src/1d/samples/euler1d/input_fd simulation_euler1d_FD2/input

euler1d_SD2:
	mpicxx src/1d/samples/euler1d/*.cc src/1d/sample_mains/main_1d_SD2.cc -Iinclude/ -Ilib/ -Linclude/ -Llib/ -lcp_1d_SD2 -larray -o euler1d -O3

	mkdir simulation_euler1d_SD2
	mkdir simulation_euler1d_SD2/mesh_files
	mkdir simulation_euler1d_SD2/rho_frames
	mkdir simulation_euler1d_SD2/rho_files
	mkdir simulation_euler1d_SD2/u1_frames
	mkdir simulation_euler1d_SD2/u1_files
	mkdir simulation_euler1d_SD2/p_frames
	mkdir simulation_euler1d_SD2/p_files
	mkdir simulation_euler1d_SD2/t_files

	mv euler1d simulation_euler1d_SD2
	cp src/1d/samples/euler1d/input_sd simulation_euler1d_SD2/input

MHD1d_FD2:
	mpicxx src/1d/samples/MHD1d/*.cc src/1d/sample_mains/main_1d_FD2.cc -Iinclude/ -Ilib/ -Linclude/ -Llib/ -lcp_1d_FD2 -larray -o MHD1d -O3

	mkdir simulation_MHD1d_FD2
	mkdir simulation_MHD1d_FD2/mesh_files
	mkdir simulation_MHD1d_FD2/u_files
	mkdir simulation_MHD1d_FD2/t_files
	mkdir simulation_MHD1d_FD2/rho_files
	mkdir simulation_MHD1d_FD2/u1_files
	mkdir simulation_MHD1d_FD2/u2_files
	mkdir simulation_MHD1d_FD2/u3_files
	mkdir simulation_MHD1d_FD2/b2_files
	mkdir simulation_MHD1d_FD2/b3_files
	mkdir simulation_MHD1d_FD2/p_files
	mkdir simulation_MHD1d_FD2/odd_files
	mkdir simulation_MHD1d_FD2/rho_frames
	mkdir simulation_MHD1d_FD2/u1_frames
	mkdir simulation_MHD1d_FD2/u2_frames
	mkdir simulation_MHD1d_FD2/u3_frames
	mkdir simulation_MHD1d_FD2/b2_frames
	mkdir simulation_MHD1d_FD2/b3_frames
	mkdir simulation_MHD1d_FD2/p_frames

	mv MHD1d simulation_MHD1d_FD2
	cp src/1d/samples/MHD1d/input_fd simulation_MHD1d_FD2/input

MHD1d_SD2:
	mpicxx src/1d/samples/MHD1d/*.cc src/1d/sample_mains/main_1d_SD2.cc -Iinclude/ -Ilib/ -Linclude/ -Llib/ -lcp_1d_SD2 -larray -o MHD1d -O3

	mkdir simulation_MHD1d_SD2
	mkdir simulation_MHD1d_SD2/mesh_files
	mkdir simulation_MHD1d_SD2/u_files
	mkdir simulation_MHD1d_SD2/t_files
	mkdir simulation_MHD1d_SD2/rho_files
	mkdir simulation_MHD1d_SD2/u1_files
	mkdir simulation_MHD1d_SD2/u2_files
	mkdir simulation_MHD1d_SD2/u3_files
	mkdir simulation_MHD1d_SD2/b2_files
	mkdir simulation_MHD1d_SD2/b3_files
	mkdir simulation_MHD1d_SD2/p_files
	mkdir simulation_MHD1d_SD2/rho_frames
	mkdir simulation_MHD1d_SD2/u1_frames
	mkdir simulation_MHD1d_SD2/u2_frames
	mkdir simulation_MHD1d_SD2/u3_frames
	mkdir simulation_MHD1d_SD2/b2_frames
	mkdir simulation_MHD1d_SD2/b3_frames
	mkdir simulation_MHD1d_SD2/p_frames

	mv MHD1d simulation_MHD1d_SD2
	cp src/1d/samples/MHD1d/input_sd simulation_MHD1d_SD2/input

scalar2d_FD2:
	mpicxx src/2d/samples/scalar2d/*.cc src/2d/sample_mains/main_2d_FD2.cc -Iinclude/ -Ilib/ -Linclude/ -Llib/ -lcp_2d_FD2 -larray -o scalar2d -O3

	mkdir simulation_scalar2d_FD2
	mkdir simulation_scalar2d_FD2/mesh_files
	mkdir simulation_scalar2d_FD2/u_frames
	mkdir simulation_scalar2d_FD2/u_files
	mkdir simulation_scalar2d_FD2/t_files
	mkdir simulation_scalar2d_FD2/odd_files

	mv scalar2d simulation_scalar2d_FD2
	cp src/2d/samples/scalar2d/input_fd simulation_scalar2d_FD2/input

scalar2d_SD2:
	mpicxx src/2d/samples/scalar2d/*.cc src/2d/sample_mains/main_2d_SD2.cc -Iinclude/ -Ilib/ -Linclude/ -Llib/ -lcp_2d_SD2 -larray -o scalar2d -O3

	mkdir simulation_scalar2d_SD2
	mkdir simulation_scalar2d_SD2/mesh_files
	mkdir simulation_scalar2d_SD2/u_frames
	mkdir simulation_scalar2d_SD2/u_files
	mkdir simulation_scalar2d_SD2/t_files

	mv scalar2d simulation_scalar2d_SD2
	cp src/2d/samples/scalar2d/input_sd simulation_scalar2d_SD2/input

euler2d_FD2:
	mpicxx src/2d/samples/euler2d/*.cc src/2d/sample_mains/main_2d_FD2.cc -Iinclude/ -Ilib/ -Linclude/ -Llib/ -lcp_2d_FD2 -larray -o euler2d -O3

	mkdir simulation_euler2d_FD2
	mkdir simulation_euler2d_FD2/mesh_files
	mkdir simulation_euler2d_FD2/rho_frames
	mkdir simulation_euler2d_FD2/rho_files
	mkdir simulation_euler2d_FD2/u_frames
	mkdir simulation_euler2d_FD2/u1_files
	mkdir simulation_euler2d_FD2/u2_files
	mkdir simulation_euler2d_FD2/p_frames
	mkdir simulation_euler2d_FD2/p_files
	mkdir simulation_euler2d_FD2/t_files
	mkdir simulation_euler2d_FD2/odd_files

	mv euler2d simulation_euler2d_FD2
	cp src/2d/samples/euler2d/input_fd simulation_euler2d_FD2/input

euler2d_SD2:
	mpicxx src/2d/samples/euler2d/*.cc src/2d/sample_mains/main_2d_SD2.cc -Iinclude/ -Ilib/ -Linclude/ -Llib/ -lcp_2d_SD2 -larray -o euler2d -O3

	mkdir simulation_euler2d_SD2
	mkdir simulation_euler2d_SD2/mesh_files
	mkdir simulation_euler2d_SD2/rho_frames
	mkdir simulation_euler2d_SD2/rho_files
	mkdir simulation_euler2d_SD2/u_frames
	mkdir simulation_euler2d_SD2/u1_files
	mkdir simulation_euler2d_SD2/u2_files
	mkdir simulation_euler2d_SD2/p_frames
	mkdir simulation_euler2d_SD2/p_files
	mkdir simulation_euler2d_SD2/t_files

	mv euler2d simulation_euler2d_SD2
	cp src/2d/samples/euler2d/input_sd simulation_euler2d_SD2/input

MHD2d_FD2:
	mpicxx src/2d/samples/MHD2d/*.cc src/2d/sample_mains/main_2d_FD2.cc -Iinclude/ -Ilib/ -Linclude/ -Llib/ -lcp_2d_FD2 -larray -o MHD2d -O3

	mkdir simulation_MHD2d_FD2
	mkdir simulation_MHD2d_FD2/mesh_files
	mkdir simulation_MHD2d_FD2/rho_frames
	mkdir simulation_MHD2d_FD2/rho_files
	mkdir simulation_MHD2d_FD2/u_frames
	mkdir simulation_MHD2d_FD2/u1_files
	mkdir simulation_MHD2d_FD2/u2_files
	mkdir simulation_MHD2d_FD2/u3_files
	mkdir simulation_MHD2d_FD2/b_frames
	mkdir simulation_MHD2d_FD2/b1_files
	mkdir simulation_MHD2d_FD2/b2_files
	mkdir simulation_MHD2d_FD2/b3_files
	mkdir simulation_MHD2d_FD2/p_frames
	mkdir simulation_MHD2d_FD2/p_files
	mkdir simulation_MHD2d_FD2/t_files
	mkdir simulation_MHD2d_FD2/odd_files

	mv MHD2d simulation_MHD2d_FD2
	cp src/2d/samples/MHD2d/input_fd simulation_MHD2d_FD2/input

MHD2d_SD2:
	mpicxx src/2d/samples/MHD2d/*.cc src/2d/sample_mains/main_2d_SD2.cc -Iinclude/ -Ilib/ -Linclude/ -Llib/ -lcp_2d_SD2 -larray -o MHD2d -O3

	mkdir simulation_MHD2d_SD2
	mkdir simulation_MHD2d_SD2/mesh_files
	mkdir simulation_MHD2d_SD2/rho_frames
	mkdir simulation_MHD2d_SD2/rho_files
	mkdir simulation_MHD2d_SD2/u_frames
	mkdir simulation_MHD2d_SD2/u1_files
	mkdir simulation_MHD2d_SD2/u2_files
	mkdir simulation_MHD2d_SD2/u3_files
	mkdir simulation_MHD2d_SD2/b_frames
	mkdir simulation_MHD2d_SD2/b1_files
	mkdir simulation_MHD2d_SD2/b2_files
	mkdir simulation_MHD2d_SD2/b3_files
	mkdir simulation_MHD2d_SD2/p_frames
	mkdir simulation_MHD2d_SD2/p_files
	mkdir simulation_MHD2d_SD2/t_files

	mv MHD2d simulation_MHD2d_SD2
	cp src/2d/samples/MHD2d/input_sd simulation_MHD2d_SD2/input

MHD2d_reconnection:
	mpicxx src/2d/samples/MHD2d/*.cc src/2d/sample_mains/main_2d_SD2.cc -Iinclude/ -Ilib/ -Linclude/ -Llib/ -lcp_2d_SD2 -larray -o MHD2d -O3 -std=c++11

	mkdir -p simulation_MHD2d_reconnection_SD2
	mkdir -p simulation_MHD2d_reconnection_SD2/mesh_files
	mkdir -p simulation_MHD2d_reconnection_SD2/rho_frames
	mkdir -p simulation_MHD2d_reconnection_SD2/rho_files
	mkdir -p simulation_MHD2d_reconnection_SD2/u_frames
	mkdir -p simulation_MHD2d_reconnection_SD2/u1_files
	mkdir -p simulation_MHD2d_reconnection_SD2/u2_files
	mkdir -p simulation_MHD2d_reconnection_SD2/u3_files
	mkdir -p simulation_MHD2d_reconnection_SD2/b_frames
	mkdir -p simulation_MHD2d_reconnection_SD2/b1_files
	mkdir -p simulation_MHD2d_reconnection_SD2/b2_files
	mkdir -p simulation_MHD2d_reconnection_SD2/b3_files
	mkdir -p simulation_MHD2d_reconnection_SD2/p_frames
	mkdir -p simulation_MHD2d_reconnection_SD2/p_files
	mkdir -p simulation_MHD2d_reconnection_SD2/t_files
	mkdir -p simulation_MHD2d_reconnection_SD2/jz_frames
	mkdir -p simulation_MHD2d_reconnection_SD2/b_on_rho_frames
	mkdir -p simulation_MHD2d_reconnection_SD2/bz_frames

	mv MHD2d simulation_MHD2d_reconnection_SD2/
	cp src/2d/samples/MHD2d/input_reconnection simulation_MHD2d_reconnection_SD2/input