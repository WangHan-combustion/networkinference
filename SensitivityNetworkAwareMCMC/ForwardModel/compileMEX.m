% Mac
%mex main.cpp SOFCModel.cpp SOFCAnode.cpp SOFCAnodeIntegrator.cpp tools.cpp -I/Applications/Cantera/include -L/Applications/Cantera/lib -L/usr/local/lib -L/usr/local/lib -L/usr/local/gfortran/lib -luser -loneD -lzeroD -lequil -lkinetics -ltransport -lthermo -lctnumerics -lctmath -ltpx -lctspectra -lconverters -lctbase -lsundials_cvodes -lsundials_nvecserial -lpthread -lctf2c -lctcxx -lgfortran -lgsl -lgslcblas -llapack -lm -lblas

% Desktop
mex main.cpp SOFCModel.cpp SOFCAnode.cpp SOFCAnodeIntegrator.cpp tools.cpp -I/usr/local/cantera/include -I/usr/include/mpich2 -L/usr/local/cantera/lib -L/usr/local/lib -L/usr/local/gfortran/lib -luser -loneD -lzeroD -lequil -lkinetics -ltransport -lthermo -lctnumerics -lctmath -ltpx -lctspectra -lconverters -lctbase -lsundials_cvodes -lsundials_nvecserial -lpthread -lctf2c -lctcxx -llapack -lblas -lm

% Matlab Cluster
%mex main.cpp SOFCModel.cpp SOFCAnode.cpp SOFCAnodeIntegrator.cpp tools.cpp -I/master/home/nikhilg1/localDryReforming/include -L/master/home/nikhilg1/localDryReforming/lib -L/usr/lib/gcc/x86_64-linux-gnu -L/master-usr/lib -L/usr/share/doc -luser -loneD -lzeroD -lequil -lkinetics -ltransport -lthermo -lctnumerics -lctmath -ltpx -lctspectra -lconverters -lctbase -lsundials_cvodes -lsundials_nvecserial -lpthread -lctf2c -lctcxx -llapack -lm -lblas