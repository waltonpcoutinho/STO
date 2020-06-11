List of required packages

g++:   try to install v.7

gfortran: sudo apt-get install libgfortran3

AMPL:  install on /opt/AMPL/
       instal amplapi inside /opt/AMPL/
       add ampl.lic inside /opt/AMPL/
       all solvers' binaries go inside /opt/AMPL/
       all license and parameter files go inside /opt/AMPL/
       export path to AMPL binary: export PATH=$PATH:/opt/AMPL
       export path to AMPL-API library: 
       export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/AMPL/amplapi/lib/

CPLEX: install on default location and edit makefile accordingly
       preferred version: 12.4

ADOLC: install automake first 'sudo apt-get install automake'
       follow instructions from https://github.com/coin-or/ADOL-C
       export path to ADOL-C library
       export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/walton/adolc_base/lib64/

       
WORHP: get from https://worhp.de/
       export path to WORHP library
       export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/WORHP/worhp_1.13-2_linux/lib/
       add libworhp.so to /opt/AMPL/
       export path to WORHP parameters and license files
          WORHP_PARAM_FILE=/opt/AMPL/worhp.xml
          export WORHP_PARAM_FILE
          WORHP_LICENSE_FILE=/opt/AMPL/worhp.lic
          export WORHP_LICENSE_FILE

IPOPT: get from https://ampl.com/products/solvers/open-source/
       add binary to /opt/AMPL/
       
