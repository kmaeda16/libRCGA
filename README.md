
*****************
  libRCGA-1.2.2
*****************


Thank you for using libRCGA!

libRCGA is a C library for real-coded genetic algorithms (RCGAs). Currently, two RCGAs, UNDX/MGG and REXstar/JGG, are implemented. For constrained optimization problems, the stochastic ranking can be used. RCGAs are paralleled by MPI. For details, see the original paper: Kazuhiro Maeda, Fred C. Boogerd, and Hiroyuki Kurata, libRCGA: a C library for real-coded genetic algorithms for rapid parameter estimation of kinetic models, IPSJ Transactions on Bioinformatics, 11: 31-40, 2018 (https://www.jstage.jst.go.jp/article/ipsjtbio/11/0/11_31/_article/-char/en).


License:

libRCGA is distributed under GNU General Public License v3.0. For academic usage, libRCGA is free. For other usages, please contact the author(s).


Release Note:

Jan 11 2018: libRCGA-1.0 released.
Jan 14 2018: libRCGA-1.1 released. Bug fixed.
Aug  8 2018: libRCGA-1.2 released. Variable name changed. Exit codes changed.
Sep  5 2018: libRCGA-1.2.1 released. Citation information corrected.
Sep 14 2018: libRCGA-1.2.2 released. Article URL added.

Requirements:

- UNIX-like operating systems. The authors tested Cygwin on Windows 7 and 10, SUSE Linux Enterprise Server 11 SP4, Red Hat Enterprise Linux Server 6.9, and macOS Sierra.
- Modern C compilers. The authors tested GCC-4.2.1, GCC-4.8.2, and ICC-16.0.1.
- MPI. The authors tested OpenMPI-1.8.4, OpenMPI-2.1.1, OpenMPI-2.1.2, MPICH2-1.4, and SGI MPI. To use parallel versions of RCGAs, MPI needs to be installed. For serial versions only, MPI is not required.


Installation:

1. Extract libRCGA-x.x.x.tar.gz anywhere you like to install libRCGA by typing "tar xzf libRCGA-x.x.x.tar.gz" at your command prompt.
2. Change the current directory to "libRCGA-x.x.x" by typing "cd libRCGA-x.x.x".
3. Modify "Makefile" if you need, and type "make". If you use GCC with -O2, then no modification is needed.
4. If you find libundxmgg_serial.a, libundxmgg_parallel.a, librexstarjgg_serial.a, librexstarjgg_parallel.a, in the directory "lib", libRCGA is successfully installed. If you do not have MPI on your machine, *_parallel.a cannot be created. Still, you can use the serial version of libRCGA, i.e. *_serial.a.


Quick Start:

Four example source codes (example_*_*.c) in the directory "example" illustrate how to use libRCGA and can be used as a template for your program. To test whether the installed libRCGA works, please follow the directions provided in example_*_*.c, and compile and execute these examples.


Uninstallation:

Installation of libRCGA does not change any files or directories outside libRCGA-x.x.x. Thus, all you have to do for uninstallation is to delete libRCGA-x.x.x by typing "rm -r libRCGA-x.x.x".


Citation:

Kazuhiro Maeda, Fred C. Boogerd, and Hiroyuki Kurata, libRCGA: a C library for real-coded genetic algorithms for rapid parameter estimation of kinetic models, IPSJ Transactions on Bioinformatics, 11: 31-40, 2018


Good luck! Any suggestions and bug reports are welcome. Please contact KM.


-------------------------------
Kazuhiro Maeda
Kyushu Institute of Technology
