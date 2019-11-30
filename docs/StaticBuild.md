Some ressources on how to have either fully static build or limited runtime/dynamic library dependencies.

Typically, this is useful for binary distribution/package maker, or in "big iron"-type HPC cluster
 where only a limited set of  dynamic libraries are installed/desirable. Can also be nice for small
 CLI frontend that may benefit from being disseminated as a standalone binary.
 
 
# Some useful commands

## Finding the runtime dependencies of an executable/library

The first move is examining file type with file, distinguished statically linked executable from other executables.

    $ file ./smoldock
    smoldock: ELF 64-bit LSB executable, x86-64, version 1 (GNU/Linux), statically linked, for GNU/Linux 3.2.0, with debug_info, not stripped
    $ file ./mpicalibrator
    mpicalibrator: ELF 64-bit LSB pie executable, x86-64, version 1 (GNU/Linux), dynamically linked, interpreter /lib64/ld-linux-x86-64.so.2, for GNU/Linux 3.2.0, with debug_info, not stripped

Stripping debug info may drastically reduce the disk size of the binary but : it has little impact on actual RAM footprint,
and the debug info allow useful strack traces and debugging in the wild. I think even optimized build can
benefit (relwithdebinfo in CMake parlance), and it facilitates greatly bug reporting.

If you have a static executable, you're already good to go. If not, you will want to examine which dynamic libraries are needed
for this binary :

    $ ldd mpicalibrator
        linux-vdso.so.1 (0x00007ffec3373000)
        libpthread.so.0 => /lib64/libpthread.so.0 (0x00007fc774ece000)
        libmpi_cxx.so.20 => /usr/lib64/libmpi_cxx.so.20 (0x00007fc774cb3000)
        libmpi.so.20 => /usr/lib64/libmpi.so.20 (0x00007fc7749be000)
        libgfortran.so.5 => /usr/lib/gcc/x86_64-pc-linux-gnu/8.2.0/libgfortran.so.5 (0x00007fc77454b000)
        libstdc++.so.6 => /usr/lib/gcc/x86_64-pc-linux-gnu/8.2.0/libstdc++.so.6 (0x00007fc774144000)
        libm.so.6 => /lib64/libm.so.6 (0x00007fc773dca000)
        libgcc_s.so.1 => /usr/lib/gcc/x86_64-pc-linux-gnu/8.2.0/libgcc_s.so.1 (0x00007fc773bb2000)
        libc.so.6 => /lib64/libc.so.6 (0x00007fc7737e7000)
        /lib64/ld-linux-x86-64.so.2 (0x00007fc77589b000)
        libopen-rte.so.20 => /usr/lib64/libopen-rte.so.20 (0x00007fc77355c000)
        libopen-pal.so.20 => /usr/lib64/libopen-pal.so.20 (0x00007fc7732c8000)
        libdl.so.2 => /lib64/libdl.so.2 (0x00007fc7730c4000)
        librt.so.1 => /lib64/librt.so.1 (0x00007fc772ebc000)
        libutil.so.1 => /lib64/libutil.so.1 (0x00007fc772cb9000)
        libhwloc.so.5 => /usr/lib64/libhwloc.so.5 (0x00007fc772a7d000)
        libevent-2.1.so.6 => /usr/lib64/libevent-2.1.so.6 (0x00007fc77282a000)
        libevent_pthreads-2.1.so.6 => /usr/lib64/libevent_pthreads-2.1.so.6 (0x00007fc772627000)
        libquadmath.so.0 => /usr/lib/gcc/x86_64-pc-linux-gnu/8.2.0/libquadmath.so.0 (0x00007fc7723e7000)
        libz.so.1 => /lib64/libz.so.1 (0x00007fc7721d0000)
        libnuma.so.1 => /usr/lib64/libnuma.so.1 (0x00007fc771fc5000)
        libudev.so.1 => /lib64/libudev.so.1 (0x00007fc771da0000)
        libpciaccess.so.0 => /usr/lib64/libpciaccess.so.0 (0x00007fc771b96000)
        libxml2.so.2 => /usr/lib64/libxml2.so.2 (0x00007fc771835000)
        libicuuc.so.63 => /usr/lib64/libicuuc.so.63 (0x00007fc771471000)
        libicudata.so.63 => /usr/lib64/libicudata.so.63 (0x00007fc76f883000)
    $ objdump -p mpicalibrator | grep NEEDED
      NEEDED               libpthread.so.0
      NEEDED               libmpi_cxx.so.20
      NEEDED               libmpi.so.20
      NEEDED               libgfortran.so.5
      NEEDED               libstdc++.so.6
      NEEDED               libm.so.6
      NEEDED               libgcc_s.so.1
      NEEDED               libc.so.6
      NEEDED               ld-linux-x86-64.so.2
      
You have two different output for those commands. ldd "fake loads" the given binary just like the loader would if executed,
 and resolves the indirect dependencies (ie dynamic library which requires other dynamic libraries). objdump -p grepping for NEEDED
 on the other hand only shows the direct dependencies, the library that the executable explicitly ask for.
 
Generally speaking, we only care about direct dependencies, because we assume that if these are correctly installed, their own
dependencies have been pulled in and installed too. If not, the system would be broken anyway. Furthermore, indirect dependencies
may change with version/implementations, and we should not be aware of that.
 
For example, libmpi.so.20 is directly required. The
OpenMPI implementation is pulling a lot of stuff like udev, hwloc, netlink libraries. But MPICH or other competing implementation may
use other hardware abstraction layer. We just depend on the libmpi.so and do not care about the pulled-in libraries.

# The fully-static to fully-dynamic continuum

We could dynamically link all the libraries we are using. It would be fine for package-based distribution,
as the package manager would install all the dependencies without problem. However, for other release channels
this is suboptimal : we would need to ship a lot of .so file with the binary, and possibly have wrapper script to
change LD_LIBRARY_PATH/LD_PRELOAD and other ugly hack.

On the other hand, fully static binary have major drawbacks. If you are statically linked to a certain version of the glic,
you will still need to have the same version installed on the deployment system to use rather common functions like
dlopen() or gethostbyname(). Statically linking only works when the binary is really "a kingdom unto itself", and only
crunch numbers, write to console and files, and that's it.

In my opinion, a nice compromise is to statically link the application-specific libraries, and leave
the system or general libraries as dynamic link. For example, we tuck in our table pretty printing library (brpinter),
but we still dynamically load the glibc, pthread or libstdc++.

For MPI application, is it also necessary to depend on a dynamic libmpi, because the version/implementation may change
in production. (MPICH has ABI compatiblity for example, requiring no recompilation and only switching of the .so file)

For general-yet-not-default libraries like Boost, the choice is more difficult. I prefer to statically link them but that's
a personal taste.

For python modules (PySmoldock.so), we need to dynically load the python libraries obviously. If we put Boost as static dependency,
we can mostly just ship the .so and users can import it without fuss. So that's what Smoldock is doing.


# Statically-linked C++ frontends

You can very easily statically link the C++ CLI frontend, or your dedicated front-end.

Smoldock provides libsmoldock.a, which contain a partially-linked object file. This is the result of a special
invokation of the linker, which will proceed as if it was linking an executable for the given libraries. This means that
given the libsmoldock.a, you can link your own program without worrying about the internal dependencies of Smoldock itself :
for example you don't need RDKit, just the smoldock headers and libsmoldock.a.

This is used to build the "smoldock" front-end :

    add_executable(smoldock Frontends/main.cpp)
    set_target_properties(smoldock PROPERTIES LINK_FLAGS " -Wl,--no-undefined -Wl,--as-needed -static  -Wl,-allow-multiple-definition   ")
    TARGET_LINK_LIBRARIES(smoldock libsmoldock)
    TARGET_LINK_LIBRARIES(smoldock pthread)
    target_link_libraries(smoldock boost_log bprinter tbb_static)
        
We link with libsmoldock, and then we link with the libraries needed by the front-end itself : if the front end was not
calling into Boost log, nor tbb or the table pretty-printer library, we could just link with libsmoldock.

There is also the possiblity, evoked above, to statically link libsmoldock but keep some/all the dependencies of the front end
as dynamic libraries. This is done for the smoldock_dyn frontend which is exactly the same except for the CMake target :


    add_executable(smoldock_dyn Frontends/main.cpp)
    set_target_properties(smoldock_dyn PROPERTIES LINK_FLAGS " -Wl,--no-undefined -Wl,--as-needed    ")
    TARGET_LINK_LIBRARIES(smoldock_dyn libsmoldock)
    TARGET_LINK_LIBRARIES(smoldock_dyn pthread)
    target_link_libraries(smoldock_dyn -Wl,-Bstatic boost_log bprinter tbb_static  -Wl,-Bdynamic)

Which has the following direct dependencies :

    $ objdump -p ./smoldock_dyn | grep NEEDED
      NEEDED               libpthread.so.0
      NEEDED               libstdc++.so.6
      NEEDED               libm.so.6
      NEEDED               libgcc_s.so.1
      NEEDED               libc.so.6
      NEEDED               ld-linux-x86-64.so.2


You can follow the same modus operandi for your own C++ frontends.

# Semi-static python modules

As discussed above, it would be nice to make the PySmoldock.so python module less dependant on dependencies of
Smoldock itself.

However, there is stark difference between making an executable and a shared library in term of linking. The
executable does not neet to be position independant : granted you dont use -fPIE, the .text and .data sections containing the 
machine code and variables will have a fixed adress and all call/memory access will be trough absolute adressing.

In terms of linking, this means the code is position-dependant, the contrary of position independant. But a shared library can and will
be mapped anywhere in the adress space, and thus need position independant adressing.

We achieve position independant adressing by adding -fPIC to the compile options. This is independant from static/dynamic linkage :
it is wholly possible to make a static libsmoldock.a with position independant code (compiled with -fPIC). However, it is quite
rare for library to be provided as static .a archive **and** compiled with -fPIC. Thus we will have a problem when trying to tuck in
our dependencies : either we use dynamic libraries (with the corresponding problem of distributing .so) or we try to statically link (-static and variant)
but fail to get an executable (the linker will complain that static-linking-type absolute adressing/relocation are not emitable in a dynamically
loaded library, as it should).

So we mostly cannot use our distribution package for our dependencies : we have to build our own static libraries from source,
taking care to add the '-fPIC' to the compilation flags. It is unclear to me whether we would need to rebuild dependencies of
dependencies in this way, but most of Smoldock dependencies are somewhat free standing.

At this point I have not completly solved this problem (ie I am not able to get a python module .so with minimal dynamic dependencies), but
I am working on it and some useful work on that front follows :

## Boost 
This seems to work :

    ./bootstrap.sh --prefix=/home/builder/StaticPic/Boost/

    ./b2 --prefix=/home/builder/StaticPic/Boost/ --build-type=complete --build-dir=/home/builder/boost_1_69_0/build/ --layout=tagged cxxflags=" -fPIC " cflags=" -fPIC " --ignore-site-config

I still have questions regarding how boost python is built in this scenario.

## RDKit

