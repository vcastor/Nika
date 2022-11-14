# Welcome to HFMP2DFT :rocket:

Nika Is Not an Known Acronym

This is the README file for NIKA project, this program can compute
the HF, MP2 and DFT calculations.

## Compiling

Just type: `make` on your terminal at the directory where is README file.

If that didn't work correctly I'll give you answer of the problems that
you might have.

First of all you need a fortran compiler, I estrictly recomend gfortran, even
if its written on C/C++, is more robust and extended use than other compilers. Intel Fortran
Compiler, NVIDIA nvfortran and AMD flang are developed just
specifically for their owns architectures.

Also gfortran is free and can be used on Linux, macOS, and its which
we used :grimacing:. However, you can change it for your favorite compiler in the
`Makefile`, if you don't have any compiler, once again, I recomend GNU Compiler
Collection (its not only fortran compiler, can also compile: `C`, `C++`,
`Objective-C`, `Ada`, `Go`, and `D`).

<hr style="border:2px solid gray"> </hr>

The program was developed using
[LAPACK](http://www.netlib.org/lapack/#_lapack_version_3_10_0) library, if your
computer doesn't have the library you will have an error like:

>`"_dsygvd_", referenced from:`
> `     _MAIN__ in ccZlxZJJ.o`

The problem can be fixed just intalling LAPACK library. In case that
you have it but doesn't compailing neither, maybe you have it in not
the usual directory that I'm supposing, just modify the `Makefile` at line
where we call the library in the compaling.

<hr style="border:2px solid gray"> </hr>

Fundamental Theorem of Computing: turn off, turn on and try again. That would
fix also the error of LAPACK if you just installed on the same terminal session
and maybe is necessary to restart the `$PATH` turning off/on the terminal
session.

## How run a calculation :wink:

On the same directory is a laucher call it `launcher.sh` (I have a lot of
creativity I know :blush:). The laucnher would have the excecutable permissions, if it
doesn't have give it with the next command:

> `chmod +x launcher.sh`

To send the calculation just give the input as an argument of the launcher.
Example:

> `$ ./launcher.sh H2.inp &`

The `&` just to not freez the terminal, its doesn't matter the name, but the
extension should be `inp`, that becuase the launcher will create the output
file with the same name but with the extention `out` and a scratch file as
`moreout`.

## To see the source code :beers:

If you want to see the source code, take acount that will be easer to read if
you have a code editor as Sublime Text, Atom, Visual Studio Code, ... or just
vi/vim with **syn on**, that to have highlighted syntaxes on reserved words.

