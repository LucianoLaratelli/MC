# NVT Monte Carlo

This is a program I wrote in late 2016, once I had started working at the
[space lab](http://drbrian.space/) as an undergraduate. Dr. Space had a list of exercises for new students to his lab that introduced us to the basics of computational science; this program is an extension of those initial exercises.

I have since cleaned it up a little bit (mostly formatting and removing
particularly poor code snippets.) An "NVT" simulation attempts to describe the
behavior of a group of atoms/molecules/etc. (we referred to them as
`atomecules`!) under conditions of constant volume, constant temperature, and
with a fixed particle count. These conditions together make up the NVT, or
[`canonical`](https://en.wikipedia.org/wiki/Canonical_ensemble) ensemble.



## Using the program

Most any old C compiler should work just fine. Compile the main NVT program with
`make`. The program also requires a starting configuration, which
`startgenerator.c` provides. Compile that as `gcc startgenerator.c -lm`.

The startgenerator executable will walk you through generating an input file for
the main NVT program. After that, run the simulation as `./NVT`. This will
prompt for an amount of iterations. As of this writing (2020-06-13), my CPU
(Ryzen 5 3600x) can run 50,000 iterations of 125 Argon atoms in a little under
fifteen minutes.

The main NVT program outputs three files: `positions.xyz`, `energies.dat`, and
`free_energies.dat`. At the time, I would commonly visualize `xyz` files with
[VMD](http://www.ks.uiuc.edu/Research/vmd/). This still appears to be a good
solution; after having installed it, run `vmd positions.xyz`. I recommend
adjusting the Drawing Method as follows: `Graphics->Representations...->Drawing
Method->VDW->Apply`. This will provide a nice atomistic view of the simulation:

![image](equil.png)
