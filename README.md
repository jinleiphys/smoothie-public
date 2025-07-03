# smoothie
 **smoothie**(**S**cattering **M**odel of **O**ptical **O**perator **Th**eory for **I**chimura-Austern-Vincent **E**quations) is a computer code developed by Jin Lei and Antonio M. Moro to perform the Non-elastic breakup calculations in inclusive breakup reactions using the formalism of Ichimura, Austern and Vincent. It uses Fortran 95 and runs on Linux and Windows machines.   

Smoothie assumes an inclusive process is of the form :

*a(=b+x) + A -> b + B*,  where *B= (x+A)**


## Input description for DWBA-NEB calculations
### &GLOBAL namelist: hcm, lmax, elab, thmin, thmax, thinc, icf, cutl, lmin, nx, rmax, nr, prior, printf, jtmax, jtmin, dwba

- hcm: radial step for scattering (bound state) wave functions when solving the Schordinger equations.
- lmax: maximum orbital angular momentum for a+A and b+B channels
- elab: lab energy of the projectile "a"
- thmin, thmax: minimum, maximum angle of outgoing b particle
- thinc: angular step
- icf: logical variable, if icf=T, perform  <a href="https://www.codecogs.com/eqnedit.php?latex=\langle&space;\psi_x&space;|&space;W_{fus}|&space;\psi_x\rangle" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\langle&space;\psi_x&space;|&space;W_{fus}|&space;\psi_x\rangle" title="\langle \psi_x | W_{fus}| \psi_x\rangle" /></a>
- cutl: maximum orbital angular momentum for x+A channel. If not defined in the input file, cutl=lmax
- nx: grid number for the angular integration by using the Gaussian quadrature method
- rmax: maximum number of the radial distance for all channels
- nr: grid number for the radial integration by using the Gaussian quadrature method
- prior: logical variable, if prior=T, use the IAV-prior form, otherwise use IAV-post form
- printf: logical variable, if printf=T, print the wave function for x-A channel and the source term. Default value printf=F
- jtmax: maximum value of total J used in the calculation. If not defined in the input, jtmax=2*lmax
- jtmin: minimum value of total J used in the calculation. If not defined in the input, jtmin=0
- dwba: control variable for different numerical method for the DWBA calculation.

      1. dwba=1, DWBA calculation without intrinsic spins and choosing rbx as integration variable
      2. dwba=2, DWBA calculation without intrinsic spins and choosing rb as integration variable
      3. dwba=3, DWBA calculation with spins and choosing rb as integration variable
      4. dwba=4, same as dwba=1, but use Lagrange-mesh method

### &SYSTEM namelist: namep, massp, zp, jp, namet, masst, zt, jt, nameb, massb, zb, jb, namex, massx, zx, jx, qval, neb, lbx, sbx,nodes
- namep, namet, nameb, namex: characters, defines the name of each particle
- massp, masst, nassb, massx: mass of each particle
- zp, zt, zb, zx: charge of each particle
- jp, jt, jb, jx: spin of each particle
- sbx: coupled jb+jx spin
- qval: bound state energy of the b-x system
- neb: number of energy points for outgoing b particle. If not specified, the code use neb=nint((ecmbmax-ecmbmin)/ecmbh) + 1
- lbx: orbital angular momentum of b-x inside the projectile bound state
- nodes: nodes of the bound state wave function of the b-x system
- gswf: if different from blank, filename to read projectile ground-state wavefunction from external file. 
  Format must be: 1st line: custom comment; 2nd line: number of radial points; next lines: radius and wavefunction (including "r" factor)

### &OUTPUT namelist: ecmbmin,ecmbmax,ecmbh,bin,wbin
- ecmbmin: minimum value of outgoing b energy in b - B relative coordinate
- ecmbmax: maximum value of outgoing b energy in b - B relative coordinate
- ecmbh: energy step
- bin: If bin/=0, integer number of bins for bin method for b - B channel.
- wbin: bin width 

### &POTENTIAL namelist: kp1, kp2, ptype, a1, a2, rc, uv, av, rv, uw, aw, rw, vsov, rsov, asov, vsow, rsow, asow, vd, avd, rvd, wd, awd, rwd, renv, renw, renvd, renwd
- kp1: characters, only equals to the following values

      1. "a" : define the potential for *a+A* channel
      2. "b" : define the potential for *b+B* channel
      3. "t" : define the potential for *b+A* channel  
      4. "x" : define the potential for *x+A* channel
      5. "f" : define the fusion potential for *x+A* channel when icf=T

- kp2: index, normally set to kp2=1

- a1, a2, rc, uv, av, rv, uw, aw, rw, vsov, rsov, asov, vsow, rsow, asow, vd, avd, rvd, wd, awd, rwd : potential parameters
- renv, renw, renvd, renwd, revls, renwls: renormalization factors for WS potentials

- ptype: defines different types of potential

      1. Woods-Saxon potential.  Potential parameters are defined through the variables: 
         * Volume WS: uv, av, rv, uw, aw, rw
         * Spin-orbit: vsov, rsov, asov, vsow, rsow, asow
         * Surface (derivative) WS: vd, avd, rvd, wd, awd, rwd
      2. Gaussian potential, defined through the parameters: uv,rv,av,uw,rw,aw
      3. Nucleon-nucleus Koning-Delaroche (KD02)
      4. Nucleon-nucleus by Varner et al (CH89)
      5. Nucleon-nucleus Beccetti-Greenless (bgPN)
      6. Deuteron OMP for 12< A< 209  "yyq06"  from Y. Han, Phys Rev C 74 044615,  removing the surface term of imaginary potential
      7. NOT USED
      8. Nucleon-actinide OMP "yyh10" from Y.Han, Phys.Rev. C 81, 024616 (2010)
      9. 3H/3He from Beccetti-Greenless (BG69)
      10. Deuteron OMP for 12< A< 209  "yyq06"  from Y. Han, Phys Rev C 74 044615
      11. Akyuz-Winther potential
      12. n-9Be energy-dependent potential by Angela Bonaccorso 
      13. alpha global potential International taken from Journal of Modern Physics EVol. 24, No. 12 (2015) 1550092
      14. BG69 with uw*0.8
      15. Read table with OMP parameters as a function of energy (originally for Morillon's DOM potential for n-64Zn)
      16. Guo triton potential (old) for 3H+120Sn
      17. Guo triton potential (new) for 3H+120Sn
      19. 3He/3H GDP08 global potential from Pang et al, Phys. Rev. C 79, 024615
      20. Avrigeanu alpha potential (2010)
      41. read in potential from fort.41
      42. read in potential from fort.42
      43. read in potential from fort.43
      44. read in potential from fort.44


## Compilation
To Install the **smoothie** program
- Edit and customize the make.inc file in the main directory
- Enter the smoothie subdirectory
- make
To clean old files and libraries type:

`prompt> make clean`

## Test examples 

### 93Nb(d,pX) 

NAMELIST

&GLOBAL      hcm=0.05  lmax=23  elab=25.5 thmin=0. thmax=180. prior=t printf=t dwba=4 cutl=12
             thinc=1   nx=34 rmax=50   nr=120 /


&SYSTEM     namep='d'     massp=2.       zp=1.0    jp=0.0
            namet='93Nb'  masst=93.0     zt=41.0   jt=0.0  
            nameb='p'     massb=1        zb=1      jb=0.0  neb=1
            namex='n'     massx=1.0086   zx=0.0    jx=0.0  lbx=0    nodes=1 qval=-2.224 /

&OUTGOING   ecmbmin=10 ecmbmax=10 ecmbh=1  bin=0  wbin=0.0 /

&OUTGOING /


&POTENTIAL  kp1='a'  kp2=1 ptype=10 a1=0 a2=93
            /


&POTENTIAL  kp1='b'  kp2=1 ptype=3 a1=0 a2=94
           /


&POTENTIAL  kp1='x'  kp2=1 ptype=3 a1=0 a2=93
            /

&POTENTIAL  kp1='p'  kp2=1 ptype=2 a1=1 a2=1 rc=1.5
            uv=72.15 av=1.484
            /

&POTENTIAL  kp1='t'  kp2=1 ptype=3 a1=0 a2=93
           /

&POTENTIAL /
