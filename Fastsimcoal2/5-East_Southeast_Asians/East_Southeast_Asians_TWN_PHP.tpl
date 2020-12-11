//Parameters for the coalescence simulation program : fastsimcoal2.exe
11 samples to simulate
//Population effective sizes (number of genes)
5555
6774
1991
11765
58036
10205
7440
N_HAN
NPNG$
N_Fil
N_Paiwan
//Samples sizes
0
0
0
0
0
0
0
4
10
4
10
//Growth rates : negative growth implies population expansion RatioASW logRASW
0
0
0
0
0
0
0
0
0
0
0
//Number of migration matrices : 0 implies no migration between demes
5
// Matrix 0
0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 1.637439119815891e-04 0 0 0 0 0
0 0 0 0 1.637439119815891e-04 0 7.486981675412405e-04 0 0 0 0
0 0 0 0 0 7.486981675412405e-04 0 7.883476124205650e-06 0 0 0
0 0 0 0 0 0 7.883476124205650e-06 0 0 0 MIG107
0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 MIG109
0 0 0 0 0 0 0 MIG710 0 MIG910 0
// Matrix 1
0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 1.637439119815891e-04 0 0 0 0 0
0 0 0 0 1.637439119815891e-04 0 7.486981675412405e-04 0 0 0 0
0 0 0 0 0 7.486981675412405e-04 0 7.883476124205650e-06 0 0 0
0 0 0 0 0 0 7.883476124205650e-06 0 0 0 MIGa107
0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 MIGa107 0 0 0
// Matrix 2
0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 1.637439119815891e-04 0 0 0 0 0
0 0 0 0 1.637439119815891e-04 0 7.486981675412405e-04 0 0 0 0
0 0 0 0 0 7.486981675412405e-04 0 7.883476124205650e-06 0 0 0
0 0 0 0 0 0 7.883476124205650e-06 0 MIGa78 0 0
0 0 0 0 0 0 0 MIGa78 0 0 0
0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0
//Matrix 3
0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 1.637439119815891e-04 0 0 0 0 0
0 0 0 0 1.637439119815891e-04 0 7.486981675412405e-04 0 0 0 0
0 0 0 0 0 7.486981675412405e-04 0 0 5.207319988995006e-04 0 0
0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 5.207319988995006e-04 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0
//Matrix 4
0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0
//historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index
19  historical event
TD_TaiwFil  9 10 1 resizeSEA 0 1
TBOT$ 9 9 0 BOTFIL$ 0 1
TBOT$ 10 10 0 BOTTAIW$ 0 1
TD_SEAHAN   10 7 1 resizeHAN 0 2
1787 8 1 3.607697028777647e-02 1 0 keep
1662 7 6 1 6.065439378685459e+00 0 3
1840 6 3 3.587077613641241e-03 1 0 keep
1901 6 6 0 8.246982666635051e-02 0 4
2008 6 5 1 1 0 4
1933 8 8 0 RESIZENG 0 4
2041 8 5 1 1 0 4
2154 5 3 2.213573856721898e-02 1 0 4
2226 5 5 0 8.219752568672678e-02 0 4
2333 5 5 0 1.216581632653061e+01 0 4
4693 5 4 1 8.550070053830838e-01 0 4
14348 0 1 1 1 0 4
4302 2 3 1 1 0 keep
17380 1 3 1 1 0 4
21346 4 3 1 3.480902146234995e+00 0 4
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
FREQ 1 0 1.25e-8 OUTEXP
