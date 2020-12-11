//Parameters for the coalescence simulation program : fastsimcoal2.exe
12 samples to simulate
//Population effective sizes (number of genes)
N_Daltai
N_ghost_D
N_Nvindija
N_ghost_V
54244
9538
N_Eu
N_As
N_PNG
N_Bun
N_Kun
N_Paiw
//Samples sizes
2 2800
0
2 2000
0
0
0
4
4
0
10
10
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
0
//Number of migration matrices : 0 implies no migration between demes
5
// Matrix 0
0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 1.59649e-04 0 0 0 0 0 0
0 0 0 0 1.59649e-04 0 MIG56 0 0 0 0 0
0 0 0 0 0 MIG56 0 MIG67 0 0 0 0
0 0 0 0 0 0 MIG67 0 0 0 0 MIG711
0 0 0 0 0 0 0 0 0 0 0 MIG811
0 0 0 0 0 0 0 0 MIGBUN$ 0 0 0
0 0 0 0 0 0 0 0 MIGKUN$ 0 0 0
0 0 0 0 0 0 0 MIG711 MIG811 0 0 0
//Matrix 1
0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 1.59649e-04 0 0 0 0 0 0
0 0 0 0 1.59649e-04 0 MIG56 0 0 0 0 0
0 0 0 0 0 MIG56 0 MIG67 0 0 0 0
0 0 0 0 0 0 MIG67 0 0 0 0 MIG711
0 0 0 0 0 0 0 0 0 0 0 MIG811
0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 MIG711 MIG811 0 0 0
// Matrix 1
0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 1.59649e-04 0 0 0 0 0 0
0 0 0 0 1.59649e-04 0 MIG56 0 0 0 0 0
0 0 0 0 0 MIG56 0 MIG67 0 0 0 0
0 0 0 0 0 0 MIG67 0 MIG78 0 0 0
0 0 0 0 0 0 0 MIG78 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0
// Matrix 2
0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 1.59649e-04 0 0 0 0 0 0
0 0 0 0 1.59649e-04 0 MIG56 0 0 0 0 0
0 0 0 0 0 MIG56 0 MIGa67 0 0 0 0
0 0 0 0 0 0 MIGa67 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0
// Matrix 3
0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0
//historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index
18  historical event
100 9 8 1 1 0 1
100 10 8 1 1 0 1
TD_AsTaiw 11 7 1 1 0 2
TAdm_Den 8 1 adDen 1 0 keep
TD_AsPgAus 8 7 1 ResAsPgAus 0 3
TAdm_NeaEu 6 3 adNeaEu 1 0 keep
TBotEu 6 6 0 ResBotEu 0 4
TD_OoAEu 6 5 1 1 0 4
TBotOcAs 7 7 0 ResBotOcAs 0 4
TD_OoAOcAs 7 5 1 1 0 4
TAdm_NeaAnc 5 3 adNeast 1 0 4
TendstBot 5 5 0 ResstBot 0 4
TstBot 5 5 0 ResEndstBot 0 4
4386 5 4 1 ResMH 0 4
TDeni$ 0 1 1 1 0 4
TNea$ 2 3 1 1 0 keep
TArch$ 1 3 1 1 0 4
TD_HA 4 3 1 ResAnc 0 4
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
FREQ 1 0 1.25e-8 OUTEXP
