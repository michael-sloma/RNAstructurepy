all:
	@echo 'RNAstructure repository testing started...'
	@echo
	make AllSub;
	make bifold;
	make bifold-smp;
	make bipartition;
	make bipartition-smp;
	make CircleCompare;
	make ct2dot;
	make dot2ct;
	make draw;
	make DuplexFold;
	make dynalign;
	make dynalign-smp;
	make dynalign_ii;
	make dynalign_ii-smp;
	make DynalignDotPlot;
	make efn2;
	make efn2-smp;
	make EnergyPlot;
	make EnsembleEnergy;
	make Fold;
	make Fold-smp;
	make MaxExpect;
	make multilign;
	make multilign-smp;
#	make NAPSS; #disable during NAPSS enhancements
	make oligoscreen;
	make oligoscreen-smp;
	make oligowalk;
	make PARTS;
	make partition;
	make partition-smp;
	make ProbabilityPlot;
	make ProbablePair;
	make ProbKnot;
	make refold;
	make RemovePseudoknots;
	make scorer;
	make ShapeKnots;
	make stochastic;
	make stochastic-smp;
	make TurboFold;
	make TurboFold-smp;
	@echo 'RNAstructure repository testing finished.'

cuda:
	@echo 'RNAstructure cuda testing started...'
	@echo
	make fold-cuda;
	make partition-cuda;
	@echo 'cuda testing finished.'

AllSub: tests/AllSub
tests/AllSub:
	@echo 'AllSub testing started...'
	@./testFiles/NixScript AllSub
	@echo 'AllSub testing finished.'
	@echo

bifold: tests/bifold
tests/bifold:
	@echo 'bifold testing started...'
	@./testFiles/NixScript bifold
	@echo 'bifold testing finished.'
	@echo

bifold-smp: tests/bifold-smp
tests/bifold-smp:
	@echo 'bifold-smp testing started...'
	@./testFiles/NixScript bifold
	@echo 'bifold-smp testing finished.'
	@echo


bipartition: tests/bipartition
tests/bipartition:
	@echo 'bipartition testing started...'
	@./testFiles/NixScript bipartition
	@echo 'bipartition testing finished.'
	@echo

bipartition-smp: tests/bipartition-smp
tests/bipartition-smp:
	@echo 'bipartition-smp testing started...'
	@./testFiles/NixScript bipartition
	@echo 'bipartition-smp testing finished.'
	@echo


CircleCompare: tests/CircleCompare
tests/CircleCompare:
	@echo 'CircleCompare testing started...'
	@./testFiles/NixScript CircleCompare
	@echo 'CircleCompare testing finished.'
	@echo

ct2dot: tests/ct2dot
tests/ct2dot:
	@echo 'ct2dot testing started...'
	@./testFiles/NixScript ct2dot
	@echo 'ct2dot testing finished.'
	@echo

dot2ct: tests/dot2ct
tests/dot2ct:
	@echo 'dot2ct testing started...'
	@./testFiles/NixScript dot2ct
	@echo 'dot2ct testing finished.'
	@echo

draw: tests/draw
tests/draw:
	@echo 'draw testing started...'
	@./testFiles/NixScript draw
	@echo 'draw testing finished.'
	@echo

DuplexFold: tests/DuplexFold
tests/DuplexFold:
	@echo 'DuplexFold testing started...'
	@./testFiles/NixScript DuplexFold
	@echo 'DuplexFold testing finished.'
	@echo

dynalign: tests/dynalign
tests/dynalign:
	@echo 'dynalign testing started...'
	@./testFiles/NixScript dynalign
	@echo 'dynalign testing finished.'
	@echo

dynalign-smp: tests/dynalign-smp
tests/dynalign-smp:
	@echo 'dynalign-smp testing started...'
	@./testFiles/NixScript dynalign-smp
	@echo 'dynalign-smp testing finished.'
	@echo


dynalign_ii: tests/dynalign_ii
tests/dynalign_ii:
	@echo 'dynalign_ii testing started...'
	@./testFiles/NixScript dynalign_ii
	@echo 'dynalign_ii testing finished.'
	@echo

dynalign_ii-smp: tests/dynalign_ii-smp
tests/dynalign_ii-smp:
	@echo 'dynalign_ii-smp testing started...'
	@./testFiles/NixScript dynalign_ii-smp
	@echo 'dynalign_ii-smp testing finished.'
	@echo


DynalignDotPlot: tests/DynalignDotPlot
tests/DynalignDotPlot:
	@echo 'DynalignDotPlot testing started...'
	@./testFiles/NixScript DynalignDotPlot
	@echo 'DynalignDotPlot testing finished.'
	@echo

efn2: tests/efn2
tests/efn2:
	@echo 'efn2 testing started...'
	@./testFiles/NixScript efn2
	@echo 'efn2 testing finished.'
	@echo

efn2-smp: tests/efn2-smp
tests/efn2-smp:
	@echo 'efn2-smp testing started...'
	@./testFiles/NixScript efn2
	@echo 'efn2-smp testing finished.'
	@echo


EnergyPlot: tests/EnergyPlot
tests/EnergyPlot:
	@echo 'EnergyPlot testing started...'
	@./testFiles/NixScript EnergyPlot
	@echo 'EnergyPlot testing finished.'
	@echo

EnsembleEnergy: tests/EnsembleEnergy
tests/EnsembleEnergy:
	@echo 'EnsembleEnergy testing started...'
	@./testFiles/NixScript EnsembleEnergy
	@echo 'EnsembleEnergy testing finished.'
	@echo

Fold: tests/Fold
tests/Fold:
	@echo 'Fold testing started...'
	@./testFiles/NixScript Fold
	@echo 'Fold testing finished.'
	@echo


Fold-smp: tests/Fold-smp
tests/Fold-smp:
	@echo 'Fold-smp testing started...'
	@./testFiles/NixScript Fold
	@echo 'Fold-smp testing finished.'
	@echo


fold-cuda: tests/Fold-cuda
tests/Fold-cuda:
	@echo 'fold-cuda testing started...'
	@./testFiles/NixScript fold-cuda
	@echo 'fold-cuda testing finished.'
	@echo

MaxExpect: tests/MaxExpect
tests/MaxExpect:
	@echo 'MaxExpect testing started...'
	@./testFiles/NixScript MaxExpect
	@echo 'MaxExpect testing finished.'
	@echo

multilign: tests/multilign
tests/multilign:
	@echo 'multilign testing started...'
	@./testFiles/NixScript multilign
	@echo 'multilign testing finished.'
	@echo

multilign-smp: tests/multilign-smp
tests/multilign-smp:
	@echo 'multilign-smp testing started...'
	@./testFiles/NixScript multilign-smp
	@echo 'multilign-smp testing finished.'
	@echo

NAPSS: tests/NAPSS
tests/NAPSS:
	@echo 'NAPSS testing started...'
	@./testFiles/NixScript NAPSS
	@echo 'NAPSS testing finished.'
	@echo

oligoscreen: tests/oligoscreen
tests/oligoscreen:
	@echo 'oligoscreen testing started...'
	@./testFiles/NixScript oligoscreen
	@echo 'oligoscreen testing finished.'
	@echo

oligoscreen-smp: tests/oligoscreen-smp
tests/oligoscreen-smp:
	@echo 'oligoscreen-smp testing started...'
	@./testFiles/NixScript oligoscreen
	@echo 'oligoscreen-smp testing finished.'
	@echo


oligowalk: tests/oligowalk
tests/oligowalk:
	@echo 'oligowalk testing started...'
	@./testFiles/NixScript oligowalk
	@echo 'oligowalk testing finished.'
	@echo

PARTS: tests/PARTS
tests/PARTS:
	@echo 'PARTS testing started...'
	@./testFiles/NixScript PARTS
	@echo 'PARTS testing finished.'
	@echo

partition: tests/partition
tests/partition:
	@echo 'partition testing started...'
	@./testFiles/NixScript partition
	@echo 'partition testing finished.'
	@echo

partition-smp: tests/partition-smp
tests/partition-smp:
	@echo 'partition-smp testing started...'
	@./testFiles/NixScript partition
	@echo 'partition-smp testing finished.'
	@echo


partition-cuda: tests/partition-cuda
tests/partition-cuda:
	@echo 'partition-cuda testing started...'
	@./testFiles/NixScript partition-cuda
	@echo 'partition-cuda testing finished.'
	@echo

ProbabilityPlot: tests/ProbabilityPlot
tests/ProbabilityPlot:
	@echo 'ProbabilityPlot testing started...'
	@./testFiles/NixScript ProbabilityPlot
	@echo 'ProbabilityPlot testing finished.'
	@echo

ProbKnot: tests/ProbKnot
tests/ProbKnot:
	@echo 'ProbKnot testing started...'
	@./testFiles/NixScript ProbKnot
	@echo 'ProbKnot testing finished.'
	@echo

ProbablePair: tests/ProbablePair
tests/ProbablePair:
	@echo 'ProbablePair testing started...'
	@./testFiles/NixScript ProbablePair
	@echo 'ProbablePair testing finished.'
	@echo

refold: tests/refold
tests/refold:
	@echo 'refold testing started...'
	@./testFiles/NixScript refold
	@echo 'refold testing finished.'
	@echo

RemovePseudoknots: tests/RemovePseudoknots
tests/RemovePseudoknots:
	@echo 'RemovePseudoknots testing started...'
	@./testFiles/NixScript RemovePseudoknots
	@echo 'RemovePseudoknots testing finished.'
	@echo

scorer: tests/scorer
tests/scorer:
	@echo 'scorer testing started...'
	@./testFiles/NixScript scorer
	@echo 'scorer testing finished.'
	@echo

ShapeKnots: tests/ShapeKnots
tests/ShapeKnots:
	@echo 'ShapeKnots testing started...'
	@./testFiles/NixScript ShapeKnots
	@echo 'ShapeKnots testing finished.'
	@echo

stochastic: tests/stochastic
tests/stochastic:
	@echo 'stochastic testing started...'
	@./testFiles/NixScript stochastic
	@echo 'stochastic testing finished.'
	@echo

stochastic-smp: tests/stochastic-smp
tests/stochastic-smp:
	@echo 'stochastic-smp testing started...'
	@./testFiles/NixScript stochastic
	@echo 'stochastic-smp testing finished.'
	@echo


TurboFold: tests/TurboFold
tests/TurboFold:
	@echo 'TurboFold testing started...'
	@./testFiles/NixScript TurboFold
	@echo 'TurboFold testing finished.'
	@echo

TurboFold-smp: tests/TurboFold-smp
tests/TurboFold-smp:
	@echo 'TurboFold-smp testing started...'
	@./testFiles/NixScript TurboFold-smp
	@echo 'TurboFold-smp testing finished.'
	@echo

runTests: all
testing: all
tests: all
