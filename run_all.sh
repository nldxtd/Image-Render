ROUNDS=2000
PHOTONS=1000000

mkdir -p output
time bin/sppm testcases/scene_a.txt output/scene_a sppm $ROUNDS $PHOTONS $ROUNDS
