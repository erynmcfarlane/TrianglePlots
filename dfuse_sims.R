#### everyone hated that I was trying to do this in slim, so we're going to try to do this using dfuse instead.

###NB, this isn't R code, I'm just writing it hear for git purposes

### all new dfuse simulations to look at how assymetry affects consisency in SNPs that are or aren't under selection ####
### all of these are 3 deme cases, rather than 11 deme cases


###One deme, no selection, migration 0.01
/gscratch/buerkle/data/incompatible/dfuse_src/dfuse -d One_deme.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_1.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1.txt -r 20 -g 100 -m 0.01 -c 0 -o onedeme_m0.01_c0_dmi

###One deme, no selection, migration 0.2
/gscratch/buerkle/data/incompatible/dfuse_src/dfuse -d One_deme.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_1.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1.txt -r 20 -g 100 -m 0.2 -c 0 -o onedeme_m0.02_c0_dmi
