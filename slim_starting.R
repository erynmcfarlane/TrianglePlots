###following https:##rdinnager.github.io/slimr/index.html and using old code from McFarlane et al. 2021 Mol Ecol
##if (!require("devtools")) install.packages(devtools) - need devtools
##devtools::install_github("rdinnager/slimr") - install slimr
#if (!require(adegenet)) install.packages("adegenet")
library(slimr)
library(adegenet)
slimr::slim_is_avail() ## check you have slim R

### we want to make a 9 x 9 matrix of populations with which to test traingle plots vs heterozygosity - need to pump out vcf files with which to put into entropy


slim_script(
  slim_block(initialize(),
             {###setting chromosome length
               defineConstant("L", 1e6);
               ###make sure I'm getting nucleotides
               initializeSLiMOptions(nucleotideBased=T);
               ###get the ancestral nucleotides to statsrt with
               initializeAncestralNucleotides(randomNucleotides(L));
               ###recombintation rate uniform across the chromosomes
               initializeRecombinationRate(1e-8);
               ## set the overall mutation rate
               initializeMutationTypeNuc("m1", 0.5, "f", 0.0); 
               ## g1 genomic element type: uses m1 for all mutations
               initializeGenomicElementType("g1", m1, 1.0, mmJukesCantor(1e-8));
               ## uniform chromosome of length 100 kb
               initializeGenomicElement(g1, 0, L-1);
            }),
  slim_block(1,
             {
               sim.addSubpop("p1", 500);
               sim.addSubpop("p2", 500);
              
             }),
  slim_block(2000, late(), 
             {
             p1.setMigrationRates(p2, 0.2);### this is where you can change the rate of migration, can also change the magnitude
             p2.setMigrationRates(p1, 0.2);
            #calcFST(p1.genomes, p2.genomes); ###not working
  }),
  slim_block(2100, late(), ### how long does the simulation run? This can also be divergence time
             {
               g=c(p1.individuals, p2.individuals);
               g.genomes.outputVCF(filePath="test.vcf", simplifyNucleotides=T);
               sim.simulationFinished();
             }))->slim_migration0.2_100gens

sr<-slim_run(slim_migration0.2_100gens, capture_output = "log", keep_all_output = TRUE)


###after this, I want to pull Fst
### I want to compare all of the Fsts, and I want to make sure that we're always comparing Fst to Fst
### Do I want to vary Fst?