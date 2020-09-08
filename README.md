# dnaseq_wdl

This repo follows the progress of transfering the GenPipes' DNASeq mugqic protocol into WDL to use on cloud infrastructure.


A WDL scheduler has been added to GenPipes on the *wdl_scheduler* branch


## Tests on Beluga:
A test environment has been set up on Beluga at:  
`/home/rdali/projects/rrg-bourqueg-ad/rdali/C3G/projects/dnaseqTestSet`



### need interactive node to run Java on Beluga:
`salloc -N 1 -n 1 --mem 10Gb --account=rrg-bourqueg-ad --time=2:00:00`  
`module load mugqic/java`

## make json input from GenPipes readset file:

`python utils/readset2json.py -r readset.dnaseq.full.txt -n DnaSeq`

## merge the parameter json with the readset json:

`jq -s '.[0] * .[1]' dnaseq_genpipes.json myReadsetFile.json > dnaseq_genpipes_readset.dnaseq.full.json `

## validate wdl script:
`java -jar $wdltool validate dnaseq_genpipes.wdl`

## run command through SLURM:

sbatch --account=rrg-bourqueg-ad --time=24:00:00 -N 1 -n 1 --mem 10G -o slurm.cromwellTest_test --wrap="module load mugqic/java/openjdk-jdk1.8.0_72;
java -Dconfig.file=./cromwell.conf -jar $cromwell run dnaseq_genpipes.wdl --inputs dnaseq_genpipes_readset.dnaseq.full.json --options ./options.json"

Where:  
wdltool=/home/rdali/tools/womtool-51.jar  
cromwell=/home/rdali/tools/cromwell-51.jar  

