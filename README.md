# dnaseq_wdl

This repo follows the progression of transfering the GenPipes' DNASeq mugqic protocol into WDL to use on cloud infrastructure.


A WDL scheduler has been added to GenPipes on the *wdl_scheduler* branch


## Tests on Beluga:
A test environment has been set up on Beluga at:  
`/home/rdali/projects/rrg-bourqueg-ad/rdali/C3G/projects/dnaseqTestSet`

## Create wdl script through GenPipes Scheduler:

`/home/rdali/projects/rrg-bourqueg-ad/rdali/C3G/repos/genpipes//pipelines/dnaseq/dnaseq.py \
-c /home/rdali/projects/rrg-bourqueg-ad/rdali/C3G/repos/genpipes//pipelines/dnaseq/dnaseq.base.ini /home/rdali/projects/rrg-bourqueg-ad/rdali/C3G/repos/genpipes//pipelines/dnaseq/dnaseq.beluga.ini \
-r readset.dnaseq.mini.txt \
-j wdl
-f > test.wdl`


### need interactive node & Java:
`salloc -N 1 -n 1 --mem 10Gb --account=rrg-bourqueg-ad --time=2:00:00`
`module load mugqic/java`

## make json input from GenPipes readset file:

`python readset2json.py -r readset.dnaseqFQBAM.txt -n DnaSeq`

## merge the parameter json with the readset json:

`jq -s '.[0] * .[1]' dnaseq_genpipes.json readset.dnaseqFQBAM.json > dnaseq_genpipes.dnaseqFQBAM.json`

## validate wdl script:
`java -jar $cromwell run dnaseq_genpipes.wdl --inputs dnaseq_genpipes.dnaseqFQBAM.json `

## run command:
`java -jar $cromwell run test.wdl`

Where:  
wdltool=/home/rdali/tools/womtool-51.jar  
cromwell=/home/rdali/tools/cromwell-51.jar  

