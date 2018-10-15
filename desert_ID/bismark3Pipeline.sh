#!/bin/sh 

# Wrapper script for submitting slurm commands
# 
#
# Usage: bash bismark2Pipeline.sh [OPTIONS] [DIRECTORY]
# Submit slurm jobs for bismark pipeline
# 
# -r <# reads>    total number of reads needed for fastqSplitter.slurm
#                 [protip: get # reads from fastqc]
# -C              do not scan directory for fastq files (good if running only portion of pipeline)
# -p              submit to dtg pipeline (-p dtg), excludes bismark alignment
# -P              submit bismark alignment to dtg, using only 1 node (-p dtg -N 1)
# -X              submit bismark alignment to dtg, using 12 nodes [WARNING: use with care!]
# -q              run bamqc
# -Q              run bamqc (deduped version)

GENOME="hg19"
# Get flags
while getopts :abc:r:pPXqQC OPT; do
	case "$OPT" in
	a) echo "-a";;
	b) echo "-b";;
	c) NCUTS=${OPTARG};;
    r) NREADS=${OPTARG};;
	p) PTT="-p dtg";;
    P) PTT2="-p dtg -n 1 -N 1 -c 16";;
	X) PTT2="-p dtg";;
	q) DOBAMQC="true";;
	Q) DOBAMQC="dedup";;
	C) CHECKFQ="false";;
	esac
done
shift $(($OPTIND-1))
DIR=${1:-"."}
cd $DIR
WD=`pwd`

# Default variables
PTT=${PTT:-""}

echo "Starting bismark parallel..."

# Check if fastqs exist (disabled with -C)
if [[ ! $CHECKFQ = "false" ]]; then
	echo "Looking for fastq files in $DIR..."
	# Check for R1 and R2 files 
	R1=($(ls ${DIR}/*R1*fastq*))
	R2=($(ls ${DIR}/*R2*fastq*))
	[[ ${#R1[@]} == 1 ]] || echo "More than one R1 file found!"
	[[ ${#R1[@]} == 1 ]] || exit 0
	[[ ${#R2[@]} == 1 ]] || echo "More than one R2 file found!"
	[[ ${#R2[@]} == 1 ]] || exit 0
	echo "Found R1 = $R1"
	echo "Found R2 = $R2"
fi

# Format $NCUTS (number of cuts to split fastqs)
if [[ ! -z "$NCUTS" ]] & [[ "$NCUTS" -gt 0 ]]; then
	NCUTP1=$(( $NCUTS + 1 ))
	echo "Splitting fastqs into $NCUTS bits and requesting $NCUTP1 nodes for alignment..."
	PTT2="-n $NCUTP1 -N $NCUTP1 -c 16"
fi

# Create folder for log files
test -e logs || mkdir logs
LOGPATH="$WD/logs"

# Location of oarallel scripts...
BISPPATH=/home/dtg-00/Groups/Salhia_lab/bin/pipelines/bismark_parallel/

# Location of alignment scripts
FASTQSPLITTER="$BISPPATH/fastqSplitter.slurm"
FASTQTRIMALL="$BISPPATH/fastqtrimall.slurm"
BISMARKBATCH="$BISPPATH/bismarkBatch.slurm"
CATALLBAM="$BISPPATH/catallbam.slurm"
SAMBLASTERDEDUP="$BISPPATH/samblaster.slurm"
BISMARKMETHYL="$BISPPATH/bismarkCallMethyl.slurm"
BISMARKCOV2CYTOSINE="$BISPPATH/cov2cytosine.slurm"
BISMARKBIGWIGS="$BISPPATH/bismarkCov2BigWigs.slurm"

# Function for extracting job ID for dependency
format_jobid () {
	local JID=`echo $@ | cut -d " " -f4`
	echo "-d $JID"
}

# Actual bismark pipeline
STEPFQSPLIT="sbatch -D $WD -o $LOGPATH/fastqsplit.out $DEPEND $PTT $FASTQSPLITTER -r $DIR $NREADS $NCUTS"
echo "$STEPFQSPLIT"
JOBID=`$STEPFQSPLIT`
DEPEND=`format_jobid $JOBID`

STEPFQTRIM="sbatch -D $WD -o $LOGPATH/fastqtrim.out $DEPEND $PTT $FASTQTRIMALL -n $DIR"
echo "$STEPFQTRIM"
JOBID=`$STEPFQTRIM`
DEPEND=`format_jobid $JOBID`

STEPBISMARK="sbatch -D $WD -o $LOGPATH/bismarkalign.out $DEPEND $PTT2 $BISMARKBATCH $DIR"
echo "$STEPBISMARK"
JOBID=`$STEPBISMARK`
DEPEND=`format_jobid $JOBID`

STEPCATALLBAM="sbatch -D $WD -o $LOGPATH/catallbam.out $DEPEND $PTT $CATALLBAM $DIR"
echo "$STEPCATALLBAM"
JOBID=`$STEPCATALLBAM`
DEPEND=`format_jobid $JOBID`

STEPSAMBLASTER="sbatch -D $WD/merged $DEPEND $PTT $SAMBLASTERDEDUP -r -d $DIR/merged"
echo "$STEPSAMBLASTER"
JOBID=`$STEPSAMBLASTER`
DEPEND=`format_jobid $JOBID`

STEPBISMARKMETH="sbatch -D $WD/merged -o $LOGPATH/bismarkmethyl.out $DEPEND $PTT $BISMARKMETHYL $DIR/merged"
echo "$STEPBISMARKMETH"
JOBID=`$STEPBISMARKMETH`
DEPEND=`format_jobid $JOBID`

STEPBISMARKCOV2CYTO="sbatch -D $WD/merged -o $LOGPATH/cov2cytosine.out $DEPEND $PTT $BISMARKCOV2CYTOSINE $DIR/merged"
echo "$STEPBISMARKCOV2CYTO"
JOBID=`$STEPBISMARKCOV2CYTO`
DEPEND=`format_jobid $JOBID`

STEPBISMARKBIGWIGS="sbatch -D $WD/merged -o $LOGPATH/bismarkbigwigs_m.out $DEPEND $PTT $BISMARKBIGWIGS -m $DIR/merged"
echo "$STEPBISMARKBIGWIGS"
JOBID=`$STEPBISMARKBIGWIGS`
DEPEND=`format_jobid $JOBID`

STEPBISMARKBIGWIGS2="sbatch -D $WD/merged -o $LOGPATH/bismarkbigwigs_u.out $DEPEND $PTT $BISMARKBIGWIGS -u $DIR/merged"
echo "$STEPBISMARKBIGWIGS2"
JOBID=`$STEPBISMARKBIGWIGS2`
DEPEND=`format_jobid $JOBID`

### Start bamqc branch
# Location of bamqc scripts
CATMDUPSORT="$BISPPATH/catmdupsort.slurm"
QUALIMAP="$BISPPATH/qualimapDtg.slurm"

if [[ "$DOBAMQC" = "true" ]]; then
	STEPCATMDUPSORT="sbatch -D $WD -o $LOGPATH/catmdupsort.out $DEPEND $PTT $CATMDUPSORT $DIR"
	echo "$STEPCATMDUPSORT"
	JOBID=`$STEPCATMDUPSORT`
	DEPEND=`format_jobid $JOBID`

	STEPQUALIMAP="sbatch -D $WD/sorted -o $LOGPATH/qualimap.out $DEPEND $PTT $QUALIMAP $DIR/sorted"
	echo "$STEPQUALIMAP"
	JOBID=`$STEPQUALIMAP`
	DEPEND=`format_jobid $JOBID`
fi

### bamqc branch #2 (deduped)
# Location of bamqc scripts
CALLSAMSORT="$BISPPATH/callsamsort.slurm"
QUALIMAP="$BISPPATH/qualimapDtg.slurm"

if [[ "$DOBAMQC" = "dedup" ]]; then
	STEPSAMSORT2="sbatch -D $WD/merged -o $LOGPATH/samsort2.out $DEPEND $PTT $CALLSAMSORT $DIR/merged"
	echo "$STEPSAMSORT2"
	JOBID=`$STEPSAMSORT2`
	DEPEND=`format_jobid $JOBID`

	STEPQUALIMAP="sbatch -D $WD/merged/sorted -o $LOGPATH/qualimap2.out $DEPEND $PTT $QUALIMAP $DIR/merged/sorted"
	echo "$STEPQUALIMAP"
	JOBID=`$STEPQUALIMAP`
	DEPEND=`format_jobid $JOBID`
fi

