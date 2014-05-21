#######################################################################################
# definition for your grid engine
#######################################################################################

QUEUE_TYPE="SGE"
QSUB_COMMON="qsub -cwd "

#QUEUE_TYPE="LSF"
#QSUB_COMMON="bsub "

#
# If you need to increase an upper limit of memory usage,
# please specify a correct option in your computational cluster.
#
# The followings are examples.
#
# The memory usage of orderedPainting depends on only the number of SNPs.
#
# This option is used only at the step of combining all orderings.
#
MAX_MEMORY=""

## for SGE, but please check options in your cluster
#MAX_MEMORY="-l s_vmem=8G -l mem_req=8" 

## for LSF
#MAX_MEMORY="-M 8000000"


#######################################################################################
# common functions
#######################################################################################

echo_fail(){
  echo "`date '+%Y%m%d_%H%S'` $*" >&2
  exit 1
}

