#######################################################################################
# definition for your grid engine
#
#   MAX_MEMORY (in GB)
#     is used only at the step of combining orderings.
#      (including bootstrapped values)
#
#     It depends on the number of SNPs (but not on the number of individuals).
#       e.g., 190,000 SNPs requires about 2.6 GB 
#     
#######################################################################################

QUEUE_TYPE="SGE"
QSUB_COMMON="qsub -cwd "
MAX_MEMORY="-l s_vmem=8G -l mem_req=8"

#QUEUE_TYPE="LSF"
#QSUB_COMMON="bsub "
#MAX_MEMORY="-M 8000000"


#######################################################################################
# common functions
#######################################################################################

echo_fail(){
  echo "`date '+%Y%m%d_%H%S'` $*" >&2
  exit 1
}

