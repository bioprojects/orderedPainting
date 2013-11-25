#! /bin/bash

echo_fail(){
  echo "`date '+%Y%m%d_%H%M%S'` $*"
  exit 1
}


LIB_DIR=lib
LOG_DIR=log
VERSION=0.0.4

if [ ! -d "${LIB_DIR}" ]; then
  mkdir ${LIB_DIR}
fi

if [ ! -d "${LOG_DIR}" ]; then
  mkdir ${LOG_DIR}
fi

CMD="chmod 755 ${LIB_DIR}/*.sh"
echo ${CMD}
eval ${CMD}

CMD="chmod 755 ${LIB_DIR}/*.pl"
echo ${CMD}
eval ${CMD}

#
# 1. ChromoPainter
#
if [ ! -f "chromopainter-${VERSION}/chromopainter" ]; then
  URL=http://www.maths.bris.ac.uk/%7Emadjl/finestructure/chromopainter-${VERSION}.tar.gz
  wget ${URL}
  if [ $? -ne 0 ]; then
    echo_fail "Download failed.  Please check the URL: ${URL}"
  fi

  tar xvfz chromopainter-${VERSION}.tar.gz
  cd chromopainter-${VERSION}

#  perl -i -pe 's/gzopen/fopen/g' ChromoPainterMutEM.c
#  perl -i -pe 's/gzprintf/fprintf/g' ChromoPainterMutEM.c
#  perl -i -pe 's/gzclose/fclose/g' ChromoPainterMutEM.c
#  perl -i -pe 's/copyprobsperlocus\.out\.gz/copyprobsperlocus.out/g' ChromoPainterMutEM.c
  
  CC=gcc
  ./configure
  make
  if [ $? -ne 0 ]; then
    echo_fail "Make of chromopainter failed"
  fi
  
  CMD="cp -f ./chromopainter ../${LIB_DIR}"
  echo ${CMD}
  eval ${CMD}
  cd ..
else
  echo "chromopainter-${VERSION}/chromopainter already exists.  Skipped."
fi

#
# 2. official scripts for ChromoPainter
#
arr_scripts=(
  makeuniformrecfile.pl
  neaverage.pl
)

for aa in ${arr_scripts[@]}
do
  if [ ! -f "${LIB_DIR}/${aa}" ]; then
    URL=http://www2.maths.bris.ac.uk/~madjl/finestructure/${aa}.zip
    wget ${URL}
    if [ $? -ne 0 ]; then
      echo_fail "Download failed.  Please check the URL: ${URL}"
    fi

    unzip ${aa}.zip
    CMD="mv -f ${aa} ./${LIB_DIR}"
    echo ${CMD}
    eval ${CMD}
  else
    echo "${aa} already exists.  Skipped."
  fi
done


#
# 3. GNU sort 8.11
#
if [ ! -d "coreutils-8.11" ]; then
  URL=http://ftp.gnu.org/gnu/coreutils/coreutils-8.11.tar.gz
  wget ${URL}
  if [ $? -ne 0 ]; then
    echo_fail "Download failed.  Please check the URL: ${URL}"
  fi

  tar xvfz coreutils-8.11.tar.gz
  cd coreutils-8.11
  CC=gcc
  ./configure
  cd ./lib
  make
  cd ../src
  make version.h
  make sort
  if [ $? -ne 0 ]; then
    echo_fail "Make of GNU sort failed"
  fi

  CMD="cp sort ../../${LIB_DIR}/"
  echo ${CMD}
  eval ${CMD}
  
  cd ../../
else
  echo "coreutils-8.11 already exists.  Skipped."
fi


#
# 4. randomize
#
cd ./${LIB_DIR}/randomize
MSG=`make`
echo ${MSG}
if [ $? -ne 0 ]; then
  echo_fail "Make of randomize failed"
fi

cd ../..


#
# 5. splitN
#
cd ./${LIB_DIR}/splitN
MSG=`make`
echo ${MSG}
if [ $? -ne 0 ]; then
  echo_fail "Make of splitN failed"
fi

cd ../..


#
# 6. postprocess
#
cd ./${LIB_DIR}/postprocess
MSG=`make`
echo ${MSG}
if [ $? -ne 0 ]; then
  echo_fail "Make of postprocess failed"
fi

cd ../..


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

