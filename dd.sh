jedi_build=/lfs/h2/emc/da/noscrub/Ting.Lei/dr-jedi-bundle-dev/jedi-bundle/build
rdas_build=/lfs/h2/emc/da/noscrub/Ting.Lei/dr-emc-rdasapp/RDASApp/build

#find "$jedi_build" "$rdas_build" \
#  -path '*vader*' -name flags.make -print \
#  -exec cat {} \;
find "$jedi_build" "$rdas_build" \
  -path '*fv3-jedi*' -name flags.make -print \
  -exec cat {} \;
