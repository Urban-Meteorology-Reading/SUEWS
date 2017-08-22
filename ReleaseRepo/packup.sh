#!/usr/bin/env bash
# a bash script for packing-up release archive.
# author: Ting Sun, ting.sun@reading.ac.uk
# version history:
# 20170803 - initial version

# clean current workspace if any leftover
rm *.zip
# get the version from Executable
fl=($(ls -d ../SUEWS_V*))
# echo "${exe[@]}"
# echo ${#fl[@]}
exe=$(basename "${fl[0]}")
echo "$exe will be packed up"
program_ver=(${exe//_V/ })
# program="${program_ver[0]}"
version="${program_ver[1]}"
# make bundle name
dir_win="Win10x64"
dir_mac="macOS"
dir_linux="Ubuntu14.04LTS"
zip_win="${exe}_${dir_win}.zip"
# echo $zip_win
zip_mac="${exe}_${dir_mac}.zip"
# echo $zip_mac
zip_linux="${exe}_${dir_linux}.zip"

# Win10x64 bundle
zip -rq ${zip_win} Input Output Manual -x .DS_*
zip -urjq ${zip_win} Executable/${dir_win}/${exe}.exe RunControl.nml
echo "${zip_win} done."
# macOS bundle
zip -rq ${zip_mac} Input Output Manual -x .DS_*
zip -urjq ${zip_mac} Executable/${dir_mac}/${exe} RunControl.nml
echo "${zip_mac} done."
# linux bundle
zip -rq ${zip_linux} Input Output Manual -x .DS_*
zip -urjq ${zip_linux} Executable/${zip_linux}/${exe} RunControl.nml
echo "${zip_mac} done."


# make dir of version, omit if exist
mkdir -p Archive/${version}

# move zip bundles to archive, overwrite if exist
mv ${zip_win} Archive/${version}/.
mv ${zip_mac} Archive/${version}/.
mv ${zip_linux} Archive/${version}/.
echo "Archives moved to Archive/${version}."
