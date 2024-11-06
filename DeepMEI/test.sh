#! /bin/bash
echo $0
base=`echo $0 |perl -npe "s/\.\/test.sh$//"|perl -F'\t' -alne 'use Cwd(getcwd,cwd);if($_=~/^\//) {print $_;} else { print getcwd()."\/".$_;}'`
echo $base
