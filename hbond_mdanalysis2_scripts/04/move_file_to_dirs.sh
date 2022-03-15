#!/bin/bash
#
#author: Hailey Bureau 
#latest edits: 19 May 2014
#
for dir in $(find . -type d) ; do cp hbpkl_post.py  "$dir" ; done
for dir in $(find . -type d) ; do cp hbpkl_job.sh  "$dir" ; done
