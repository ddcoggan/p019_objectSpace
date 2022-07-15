
subject=F117
export SUBJECTS_DIR=/mnt/HDD12TB/freesurfer/
cd /mnt/HDD12TB/freesurfer/subjects
freeview -v \
$subject/mri/T1.mgz \
$subject/mri/wm.mgz \
$subject/mri/brainmask.mgz \
$subject/mri/aseg.presurf.mgz:colormap=lut:opacity=0.2 \
-f $subject/surf/lh.white:edgecolor=blue \
$subject/surf/lh.pial:edgecolor=red \
$subject/surf/rh.white:edgecolor=blue \
$subject/surf/rh.pial:edgecolor=red