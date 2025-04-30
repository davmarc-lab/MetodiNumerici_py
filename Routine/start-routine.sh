#/bin/bash

# define directory name
name="Routine_`date +%F_%H-%M`"

if [ ! -d Sessions ] ; then
    mkdir Sessions
fi

# create the new directory
if [ ! -d Sessions/$name ] ; then
    mkdir Sessions/$name 
fi

# copy the template in the new dir
cp ./RoutineTemplate.ipynb ./Sessions/$name/Routine.ipynb
