#/bin/bash

# define directory name
name="Routine_`date +%F-%H:%M`"

# create the new directory
if [ ! -d Sessions/$name ] ; then
    mkdir Sessions/$name 
fi

# copy the template in the new dir
cp ./RoutineTemplate.ipynb ./Sessions/$name
