#!/bin/bash

#set -e
#set -u

source ./configs.conf

bash spacewrapper.sh $A $A_image
bash spacewrapper.sh $B $B_image
bash spacewrapper.sh $C $C_image
bash spacewrapper.sh $D $D_image
bash spacewrapper.sh $E $E_image
bash spacewrapper.sh $F $F_image
bash spacewrapper.sh $G $G_image
bash spacewrapper.sh $H $H_image
