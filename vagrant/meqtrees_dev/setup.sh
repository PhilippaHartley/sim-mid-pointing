# Set up environment for running MeqTrees.

export MQT=${HOME}/meqtrees
export PATH=${MQT}/tigger/TigGUI/bin:${MQT}/purr/Purr/bin:${MQT}/owlcat/Owlcat/bin:${MQT}/pyxis/Pyxis/bin:${PATH}
export PYTHONPATH=${HOME}/analysis_scripts:${MQT}/owlcat:${MQT}/pyxis:${MQT}/kittens:${MQT}/purr:${MQT}/tigger:${PYTHONPATH}
export TIMBA_PATH=${MQT}/meqtrees-timba
export MEQTREES_CATTERY_PATH=${MQT}/meqtrees-cattery/Cattery
source ${MQT}/meqtrees-timba/install/timba-init.sh
timba-setup release
