#!/bin/bash
#
# Run the casa simulator
#

##### Constants
CASA=/home/myso/opt/casa-release-4.7.2-el7/bin/casa
CASA_OPT="--nologger --nogui --log2term"
SCRIPT="$(dirname "$(readlink -f "$0")")/casa_simulator.py"

##### Functions
function usage {
    cat <<fin
Usage: run_casa_simulator.sh [[-h]] config input output

Arguments:
  config                Configuration file name
  input                 Input image file name
  output                Output image file name

Options:
  -h, --help            Help
fin
}

# Process arguments
options=""
while [ "$1" != "" ]; do
    case $1 in
        -h | --help )           usage
                                exit
                                ;;
        * )                     config=$1
                                shift
                                input=$1
                                shift
                                output=$1
                                shift
                                ;;
    esac
done

pwd >> last_casa_simulator.log
echo "$CASA $CASA_OPT -c $SCRIPT $config $input $output" >> last_casa_simulator.log
$CASA $CASA_OPT -c $SCRIPT $config $input $output

rm -rf CONF*
rm -rf cgrid_ft.im
#rm -rf casa*.log ipython*.log
