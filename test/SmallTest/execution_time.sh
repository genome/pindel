#!/bin/bash
#
# Test if the command executes within the given execution time
#

### Functions ###
function exitWithMessage() {
	echo "$*"
	exit
}

function onExit() {
	rm "${logfile}"
}

### Main code ###

logfile=$(mktemp /tmp/execution_time.XXXXXX)

# clean-up on exit
trap onExit EXIT

# Check arguments
test $# -gt 1 || exitWithMessage "$0 time command [parameters]"
maximum_time=$1
shift

(time -p $*) 2> ${logfile}
actual_time=$(($(grep user "${logfile}" | cut -c 5- | cut -d. -f1)+1))

echo -n "Time spent $(grep user ${logfile} | cut -c 5-)" >&2
test ${actual_time} -le ${maximum_time} || exitWithMessage " was longer than maximum time (${maximum_time})." >&2

echo "." >&2
exit 0

