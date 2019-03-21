#!/bin/bash

# Name: example.sh
#
# Description:
#
#   It should be possible for someone else- or you when you come back
#   in six months- to learn how to use the script by reading the
#   comments (and self-help, if provided) without reading the code.
#   Be kind to future you.
#
#   "Private" Bash functions are designated by names beginning with
#   "_fn_", and public Bash functions by "fn_".
#
# Example:
#
#   Show an example of how your script is to be invoked, for example:
#   sojinq --user_id=USER NAME --project_id=PROJECT ID
#
# Args:
#
#   If your script takes arguments during invocation, document them here, like so:
#   param1 (int): The first parameter.
#   param2 (str, optional): The secod parameter. Defaults to None.
#   If your script doesn't take any arguments, put "None".
#
# Returns:
#
#   Document any special return codes. If there aren't any, put "None."
#
# Files changed:
#
#   Document any files changed by the script. If there aren't any, put "None."
#
# Authors:
#
#   See AUTHORS in this project's root directory.
#
# Contributors:
#
#   See CONTRIBUTORS in this project's root directory.
#   See CONTRIBUTING.md in this project's root directory for how to contribute.
#
# Changelog:
#
#   See the CHANGELOG.md in this project's root directory.
#
# License:
#
#   See LICENSE.md in this project's root directory.
#

#
# ----------------------------------------------------------------------
# --- Bash Options
# ----------------------------------------------------------------------
#

# No unset vars
set -o nounset

# Exit the script if any statement returns a non-true return value.
# If you need to test the return value of a command you suspect might fail
# use the following construct: 
# ifail || ret=$? || true
set -o errexit

# Produce a failure return code if any command errors.
set -o pipefail

# Enable errtrace so that the ERR trap is also triggered when the error
# (a ccommand returning a nonzero code) occurs inside a function or a
# subshell. The context of a function or a subshell does not inherit the
# ERR trap unless the errtrace is enabled.
set -o errtrace

#
# ----------------------------------------------------------------------
# --- Global Variables
# ----------------------------------------------------------------------
# ---
# --- Conventions for global variables:
# ---
# ---   Keep the use of global variables to a minimum.
# ---   Use UPPER_CASE for naming.
# ---   Use readonly declarations to make global variables immutable.
# ---   Use global variables to replace cryptic $0, $1, and so on.
#

readonly PROGNAME=$(basename "${0}")
readonly PROGDIR=$(readlink -m "$(dirname "${0}")")
readonly BIOREALM_PROJECT_DIR="$(dirname "${PROGDIR}")"
readonly ARGS="${@}"

#
# ----------------------------------------------------------------------
# --- Standard Functions
# ----------------------------------------------------------------------
#

#
# ----------------------------------------------------------------------
# --- Perform any housekeeping required if we receive a SIGHUP, SIGINT,
# --- or SIGTERM signal, or if we are called by the _fn_error_trap.
# --- This will not be performed if we receive a SIGKILL (kill -9)
# --- since SIGKILL immediately terminates our process.
# --- For more information, see: http://linuxcommand.org/wss0160.php
# ----------------------------------------------------------------------
#
_fn_clean_up() {
  local _exit_status
  _exit_status=$1

  # Add any additiona housekeeping code to be executed upon
  # termination here.

  exit ${_exit_status}
}

# -------------------------------------
# --- Display an error message if we catch any exceptions and call the
# --- _fn_clean_up function to perform any housekeeping upon termination.
# -------------------------------------
# Globals:
#   None
# Arguments:
#   $0 - The name of the script as passed via the Bash exception trap.
#   $1 - The last line executed.
#   $2 - The error code returned.
# Returns:
#   None
# -------------------------------------
_fn_error_trap() {
  local _exit_status
  local _script_name="$0"    # equals to my script name
  local _last_line="$1"      # argument 1: last line of error occurence
  local _last_err="$2"       # argument 2: error code of last command
  
  _exit_status=1
  echo "-------------------------------------"
  echo "Caught an exception."
  echo "    File: ${_script_name}"
  echo "    Line number: ${_last_line}"
  echo "    Exit status: ${_last_err}"
  echo "-------------------------------------"
  _fn_clean_up ${_exit_status}
}

# -------------------------------------
# Display a message in the desired color.
# Globals:
#   None
# Arguments:
#   $1 - Desired text color. One of red, yellow, green, or blue.
#   $2 - The message to display.
# Returns:
#   None
# -------------------------------------
fn_msg() {
  local color="${1}"
  local msg="${2}"

  case "${color}" in
    red)    printf "\e[1;31m%s\e[0m" "${msg}";;
    yellow) printf "\e[1;33m%s\e[0m" "${msg}";;
    green)  printf "\e[1;32m%s\e[0m" "${msg}";;
    blue)   printf "\e[1;36m%s\e[0m" "${msg}";;
    *)      printf "\e[1m%s\e[0m" "${msg}";;
  esac
}

# -------------------------------------
# Display a message in the desired color followed by a new line.
# Globals:
#   None
# Arguments:
#   $1 - Desired text color. One of red, yellow, green, or blue.
#   $2 - The message to display.
# Returns:
#   None
# -------------------------------------
fn_msg_w_newline() {
  local color="${1}"
  local msg="${2}"

  fn_msg "${color}" "${msg}"
  printf "\n"
}

# -------------------------------------
# Display a message in green that begins with the name of the script and
# the word "[INFO]"
# Globals:
#   PROGNAME
# Arguments:
#   $1 - The message text to display after "[INFO]"
# Returns:
#   None
# -------------------------------------
fn_log_info() {
  local msg="${1}"

  fn_msg_w_newline green "${PROGNAME}: [INFO] ${msg}"
}

# -------------------------------------
# Display a message in yellow that begins with the name of the script and
# the word "[WARNING]"
# Globals:
#   PROGNAME
# Arguments:
#   $1 - The message text to display after "[WARNING]"
# Returns:
#   None
# -------------------------------------
fn_log_warn() {
  local msg="${1}"

  fn_msg_w_newline yellow "$PROGNAME: [WARNING] ${msg}"
}

# -------------------------------------
# Display a message in red that begins with the name of the script and
# the word "[ERROR]"
# Globals:
#   PROGNAME
# Arguments:
#   $1 - The message text to display after "[ERROR]"
# Returns:
#   None
# -------------------------------------
fn_log_error() {
  local msg="${1}"

  fn_msg_w_newline red "$PROGNAME: [ERROR] ${msg}"
}

# -------------------------------------
# Display a message in red that begins with the name of the script and
# the word "[ERROR]" then terminate the script.
# Globals:
#   PROGNAME
# Arguments:
#   $1 - The message text to display after "[ERROR]"
# Returns:
#   None
# -------------------------------------
fn_log_terminate() {
  local msg="${1}"

  fn_log_error "$PROGNAME: [ERROR] ${msg}"
  fn_log_error "Called with arguments ${*}"
  exit 1
}

# -------------------------------------
# Display a warning message and continue execution.
# Globals:
#   None
# Arguments:
#   $1 - A description of what we were doing at the time.
#   $2 - The command/function we're wrapping.
#   $3 and so on - Arguments to pass to the command/function we're wrapping.
# Returns:
#   None
# -------------------------------------
fn_wrap_and_continue_on_error() {
  local warning_msg="${1}"
  shift
  local cmd="${1}"
  shift

  $cmd "$@" || fn_log_warn "$warning_msg"
}

# -------------------------------------
# Display an error message and abort execution.
# exit 1.
# Globals:
#   None
# Arguments:
#   $1 - A description of what we were doing at the time.
#   $2 - The command/function we're wrapping.
#   $3 and so on - Arguments to pass to the command/function we're wrapping.
# Returns:
#   None
# -------------------------------------
fn_wrap_and_fail_on_error() {
  local error_msg="${1}"
  shift
  local cmd="${1}"
  shift

  $cmd "$@" || fn_log_terminate "$error_msg"
}

#
# ----------------------------------------------------------------------
# --- Custom Functions
# ----------------------------------------------------------------------
#

#
# ----------------------------------------------------------------------
# --- Main
# ----------------------------------------------------------------------
#

main() {

# Local vars
#    Example:
#    local i

# Mainline code goes here, should be almost entirely calls to functions
# defined above.

}

#
# ----------------------------------------------------------------------
# --- ADD NO CODE BELOW THIS LINE
# ----------------------------------------------------------------------
# ---
# --- Only global declarations (at the top), the trap set ups, and
# --- the call to main (below) go in the global scope.
# ---
#-----------------------------------------------------------------------
#

trap "_fn_error_trap ${LINENO} $?" ERR

trap _fn_clean_up SIGHUP SIGINT SIGTERM

main "${ARGS}"
