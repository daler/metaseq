#!/bin/bash


set -e

# Defaults ------------------------------------------------------
MINICONDA_DIR="${HOME}/miniconda"
INSTALL_DIR="${HOME}/tools"
TO_INSTALL="bedtools,samtools,tabix,ucsc"
INSTALL_MINICONDA=1
ENVNAME="metaseq-test"
LOG=/dev/stdout
VERBOSE=0
GIT_TAG=""
TRAVIS_CI=0

USE_BEDTOOLS_VERSION=2.21.0
USE_SAMTOOLS_VERSION=0.1.19
USE_TABIX_VERSION=0.2.6
MINICONDA_VERSION=3.7.0

# ---------------------------------------------------------------

while getopts "hd:i:m:Me:vg:t" x; do
    case "$x" in
        h)
            echo "
usage: $0 [options] | check

    This script helps install all the necessary prerequisites for running
    metaseq.

    When the single argument \"check\" is provided, this script simply prints
    any metaseq prerequisites it finds on the path.  Example usage:

        $0 check

    Otherwise:

    This script can install the prerequisite genomics tools (BEDTools,
    samtools, tabix, and the UCSC tools bigWigSummary, bedGraphToBigWig, and
    bigWigToBedGraph) necessary for running metaseq.  It can also install
    a standalone Python installation (Miniconda, a more streamlined version of
    Anaconda; more info at http://docs.continuum.io/anaconda/index.html) that
    includes metaseq and all its Python requirements.  This can be used
    alongside any other Python versions that may be installed.

    By default, this script will:

        1) Download and install BEDTools, samtools, tabix, bigWigSummary,
           bedGraphToBigWig, bigWigToBedGraph.  It will place them in the
           following directory:

              ${INSTALL_DIR}

        2) Download and install Miniconda into the following directory:

              ${MINICONDA_DIR}

        3) Create an isolated environment, located at:

              ${HOME}/miniconda/envs/${ENVNAME}

        4) Install Python prerequisites for metaseq in this environment

        5) Install the latest metaseq release in this environment

    Further instructions are displayed upon successful installation.

    Options
    -------

    -h           Print this help message and exit.

    -d DIR       Directory where pre-requisistes will be installed.
                 Default is ${INSTALL_DIR}.

    -i STRING    A string (usually comma-separated with no spaces) containing
                 the names of prerequisite tools to install in DIR.  Options
                 are bedtools, samtools, tabix, ucsc. Default is
                 \"${TO_INSTALL}\".

    -m DIR       Directory where miniconda (an isolated python distribution) will be installed.
                 Default is ${MINICONDA_DIR}.

    -e ENVNAME   Name of new environment that will be created with all the
                 Python prerequisites.  Default is ${ENVNAME}.

    -M           Disable installation of Miniconda.  This option is best if you
                 want to use an existing installation of the scientific Python
                 stack (scipy, numpy, pandas, matplotlib). However, even if
                 you already have these installed, a Miniconda installation is
                 completely separate and so -M is not required.  It just saves
                 some setup time.

    -v           Verbose mode.  By default, the stdout and stderr from each
                 installation step will be written to a separate log file.  If
                 -v is used, then all output will be written to stdout.

    -t           Special flag used for running tests on travis-ci.org.
                 Disables interactive prompts and expects to be called from
                 inside the metaseq source directory.

    -g TAG       By default, this script installs the latest metaseq release.
                 Use the -g flag to specify a git tag from the metaseq github
                 repository.  The latest version can be accessed using the
                 'master' tag. Use the special tag 'disable' to prevent metaseq
                 from being installed at all.


    The genomics tool versions can be changed by editing the
    'USE_<tool>_VERSION' variables at the top of this script. The versions
    currently configured are:

        bedtools      : ${USE_BEDTOOLS_VERSION}
        samtools      : ${USE_SAMTOOLS_VERSION}
        tabix         : ${USE_TABIX_VERSION}
        bigWigSummary : latest version on UCSC servers


    Example usage:

        # \"The works\": download and install all genomics tools as well as an
        # isolated environment.  Depending on your hardware and internet
        # connection speed, this could take 5 minutes or more.

            $0

        # Already have bedtools, samtools, tabix, bigWigSummary and just want
        # to create an isolated environment in 'my-env'

            $0 -i \"\" -e my-env


        # Install the genomics tools without making a new isolated environment:

            $0 -M

        # Just want the prerequisites and an isolated environment but don't
        # want to install metaseq?

            $0 -g disable

        # Output *everything* to the terminal instead of saving to separate
        # logs.  This generates lots of output but it can be useful for
        # debugging:

            $0 -v


            "
            exit 2
            ;;
        d)
            INSTALL_DIR="${OPTARG}"
            ;;
        e)
            ENVNAME="${OPTARG}"
            ;;
        i)
            TO_INSTALL="${OPTARG}"
            ;;
        m)
            MINICONDA_DIR="${OPTARG}"
            ;;
        M)
            INSTALL_MINICONDA=0
            ;;
        g)
            GIT_TAG="${OPTARG}"
            ;;
        v)
            VERBOSE=1
            ;;
        t)
            # undocumented option for running tests on travis-ci, which assumes
            # we're in the checked-out metaseq directory
            TRAVIS_CI=1
            ;;
        ?)
            echo "did not recognize option, try -h"
            exit 1
            ;;
    esac
done

# Logs first arg to stdout, with a timestamp
log ()
{
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] : $1"
}


# Mac or Linux?
uname | grep "Darwin" > /dev/null && SYSTEM_TYPE=mac || SYSTEM_TYPE=linux
log "Detected operating system: ${SYSTEM_TYPE} (based on the output of uname, which is \"$(uname)\")".

MACH_TYPE=$(uname -m)

if [[ ($MACH_TYPE != "x86_64") && ($SYSTEM_TYPE = "mac") && ( $INSTALL_MINICONDA = 1 ) ]]; then
    echo "
    Sorry, installing miniconda on 32-bit Mac OSX is not supported.  Please see
    http://conda.pydata.org/miniconda.html.  Exiting!
    "
    exit 1
fi

# Determine how to download files.  wget is (typically) installed by default on
# Linux; curl on Mac.
if [[ $SYSTEM_TYPE = "mac" ]]; then
    downloader ()
    {
        log "Downloading $1"
        curl --location $1 > $2
    }
else
    downloader ()
    {
        log "Downloading $1"
        wget $1 -O $2
    }
fi

# Functions to check installation --------------------------------------------
check_bedtools () {
    command -v bedtools > /dev/null 2>&1 \
        && echo "yes [ $(bedtools --version) ] $(which bedtools)" \
        || echo "no"
}

check_tabix () {
    command -v tabix > /dev/null 2>&1 \
        && echo "yes [ $(tabix 2>&1 > /dev/null | grep Version) ] $(which tabix)" \
        || echo "no"
}

check_samtools () {
    command -v samtools > /dev/null 2>&1 \
        && echo "yes [ $(samtools 2>&1 > /dev/null | grep Version) ] $(which samtools)" \
        || echo "no"
}

check_bigWigSummary () {
    command -v bigWigSummary > /dev/null 2>&1 \
        && echo "yes [ ] $(which bigWigSummary)" \
        || echo "no"
}

check_bigWigToBedGraph () {
    command -v bigWigToBedGraph > /dev/null 2>&1 \
        && echo "yes [ ] $(which bigWigToBedGraph)" \
        || echo "no"
}

check_bedGraphToBigWig () {
    command -v bedGraphToBigWig > /dev/null 2>&1 \
        && echo "yes [ ] $(which bedGraphToBigWig)" \
        || echo "no"
}

check_gcc () {
    command -v gcc > /dev/null 2>&1 \
        && echo "yes [ $(gcc --version 2> /dev/null | grep -E 'version|gcc') ] $(which gcc)" \
        || echo "no"
}

check_gplusplus () {
    command -v g++ > /dev/null 2>&1 \
        && echo "yes [ $(g++ --version 2> /dev/null | grep -E 'version|g\+\+') ] $(which g++)" \
        || echo "no"
}

check_python () {
    command -v python > /dev/null 2>&1 \
        && echo "yes [ $(python --version 2>&1) ] $(which python)"  \
        || echo "no"
}

check_git () {
    command -v git > /dev/null 2>&1 \
        && echo "yes [ $(git --version) ] $(which git)" \
        || echo "no"
}

check_python_package () {
    pkg="$1"
    vsn () {
        (cd $HOME && python -c "import $pkg; print ${pkg}.__version__")
    }

    pth () {
        (cd $HOME && python -c "import $pkg, os; print os.path.dirname(${pkg}.__file__)")
    }
    vsn > /dev/null 2>&1 \
        && echo "yes [ $(vsn) ] $(pth)" \
        || echo "no"
}

check_all () {

    echo "
    Compilers:
        gcc              : $(check_gcc)
        g++              : $(check_gplusplus)

    Genomics tools:
        samtools         : $(check_samtools)
        bedtools         : $(check_bedtools)
        tabix            : $(check_tabix)
        bigWigSummary    : $(check_bigWigSummary)
        bigWigToBedGraph : $(check_bigWigToBedGraph)
        bedGraphToBigWig : $(check_bedGraphToBigWig)

    Python:
        python           : $(check_python)
        scipy            : $(check_python_package scipy)
        metaseq          : $(check_python_package metaseq)
        numpy            : $(check_python_package numpy)
        matplotlib       : $(check_python_package matplotlib)
        pybedtools       : $(check_python_package pybedtools)
    "
}

if [[ "$1" == "check" ]]; then
    check_all
    exit 0
fi

# Path and log file setup ------------------------------------------

# Make the install dir and convert to absolute path
mkdir -p "${INSTALL_DIR}"
INSTALL_DIR=$(cd "${INSTALL_DIR}"; pwd)

# Paths for individual tools
BEDTOOLS_PATH="${INSTALL_DIR}/bedtools-${USE_BEDTOOLS_VERSION}"
TABIX_PATH="${INSTALL_DIR}/tabix-${USE_TABIX_VERSION}"
SAMTOOLS_PATH="${INSTALL_DIR}/samtools-${USE_SAMTOOLS_VERSION}"
BW_PATH="${INSTALL_DIR}/ucsc"

# Set up log files (or just stdout if requested)
if [[ $VERBOSE = 1 ]]; then
    BEDTOOLS_INSTALL_LOG=/dev/stdout
    SAMTOOLS_INSTALL_LOG=/dev/stdout
    TABIX_INSTALL_LOG=/dev/stdout
    BW_INSTALL_LOG=/dev/stdout
    MINICONDA_INSTALL_LOG=/dev/stdout
    ENV_INSTALL_LOG=/dev/stdout
    METASEQ_INSTALL_LOG=/dev/stdout
else
    LOG_DIR="${INSTALL_DIR}/logs"
    mkdir -p "${LOG_DIR}"
    BEDTOOLS_INSTALL_LOG="${LOG_DIR}/bedtools-installation.log"
    SAMTOOLS_INSTALL_LOG="${LOG_DIR}/samtools-installation.log"
    TABIX_INSTALL_LOG="${LOG_DIR}/tabix-installation.log"
    BW_INSTALL_LOG="${LOG_DIR}/bigWigSummary-installation.log"
    MINICONDA_INSTALL_LOG="${LOG_DIR}/miniconda-installation.log"
    ENV_INSTALL_LOG="${LOG_DIR}/miniconda-environment-installation.log"
    METASEQ_INSTALL_LOG="${LOG_DIR}/metaseq-installation.log"
fi

# Keep track of installations we found before running the script.
INITIAL=$(check_all)

cat <<EOF
Found the following prerequisites on the path:

$INITIAL

EOF

# System-specific error messages ----------------------------------------
#
#   These are too be displayed when g++/gcc are unavailable.
gplusplus_error_message () {
    if [[ $SYSTEM_TYPE = "mac" ]]; then
        cat <<EOF

    g++ is required to compile BEDTools.  On Mac, you need to install the
    correct version of XCode (https://developer.apple.com/xcode/downloads/) for
    your version of Mac OSX.

    Please install gcc and then re-run this script.

    Exiting!

EOF
    else
        cat <<EOF

    g++ is required to compile BEDTools.  On Ubuntu, please use:

        sudo apt-get install build-essential

    and then re-run this script.

    Exiting!

EOF
    fi
}

gcc_error_message () {
    if [[ $SYSTEM_TYPE = "mac" ]]; then
        cat <<EOF

    gcc is required to compile samtools and tabix.  On Mac, you need to install
    the correct version of Xcode (https://developer.apple.com/xcode/downloads/)
    for your version of Mac OSX.

    Please install gcc and then re-run this script.

    Exiting!

EOF
    else
        cat <<EOF

    gcc is required to compile samtools and tabix.  On Ubuntu, please use:

        sudo apt-get install build-essential

    and then re-run this script.

    Exiting!

EOF
    fi
}

tailinfo () {

    if [[ $VERBOSE = 1 ]]; then
        echo ""
    else
        echo "To follow progress, paste the command 'tail -f $1' without the quotes in another terminal."
    fi

}

# Installation functions -----------------------------------------------------
#
#   Each function installs the prerequisite into a subdirectory of
#   $INSTALL_DIR.  The subdirectory is named after the tool name and version
#   (see top of this script for version info).
#
#   Each function also writes the installed path $INSTALL_DIR/paths.  At the
#   end of this script, this file is sourced to show that programs have been
#   found; it also provides a single reference for paths to add to $PATH.
#
#   Functions may check for the presence of gcc or g++ as needed.
#
# ----------------------------------------------------------------------------

# Ensure we start with fresh copy of the paths file.
cat /dev/null > "${INSTALL_DIR}/paths"
echo "# Added by metaseq installation, $(date)" >> "${INSTALL_DIR}/paths"

install_bedtools () {
    # Download and install BEDTools from github releases.
    mkdir -p ${BEDTOOLS_PATH}
    ARCHIVE="${INSTALL_DIR}/bedtools-${USE_BEDTOOLS_VERSION}.tar.gz"
    (
        cd ${INSTALL_DIR} \
        && downloader \
            "https://github.com/arq5x/bedtools2/releases/download/v${USE_BEDTOOLS_VERSION}/bedtools-${USE_BEDTOOLS_VERSION}.tar.gz" "${ARCHIVE}" \
            && tar -xzf "${ARCHIVE}" --directory ${BEDTOOLS_PATH} --strip-components=1 \
            && cd "${BEDTOOLS_PATH}" \
            && make \
            && rm ${ARCHIVE} \
            && echo "export PATH=\"${BEDTOOLS_PATH}/bin:\$PATH\"" >> "${INSTALL_DIR}/paths"
    )
}

install_tabix () {
    # Download and install tabix
    if [[ check_gcc = "no" ]]; then
        gcc_error_message
        exit 1
    fi
    mkdir -p ${TABIX_PATH}
    ARCHIVE="${INSTALL_DIR}/tabix-${USE_TABIX_VERSION}.tar.bz2"
    (
        cd "${INSTALL_DIR}" \
        && downloader  "http://sourceforge.net/projects/samtools/files/tabix/tabix-${USE_TABIX_VERSION}.tar.bz2"  ${ARCHIVE} \
        && tar -xjf "${ARCHIVE}" --directory ${TABIX_PATH} --strip-components=1 \
        && cd "${TABIX_PATH}" \
        && make \
        && rm ${ARCHIVE} \
        && echo "export PATH=\"${TABIX_PATH}:\$PATH\"" >> "${INSTALL_DIR}/paths"
    )
}

install_samtools () {
    # Download and install samtools
    if [[ check_gcc = "no" ]]; then
        gcc_error_message
        exit 1
    fi
    mkdir -p ${SAMTOOLS_PATH}
    ARCHIVE="${INSTALL_DIR}/samtools-${USE_SAMTOOLS_VERSION}.tar.bz2"
    (
        downloader "https://downloads.sourceforge.net/project/samtools/samtools/${USE_SAMTOOLS_VERSION}/samtools-${USE_SAMTOOLS_VERSION}.tar.bz2" ${ARCHIVE} \
    && tar -xjf "${ARCHIVE}" --directory ${SAMTOOLS_PATH} --strip-components=1 \
    && cd "${SAMTOOLS_PATH}" \
    && make \
    && rm ${ARCHIVE} \
    && echo "export PATH=\"${SAMTOOLS_PATH}:\$PATH\"" >> "${INSTALL_DIR}/paths"
    )
}

install_ucsctools () {
    # Download and install bigWigSummary, bigWigToBedGraph, bedGraphToBigWig
    mkdir -p "${BW_PATH}" && cd "${BW_PATH}" && \
    (
        cd "${INSTALL_DIR}"
        if [[ $SYSTEM_TYPE = "linux" ]]; then
            url=http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64
        elif [[ $SYSTEM_TYPE = "mac" ]]; then
            url=http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64
        fi
        for prog in bigWigSummary bigWigToBedGraph bedGraphToBigWig; do
            downloader "$url/$prog" "${BW_PATH}/${prog}" \
            && chmod +x "${BW_PATH}/$prog" 
        done
        echo "export PATH=\"${BW_PATH}:\$PATH\"" >> "${INSTALL_DIR}/paths"

    )
}


# Genomics tools installation ---------------------------------------
# Here we're just checking for the presence of these strings in the
# $TO_INSTALL var.
#
# BEDTools installation ---------------------------------------------
echo "${TO_INSTALL}" | grep "bedtools" > /dev/null \
    && {
         log "Installing BEDTools version ${USE_BEDTOOLS_VERSION} to ${BEDTOOLS_PATH}.  $(tailinfo ${BEDTOOLS_INSTALL_LOG})" \
         && { [[ check_gplusplus = "no" ]] && { gplusplus_error_message; exit 1; } } \
         || install_bedtools > $BEDTOOLS_INSTALL_LOG 2>&1 && log "Done, see ${BEDTOOLS_INSTALL_LOG}.";
     } || log "skipping bedtools installation"

# samtools installation ---------------------------------------------
echo "${TO_INSTALL}" | grep "samtools" > /dev/null \
    && {
        log "Installing samtools version ${USE_SAMTOOLS_VERSION} to ${SAMTOOLS_PATH}. $(tailinfo ${SAMTOOLS_INSTALL_LOG})" \
        && { [[ check_gcc = "no" ]] && { gcc_error_message; exit 1; } } \
        || install_samtools > $SAMTOOLS_INSTALL_LOG 2>&1 && log "Done, see ${SAMTOOLS_INSTALL_LOG}.";
    } || log "skipping samtools installation"

# tabix installation ------------------------------------------------
echo "${TO_INSTALL}" | grep "tabix" > /dev/null \
    && {
        log "Installing tabix version ${USE_TABIX_VERSION} to ${TABIX_PATH}. $(tailinfo ${TABIX_INSTALL_LOG})" \
        && { [[ check_gcc = "no" ]] && { gcc_error_message; exit 1; } } \
        ||  install_tabix > $TABIX_INSTALL_LOG 2>&1 && log "Done, see ${TABIX_INSTALL_LOG}.";
    } || log "skipping tabix installation"

# UCSC tools installation -------------------------------------------
echo "${TO_INSTALL}" | grep "ucsc" > /dev/null \
    && {
        log "Installing UCSC tools to ${BW_PATH}. $(tailinfo ${BW_INSTALL_LOG})" \
        && install_ucsctools > $BW_INSTALL_LOG 2>&1 && log "Done, see ${BW_INSTALL_LOG}.";
    } || log "skipping UCSC tool installation"


# Miniconda installation ---------------------------------------------
if [[ "${INSTALL_MINICONDA}" == 1 ]]; then
    log "Installing Miniconda to ${MINICONDA_DIR}. $(tailinfo ${MINICONDA_INSTALL_LOG})"
    if [[ $SYSTEM_TYPE = "mac" ]]; then
        url="http://repo.continuum.io/miniconda/Miniconda-${MINICONDA_VERSION}-MacOSX-x86_64.sh"
    else
        url="http://repo.continuum.io/miniconda/Miniconda-${MINICONDA_VERSION}-Linux-x86_64.sh"
    fi

    # Only download and install if it doesn't already exist.
    if [[ ! -e ${MINICONDA_DIR} ]]; then

        install_miniconda () {
            downloader "${url}" miniconda.sh
            bash miniconda.sh -b -p "${MINICONDA_DIR}"
            hash -r
            ${MINICONDA_DIR}/bin/conda config --set always_yes yes
            ${MINICONDA_DIR}/bin/conda update conda

        }

        install_miniconda > $MINICONDA_INSTALL_LOG
        log "Done installing miniconda, see ${MINICONDA_INSTALL_LOG}."

    else
        log "${MINICONDA_DIR} exists -- using existing installation there, and environment ${ENVNAME}."
    fi

    # Set up environment
    log "Setting up isolated Python environment. $(tailinfo $ENV_INSTALL_LOG)"
    if [[ -e "${MINICONDA_DIR}/envs/${ENVNAME}" ]]; then

        # If it exists, activate then install
        log "Activating existing environment."
        source "${MINICONDA_DIR}/bin/activate" "${ENVNAME}"
        log "Installing prerequisites into existing environment. $(tailinfo ${ENV_INSTALL_LOG})' without the quotes in another terminal for details."
        conda install \
            pip ipython ipython-notebook pytables numexpr pandas biopython \
            cython matplotlib numpy scipy nose pycurl scikit-learn \
            > $ENV_INSTALL_LOG \
        && log "Done, see ${ENV_INSTALL_LOG}" \
        || { log "Error installing prerequisites, please see ${ENV_INSTALL_LOG}"; exit 1; }
    else

        # Otherwise create and then activate
        log "Installing prerequisites into a new environment. $(tailinfo ${ENV_INSTALL_LOG})"
        ${MINICONDA_DIR}/bin/conda create -n "${ENVNAME}" \
            pip ipython ipython-notebook pytables numexpr pandas biopython \
            cython matplotlib numpy scipy nose pycurl scikit-learn \
            > $ENV_INSTALL_LOG \
        && log "Done, see ${ENV_INSTALL_LOG}" \
        && log "Activating new environment." \
        && source "${MINICONDA_DIR}/bin/activate" "${ENVNAME}" \
        || { log "Error installing prerequisites, please see ${ENV_INSTALL_LOG}"; exit 1; }
    fi

    echo "# Added by metaseq installation, $(date)" > ${INSTALL_DIR}/miniconda-paths
    echo "export PATH=${MINICONDA_DIR}/bin:\$PATH" >> ${INSTALL_DIR}/miniconda-paths

    log "Done, see ${ENV_INSTALL_LOG}."
else
    log "The -M option was used, so skipping Miniconda installation."
fi


    # Metaseq installation
    log "Installing Python requirements for metaseq and metaseq itself.  $(tailinfo ${METASEQ_INSTALL_LOG})"

if [[ ${GIT_TAG} = "disable" ]]; then
    log "used -g=disable, so not installing metaseq"

elif [[ ${GIT_TAG} = "" ]]; then
    log "Installing from PyPI, $(tailinfo ${METASEQ_INSTALL_LOG})"
    pip install "metaseq" > $METASEQ_INSTALL_LOG \
    && log "Done, see ${METASEQ_INSTALL_LOG}" \
    || { log "Error installing metaseq from PyPI, see ${METASEQ_INSTALL_LOG}"; exit 1; }

else
    if [[ $(check_git) = "no" ]]; then
        echo "
        git does not appear to be installed, so the metaseq github repository
        cannot be cloned.

        Exiting!
        "
        exit 1
    fi
    log "Cloning metaseq repository"
    git clone https://github.com/daler/metaseq.git ${INSTALL_DIR}/metaseq \
    && log "Checking out ${GIT_TAG} and installing" \
    && ( cd ${INSTALL_DIR}/metaseq && git checkout $GIT_TAG && pip install . > ${METASEQ_INSTALL_LOG} ) \
    && log "Done, see ${METASEQ_INSTALL_LOG}" \
    || { log "Error installing metaseq from git clone, see ${METASEQ_INSTALL_LOG}"; exit 1; }
fi


# Each installation process wrote the installed path to this file.  Source it
# now, so we can confirm that the tools were found.  Note that at this point
# the miniconda env should still be activated.
source "${INSTALL_DIR}/paths"

README="${INSTALL_DIR}/README.txt"
cat /dev/null > ${README}

cat >> ${README} <<EOF
Results
-------
Before running this script, the following programs were detected on your path:

$INITIAL


After running this script the following programs can now be detected:

$(check_all)
EOF


[[ $SYSTEM_TYPE = "mac" ]] && profile="${HOME}/.bash_profile" || profile="${HOME}/.bashrc"


if [[ $(grep "export" ${INSTALL_DIR}/paths) ]]; then
    # If on travis-ci, we definitely don't want interactive prompts holding
    # things up....
    if [[ $TRAVIS_CI != 1 ]]; then
        cat << EOF

    I can automatically add the locations of the installed programs to the
    beginning of your \$PATH, which will allow these programs to be
    found when you open a new terminal. This is a good idea if you're not
    comfortable with editing your \$PATH.  Specifically, I can add the contents
    of ${INSTALL_DIR}/paths to the end of your ${profile} file.
    You can always delete these lines later.

    Do you want me to add these programs to your path now?
EOF
        echo -n "(yes/no): "
        read ans

        if [[ ($ans != "yes") && ($ans != "Yes") && ($ans != "YES") &&
            ($ans != "y") && ($ans != "Y") ]]; then

            # If not added automatically here,, then write instructions to the
            # README
            cat >> ${README} <<EOF
Add installed programs to PATH
------------------------------
The paths to these installed programs have been written to the file
${INSTALL_DIR}/paths.  In order to permanently install them, you will need to
add these directories to your ${profile} file.  Pasting this line in a terminal
should do the trick:

    cat ${INSTALL_DIR}/paths >> ${profile}

Then open a new terminal window to activate the changes.

EOF
        else
            # Otherwise do it and then log the changes to the README.
            cat ${INSTALL_DIR}/paths >> ${profile}
            cat >> ${README} << EOF
[ $(date) ]: The following lines were appended to the end of ${profile}:

$(cat ${INSTALL_DIR}/paths)

EOF
            log "Added paths to ${profile}"
        fi # end if y/n
    fi # end if TRAVIS_CI
fi # end if something was installed

if [[ ${INSTALL_MINICONDA} = 1 ]]; then
    if [[ $TRAVIS_CI != 1 ]]; then
        cat << EOF

    Would you like to prepend the location of the Miniconda installation to
    your \$PATH? This will make it easier to activate enivronments.  This is
    a good idea  if you're not comfortable with editing your \$PATH.
    Specifically, I can add the contents of
    ${INSTALL_DIR}/miniconda-path to the end of your ${profile} file.
    You can always delete these lines later.

    Do you want to add the Miniconda installation to your path now?
EOF
        echo -n "(yes/no): "
        read ans

        if [[ ($ans != "yes") && ($ans != "Yes") && ($ans != "YES") &&
            ($ans != "y") && ($ans != "Y") ]]; then

            # If not added automatically here,, then write instructions to the
            # README
            cat >> ${README} <<EOF


Using the miniconda installation
--------------------------------
A miniconda installation is now at ${MINICONDA_DIR}. To make this permanently
available, please paste the following line in a terminal, which will prepend
the miniconda installation to your PATH:

    cat ${INSTALL_DIR}/miniconda-paths >> ${profile}

Then open a new terminal to activate the changes.
EOF
        else
            # Otherwise do it and then log the changes to the README.
            cat ${INSTALL_DIR}/miniconda-paths >> ${profile}
            cat >> ${README} << EOF
[ $(date) ]: The following lines were appended to the end of ${profile}:

$(cat ${INSTALL_DIR}/miniconda-paths)

EOF
            log "Added miniconda path to ${profile}"
        fi # end if y/n

    cat >> ${README} <<EOF

    Using the new environment
    -------------------------
    An isolated Python environment has been created called ${ENVNAME}.

    When you are ready to use it, open a new terminal and use the command:

        source activate ${ENVNAME}

    Now you will be using the isolated Python environment; anything you do will not
    touch the system-wide installation. When you're done, use

        source deactivate

    to return to normal.

    Planning on trying out metaseq?  After activating the environment you can
    download the example data with:

        download_metaseq_example_data.py

    and follow the tutorial at https://pythonhosted.org/metaseq/example_session.html

EOF

    fi # end if travis-ci
fi # end if miniconda installed

log "

    Done!  The results of this installation have been written to:

        ${README}

    Please follow the instructions in that file to complete the installation.
"
