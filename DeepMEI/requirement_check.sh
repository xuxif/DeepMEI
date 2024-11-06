#!/bin/bash

# Function to check if a command exists
check_command() {
    command -v $1 >/dev/null 2>&1
}

# Function to check version within a range or exact match
check_version() {
    installed_version=$($1 2>&1 | grep -oP "$2" | head -1)
    
    if [[ -z $installed_version ]]; then
        echo "Could not determine the version for $1."
        return
    fi

    lower_bound=$3
    upper_bound=$4
    if [[ -z $upper_bound ]]; then
        comparison=$(python3 -c "from packaging import version; print(version.parse('$installed_version') >= version.parse('$lower_bound'))")
    else
        comparison=$(python3 -c "from packaging import version; print(version.parse('$installed_version') >= version.parse('$lower_bound') and version.parse('$installed_version') < version.parse('$upper_bound'))")
    fi
    if [ "$comparison" == "True" ]; then
        echo "$1 version $installed_version meets the requirement ($lower_bound to ${upper_bound:-})."
    else
        echo "$1 version $installed_version does not meet the requirement ($lower_bound to ${upper_bound:-})."
    fi
}

# Check samtools
if check_command "samtools"; then
    check_version "samtools" "[0-9]+\.[0-9]+" "1.15"
else
    echo "samtools is not installed."
fi

# Check samtools import and fastq
if check_command "samtools"; then
    if samtools 2>&1 | grep -q "import"; then
        echo "samtools import is available."
    else
        echo "samtools import is not available."
    fi
    
    if samtools 2>&1 | grep -q "fastq"; then
        echo "samtools fastq is available."
    else
        echo "samtools fastq is not available."
    fi
else
    echo "samtools is not installed."
fi

# Check BWA
if check_command "bwa"; then
    check_version "bwa" "[0-9]+\.[0-9]+\.[0-9]+" "0.7.17"
else
    echo "bwa is not installed."
fi

# Check bedtools
if check_command "bedtools"; then
    check_version "bedtools" "[0-9]+\.[0-9]+\.[0-9]+" "2.30"
else
    echo "bedtools is not installed."
fi

# Check RepeatMasker
if check_command "RepeatMasker"; then
    check_version "RepeatMasker" "[0-9]+\.[0-9]+\.[0-9]+" "4.1.2"
else
    echo "RepeatMasker is not installed."
fi

# Check Python version
python_version=$(python3 --version 2>&1 | grep -oP "[0-9]+\.[0-9]+\.[0-9]+")
python_required="3.8"
python_max="3.10"
if [[ -n $python_version ]]; then
    if python3 -c "import sys; from packaging import version; sys.exit(not (version.parse('$python_version') >= version.parse('$python_required') and version.parse('$python_version') < version.parse('$python_max')))" ; then
        echo "Python version $python_version meets the requirement ($python_required to $python_max)."
    else
        echo "Python version $python_version does not meet the requirement ($python_required to $python_max)."
    fi
else
    echo "Could not determine the Python version."
fi

# Check Perl version
perl_version=$(perl -v | grep -oP "v[0-9]+\.[0-9]+\.[0-9]+" | tr -d 'v')
perl_required="5.32"
if [[ -n $perl_version ]]; then
    if perl -e "use version; exit(version->parse('$perl_version') >= version->parse('$perl_required') ? 0 : 1)"; then
        echo "Perl version $perl_version meets the requirement (>= $perl_required)."
    else
        echo "Perl version $perl_version does not meet the requirement (>= $perl_required)."
    fi
else
    echo "Could not determine the Perl version."
fi

# Check for Python packages using Python to directly get the version
check_python_package_version() {
    package_version=$(python3 -c "import $1; print($1.__version__)" 2>/dev/null)
    
    if [[ -z $package_version ]]; then
        echo "Could not determine the version for Python package $1."
        return
    fi

    required_version=$2
    upper_bound=$3
    if [[ -z $upper_bound ]]; then
        comparison=$(python3 -c "from packaging import version; print(version.parse('$package_version') == version.parse('$required_version'))")
    else
        comparison=$(python3 -c "from packaging import version; print(version.parse('$package_version') >= version.parse('$required_version') and version.parse('$package_version') < version.parse('$upper_bound'))")
    fi
    if [ "$comparison" == "True" ]; then
        echo "Python package $1 version $package_version meets the requirement ($required_version to ${upper_bound:-})."
    else
        echo "Python package $1 version $package_version does not meet the requirement ($required_version to ${upper_bound:-})."
    fi
}

# Check specific Python packages
check_python_package_version "pysam" "0.17.0"
check_python_package_version "tensorflow" "2.7.0" "2.12.0"

